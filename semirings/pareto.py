"""Pareto-frontier semiring (min convention, arbitrary dimension).

Elements are antichains of `Point`s under componentwise dominance: a point p
dominates q iff p_i <= q_i for all i with at least one strict inequality.
The semiring keeps only non-dominated points; each surviving `Point` carries
the derivation backpointer assembled by `Point.__add__` during Minkowski
sums, so the antichain doubles as an argmin record.

  +    : union, pruned to Pareto-minimal
  *    : Minkowski sum (componentwise add), pruned to Pareto-minimal
  zero : empty antichain
  one  : {(0, ..., 0)}
  star : `one` whenever every point has all coordinates >= 0 (the origin
         dominates every Minkowski iterate); otherwise a bounded fixpoint
         iteration that may fail to converge for genuinely open elements
         (e.g., points with mixed-sign coordinates).

Use `make_pareto(d)` to build a class for a specific dimension. The default
`Pareto = make_pareto(2)` is exposed for the common 2D case.
"""

from semirings.base import Semiring
from semirings.convex_hull import Point, _directed_hausdorff


_STAR_MAX_ITER = 64


def _lift(p):
    return p if isinstance(p, Point) else Point(tuple(p))


def _prune_min(points):
    """Return the Pareto-minimal subset of `points`, sorted canonically.

    Dedup is keyed on coords (Point's __eq__/__hash__ already are), so
    distinct derivations producing the same coords collapse to whichever
    arrived first.
    """
    seen = {}
    for p in points:
        seen.setdefault(p.coords, p)
    pts = list(seen.values())
    if not pts:
        return ()
    keep = [p for p in pts
            if not any(_dominates(q, p) for q in pts if q.coords != p.coords)]
    return tuple(sorted(keep, key=lambda p: p.coords))


def _dominates(p, q):
    """True iff `p` strictly Pareto-dominates `q` (min convention)."""
    le = True
    strict = False
    for a, b in zip(p.coords, q.coords):
        if a > b:
            le = False
            break
        if a < b:
            strict = True
    return le and strict


def make_pareto(d, name=None):
    """Build the d-dimensional Pareto-frontier semiring class."""

    class _Pareto(Semiring):
        dim = d

        def __init__(self, points=()):
            pts = [_lift(p) for p in points]
            for p in pts:
                if p.dim != d:
                    raise ValueError(
                        f'expected {d}-dim points, got {p.dim}-dim {p.coords}'
                    )
            self.points = _prune_min(pts)

        def __repr__(self):
            return f'{type(self).__name__}({list(self.points)})'

        def __iter__(self):
            return iter(self.points)

        def __len__(self):
            return len(self.points)

        def __eq__(self, other):
            return isinstance(other, _Pareto) and self.points == other.points

        def __hash__(self):
            return hash(self.points)

        def __add__(self, other):
            if not isinstance(other, _Pareto):
                return NotImplemented
            return _Pareto(self.points + other.points)

        def __mul__(self, other):
            if not isinstance(other, _Pareto):
                return NotImplemented
            if not self.points or not other.points:
                return _Pareto.zero
            return _Pareto([p + q for p in self.points for q in other.points])

        def star(self):
            # Fast path: if every point has all non-negative coords, the
            # origin (= one) dominates every Minkowski iterate.
            if all(c >= 0 for p in self.points for c in p.coords):
                return _Pareto.one
            # Bounded fixpoint iteration for the general case.
            prev = _Pareto.zero
            for _ in range(_STAR_MAX_ITER):
                curr = _Pareto.one + self * prev
                if curr == prev:
                    return curr
                prev = curr
            raise ValueError(
                f'{type(self).__name__}.star did not converge in '
                f'{_STAR_MAX_ITER} iterations on {self!r}; this element is '
                f'in the open part of the carrier (points with mixed signs '
                f'producing iterates that drift outward indefinitely).'
            )

        def metric(self, other):
            """Hausdorff distance with Chebyshev (L_inf) ground metric.

            Returns 0 for equal antichains, +inf when exactly one side is
            empty (otherwise identity-of-indiscernibles would fail).
            """
            if not isinstance(other, _Pareto):
                return float('inf')
            if self.points == other.points:
                return 0.0
            if not self.points or not other.points:
                return float('inf')
            return max(_directed_hausdorff(self.points, other.points),
                       _directed_hausdorff(other.points, self.points))

        def plot(self, title=None, color='#1f77b4', marker_size=10,
                 show_origin=True, labels=None):
            """Plotly figure for d in {1, 2, 3}; raises otherwise.

            Hover tooltips show the full coordinate tuple plus an optional
            user-supplied label.
            """
            return _plot(self, title=title, color=color,
                         marker_size=marker_size, show_origin=show_origin,
                         labels=labels)

    _Pareto.__name__ = name or f'Pareto{d}D'
    _Pareto.__qualname__ = _Pareto.__name__
    _Pareto.zero = _Pareto([])
    _Pareto.one = _Pareto([(0,) * d])
    return _Pareto


def _plot(P, title, color, marker_size, show_origin, labels):
    try:
        import plotly.graph_objects as go
    except ImportError as e:
        raise ImportError(
            'plot() requires plotly: pip install plotly'
        ) from e

    pts = list(P.points)
    if not pts:
        fig = go.Figure()
        fig.update_layout(title=title or f'{type(P).__name__}: empty')
        return fig

    if labels is None:
        labels = [f'p{i}' for i in range(len(pts))]
    elif len(labels) != len(pts):
        raise ValueError(f'expected {len(pts)} labels, got {len(labels)}')

    d = P.dim
    if d == 1:
        return _plot_1d(go, pts, labels, title, color, marker_size, show_origin, P)
    elif d == 2:
        return _plot_2d(go, pts, labels, title, color, marker_size, show_origin, P)
    elif d == 3:
        return _plot_3d(go, pts, labels, title, color, marker_size, show_origin, P)
    raise NotImplementedError(f'plot() supports d in {{1,2,3}}; got d={d}')


def _hover(label, point):
    coords = ', '.join(f'{x:g}' for x in point.coords)
    return f'{label}<br>({coords})'


def _plot_1d(go, pts, labels, title, color, marker_size, show_origin, P):
    xs = [p.coords[0] for p in pts]
    fig = go.Figure(go.Scatter(
        x=xs, y=[0] * len(xs), mode='markers+text',
        marker=dict(size=marker_size, color=color, line=dict(width=1, color='white')),
        text=labels, textposition='top center',
        hovertext=[_hover(l, p) for l, p in zip(labels, pts)],
        hoverinfo='text', name='frontier',
    ))
    if show_origin:
        fig.add_vline(x=0, line=dict(color='gray', width=0.6, dash='dot'))
    fig.update_layout(
        title=title or f'{type(P).__name__} ({len(pts)} pts)',
        xaxis_title='cost_1', yaxis=dict(visible=False, range=[-1, 1]),
        height=200, showlegend=False,
    )
    return fig


def _plot_2d(go, pts, labels, title, color, marker_size, show_origin, P):
    pts_sorted = sorted(pts, key=lambda p: p.coords)
    xs = [p.coords[0] for p in pts_sorted]
    ys = [p.coords[1] for p in pts_sorted]
    label_map = {p: l for l, p in zip(labels, pts)}
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='lines',
        line=dict(color=color, width=1, dash='dot'),
        hoverinfo='skip', showlegend=False,
    ))
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='markers+text',
        marker=dict(size=marker_size, color=color, line=dict(width=1, color='white')),
        text=[label_map[p] for p in pts_sorted],
        textposition='top right',
        hovertext=[_hover(label_map[p], p) for p in pts_sorted],
        hoverinfo='text', name='frontier',
    ))
    if show_origin:
        fig.add_hline(y=0, line=dict(color='gray', width=0.5, dash='dot'))
        fig.add_vline(x=0, line=dict(color='gray', width=0.5, dash='dot'))
    fig.update_layout(
        title=title or f'{type(P).__name__} ({len(pts)} pts)',
        xaxis_title='cost_1', yaxis_title='cost_2',
        height=420, width=480, showlegend=False,
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


def _plot_3d(go, pts, labels, title, color, marker_size, show_origin, P):
    xs = [p.coords[0] for p in pts]
    ys = [p.coords[1] for p in pts]
    zs = [p.coords[2] for p in pts]
    fig = go.Figure(go.Scatter3d(
        x=xs, y=ys, z=zs, mode='markers+text',
        marker=dict(size=marker_size, color=color, line=dict(width=1, color='white')),
        text=labels,
        hovertext=[_hover(l, p) for l, p in zip(labels, pts)],
        hoverinfo='text', name='frontier',
    ))
    fig.update_layout(
        title=title or f'{type(P).__name__} ({len(pts)} pts)',
        scene=dict(xaxis_title='cost_1', yaxis_title='cost_2', zaxis_title='cost_3'),
        height=520, width=560, showlegend=False,
    )
    return fig


Pareto = make_pareto(2, 'Pareto')
