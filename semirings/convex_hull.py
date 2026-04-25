"""d-dimensional convex-hull semiring.

Generalizes Dyer's 2D MERT semiring (arxiv:1307.3675) to arbitrary dimension.
For d=2 both `+` (hull of union) and `*` (Minkowski sum) run in O(n+m) on
already-convex inputs via chain merging and edge-angle merging respectively.
For d∉{1,2} we fall back to qhull (scipy.spatial.ConvexHull).

Vertices for d=2 are stored CCW from the leftmost-lowest point; the linear-time
fast paths rely on this canonical orientation.
"""
import numpy as np
import scipy.spatial

from semirings.base import Semiring


class Point:
    """Point in R^d with optional derivation backpointer."""

    __slots__ = ('coords', 'd')

    def __init__(self, *coords, d=None):
        if len(coords) == 1 and isinstance(coords[0], (tuple, list, np.ndarray)):
            coords = tuple(coords[0])
        self.coords = tuple(coords)
        self.d = d

    @property
    def x(self): return self.coords[0]
    @property
    def y(self): return self.coords[1]
    @property
    def z(self): return self.coords[2]
    @property
    def dim(self): return len(self.coords)

    def __add__(self, other):
        return Point(
            *(a + b for a, b in zip(self.coords, other.coords)),
            d=(self, other),
        )

    def __eq__(self, other):
        return isinstance(other, Point) and self.coords == other.coords

    def __hash__(self):
        return hash(self.coords)

    def __repr__(self):
        return str(self.d) if self.d is not None else f'Point{self.coords}'

    def metric(self, other):
        """Chebyshev (L_inf) distance: max coordinate-wise |a - b|."""
        return max(abs(a - b) for a, b in zip(self.coords, other.coords))


class ConvexHull(Semiring):
    """Convex-hull semiring in R^d.

    Each element is a convex polytope, stored as its vertex set. `+` is hull
    of union (idempotent); `*` is Minkowski sum. `star()` diverges and is
    not defined.
    """

    __slots__ = ('points', 'dim')

    def __init__(self, points, *, _vetted=False):
        points = tuple(points)
        if _vetted:
            self.points = points
            self.dim = points[0].dim if points else 0
            return
        if not points:
            self.points = ()
            self.dim = 0
            return
        self.dim = points[0].dim
        self.points = tuple(_hull(points, self.dim))

    def __iter__(self): return iter(self.points)
    def __len__(self): return len(self.points)

    def __eq__(self, other):
        return isinstance(other, ConvexHull) and set(self.points) == set(other.points)

    def __hash__(self):
        return hash(frozenset(self.points))

    def __repr__(self):
        return f'{type(self).__name__}({list(self.points)})'

    def __add__(self, other):
        if not self.points: return other
        if not other.points: return self
        assert isinstance(other, ConvexHull)
        assert self.dim == other.dim, f"dim mismatch: {self.dim} vs {other.dim}"
        if self.dim == 2:
            return ConvexHull(_union_2d(self.points, other.points), _vetted=True)
        return ConvexHull(self.points + other.points)

    def __mul__(self, other):
        # Origin singletons act as multiplicative identity in any dim.
        if _is_origin_singleton(other): return self
        if _is_origin_singleton(self): return other
        if not self.points or not other.points: return _zero
        assert isinstance(other, ConvexHull)
        assert self.dim == other.dim, f"dim mismatch: {self.dim} vs {other.dim}"
        if self.dim == 2 and len(self.points) >= 3 and len(other.points) >= 3:
            return ConvexHull(_minkowski_2d(self.points, other.points), _vetted=True)
        if len(self.points) == 1:
            s = self.points[0]
            return ConvexHull([s + p for p in other.points], _vetted=True)
        if len(other.points) == 1:
            o = other.points[0]
            return ConvexHull([p + o for p in self.points], _vetted=True)
        return ConvexHull([a + b for a in self.points for b in other.points])

    @classmethod
    def multiplicity(cls, v, m):
        return v if m > 0 else cls.zero

    def star(self):
        raise NotImplementedError("star diverges for nontrivial hulls")

    def metric(self, other):
        """Hausdorff distance with Chebyshev (L_inf) ground metric.

        Returns 0 for equal hulls, +inf when exactly one side is empty
        (otherwise identity-of-indiscernibles would fail). For convex
        polytopes the max distance to the other set is attained at a vertex
        of self, so iterating over `self.points`/`other.points` is exact.
        """
        if not isinstance(other, ConvexHull):
            return float('inf')
        if self.points == other.points:
            return 0.0
        if not self.points or not other.points:
            return float('inf')
        return max(_directed_hausdorff(self.points, other.points),
                   _directed_hausdorff(other.points, self.points))

    def draw(self, label_fmt=None, c='k', ax=None):
        return _draw(self, label_fmt, c, ax)

    def plot(self, title=None, color='#1f77b4', marker_size=8, opacity=0.18,
             labels=None):
        """Plotly figure for d in {1, 2, 3}; raises otherwise.

        Hover tooltips show the full coordinate tuple plus an optional
        user-supplied label (one per vertex, in `self.points` order).
        """
        return _plot(self, title=title, color=color, marker_size=marker_size,
                     opacity=opacity, labels=labels)


def _is_origin_singleton(h):
    return len(h.points) == 1 and not any(h.points[0].coords)


# --- hull computation ---------------------------------------------------------

def _hull(points, d):
    if len(points) <= 1:
        return list(points)
    if d == 1:
        lo = min(points, key=lambda p: p.coords[0])
        hi = max(points, key=lambda p: p.coords[0])
        return [lo] if lo.coords == hi.coords else [lo, hi]
    if d == 2:
        seen = {}
        for p in points:
            seen.setdefault(p.coords, p)
        return _monotone_chain(sorted(seen.values(),
                                      key=lambda p: (p.coords[0], p.coords[1])))
    arr = np.array([p.coords for p in points], dtype=float)
    try:
        # No QJ here: joggling perturbs inputs and can promote interior pair-sums
        # to spurious "vertices", breaking distributivity for the semiring.
        c = scipy.spatial.ConvexHull(arr)
        return [points[i] for i in c.vertices]
    except scipy.spatial.QhullError:
        # Input lies in a lower-dim affine subspace — project and recurse.
        return _hull_lower_dim(points, d)


def _hull_lower_dim(points, d):
    """Hull of points that lie in a lower-dim affine subspace of R^d.

    SVD the centered coordinates, project to the active subspace, and run
    qhull there. We carry indices through the projection so the returned
    Points are the originals (with their backpointers intact).
    """
    seen = {}
    for p in points:
        seen.setdefault(p.coords, p)
    uniq = list(seen.values())
    if len(uniq) <= d:
        return uniq
    arr = np.array([p.coords for p in uniq], dtype=float)
    centered = arr - arr.mean(axis=0)
    _, s, vh = np.linalg.svd(centered, full_matrices=False)
    rank = int((s > 1e-9 * (s.max() if s.size else 1)).sum())
    if rank >= d or rank == 0:
        return uniq
    proj = centered @ vh[:rank].T
    if rank == 1:
        lo = int(proj[:, 0].argmin())
        hi = int(proj[:, 0].argmax())
        return [uniq[lo]] if lo == hi else [uniq[lo], uniq[hi]]
    try:
        ch = scipy.spatial.ConvexHull(proj)
        return [uniq[i] for i in ch.vertices]
    except scipy.spatial.QhullError:
        return uniq


def _monotone_chain(pts):
    """Andrew's monotone chain on (x,y)-presorted points.

    Returns CCW vertices starting from the leftmost-lowest point. Strips
    collinear points (strict turn test).
    """
    if len(pts) <= 1:
        return list(pts)
    cross = lambda o, a, b: (
        (a.coords[0] - o.coords[0]) * (b.coords[1] - o.coords[1])
        - (a.coords[1] - o.coords[1]) * (b.coords[0] - o.coords[0])
    )
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


# --- 2D fast paths ------------------------------------------------------------
#
# Both rely on the polygon being in canonical form: CCW, starting at the
# leftmost-lowest vertex. With this layout the polygon's x-coords are bitonic
# (ascending to the rightmost vertex, then descending), which lets us split it
# into two presorted-by-x chains in linear time.

def _sorted_xy_from_polygon(P):
    """Vertices of a canonical CCW polygon, sorted by (x, y). O(n)."""
    if len(P) <= 2:
        return sorted(P, key=lambda p: (p.coords[0], p.coords[1]))
    n = len(P)
    # Largest x with tiebreak largest y — the end of the lower chain when a
    # vertical edge sits at the right side of the hull.
    k = max(range(n), key=lambda i: (P[i].coords[0], P[i].coords[1]))
    lower = list(P[:k + 1])                       # x ascending
    upper = list(reversed(P[k:] + (P[0],)))       # x ascending
    return _merge_xy(lower, upper)


def _merge_xy(A, B):
    """Merge two lists pre-sorted by (x, y), deduping coord-equal points."""
    out = []
    i, j, nA, nB = 0, 0, len(A), len(B)
    while i < nA and j < nB:
        ka, kb = A[i].coords, B[j].coords
        if ka < kb:
            out.append(A[i]); i += 1
        elif ka > kb:
            out.append(B[j]); j += 1
        else:
            out.append(A[i]); i += 1; j += 1
    out.extend(A[i:]); out.extend(B[j:])
    return out


def _union_2d(P, Q):
    """Hull of P ∪ Q for canonical CCW polygons P, Q. O(n+m)."""
    return _monotone_chain(_merge_xy(_sorted_xy_from_polygon(P),
                                     _sorted_xy_from_polygon(Q)))


def _minkowski_2d(P, Q):
    """Minkowski sum of two CCW polygons (≥3 vertices each).

    Berg et al. edge-angle merge: O(n+m). The merge wants each polygon to
    start at its lowest-leftmost vertex; canonical form is leftmost-lowest,
    so we rotate at both ends.
    """
    P = _rotate_to(list(P), key=lambda p: (p.coords[1], p.coords[0]))
    Q = _rotate_to(list(Q), key=lambda p: (p.coords[1], p.coords[0]))
    n, m = len(P), len(Q)
    R = []
    i, j = 0, 0
    while i < n or j < m:
        R.append(P[i % n] + Q[j % m])
        if i == n:
            j += 1
        elif j == m:
            i += 1
        else:
            ep0 = P[(i + 1) % n].coords[0] - P[i].coords[0]
            ep1 = P[(i + 1) % n].coords[1] - P[i].coords[1]
            eq0 = Q[(j + 1) % m].coords[0] - Q[j].coords[0]
            eq1 = Q[(j + 1) % m].coords[1] - Q[j].coords[1]
            cross = ep0 * eq1 - ep1 * eq0
            if cross > 0:   i += 1
            elif cross < 0: j += 1
            else:           i += 1; j += 1
    return _rotate_to(R, key=lambda p: (p.coords[0], p.coords[1]))


def _rotate_to(P, *, key):
    k = min(range(len(P)), key=lambda i: key(P[i]))
    return P[k:] + P[:k]


# --- visualization ------------------------------------------------------------

def _draw(hull, label_fmt, c, ax):  # pragma: no cover
    import pylab as pl
    if hull.dim not in (1, 2, 3):
        raise ValueError(f"draw only supports dim ∈ {{1,2,3}} (got {hull.dim})")
    if not hull.points:
        return
    label = label_fmt if label_fmt is not None else (lambda p: str(p.d))
    coords = np.array([p.coords for p in hull.points])

    if hull.dim == 1:
        if ax is None:
            ax = pl.gca()
        xs = coords[:, 0]
        if len(xs) >= 2:
            ax.plot([xs.min(), xs.max()], [0, 0],
                    c=c, alpha=0.4, lw=1, zorder=1)
        ax.scatter(xs, np.zeros(len(xs)), c=c, alpha=0.7, zorder=2)
        for p in hull.points:
            ax.text(p.coords[0], 0, label(p))
        ax.set_yticks([])

    elif hull.dim == 2:
        if ax is None:
            ax = pl.gca()
        if len(coords) >= 3:
            sc = scipy.spatial.ConvexHull(coords, qhull_options='QJ Pp')
            for s in sc.simplices:
                ax.plot(coords[s, 0], coords[s, 1],
                        c=c, alpha=0.5, lw=0.5, zorder=1)
        elif len(coords) == 2:
            ax.plot(coords[:, 0], coords[:, 1], c=c, alpha=0.5, lw=0.5)
        ax.scatter(coords[:, 0], coords[:, 1], c=c, alpha=0.7, zorder=2)
        for p in hull.points:
            ax.text(p.coords[0], p.coords[1], label(p))
        ax.set_aspect('equal')

    else:  # 3D
        if ax is None:
            ax = pl.gcf().add_subplot(111, projection='3d')
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=c, alpha=0.7)
        if len(coords) >= 4:
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            sc = scipy.spatial.ConvexHull(coords, qhull_options='QJ Pp')
            faces = [coords[s] for s in sc.simplices]
            ax.add_collection3d(Poly3DCollection(
                faces, alpha=0.2, edgecolor=c, facecolor=c, linewidth=0.5))
        for p in hull.points:
            ax.text(p.coords[0], p.coords[1], p.coords[2], label(p))


# --- metric -------------------------------------------------------------------

def _directed_hausdorff(A, B):
    """sup_{a in A} inf_{b in B} a.metric(b) over `Point` iterables."""
    return max(min(a.metric(b) for b in B) for a in A)


# --- plotly visualization -----------------------------------------------------

def _plot(hull, title, color, marker_size, opacity, labels):
    try:
        import plotly.graph_objects as go
    except ImportError as e:
        raise ImportError('plot() requires plotly: pip install plotly') from e

    pts = list(hull.points)
    if not pts:
        fig = go.Figure()
        fig.update_layout(title=title or f'{type(hull).__name__}: empty')
        return fig

    if labels is None:
        labels = [f'p{i}' for i in range(len(pts))]
    elif len(labels) != len(pts):
        raise ValueError(f'expected {len(pts)} labels, got {len(labels)}')

    d = hull.dim
    if d == 1: return _plot_1d(go, pts, labels, title, color, marker_size, hull)
    if d == 2: return _plot_2d(go, pts, labels, title, color, marker_size, opacity, hull)
    if d == 3: return _plot_3d(go, pts, labels, title, color, marker_size, opacity, hull)
    raise NotImplementedError(f'plot() supports d in {{1,2,3}}; got d={d}')


def _hover(label, point):
    return f'{label}<br>({", ".join(f"{x:g}" for x in point.coords)})'


def _plot_1d(go, pts, labels, title, color, marker_size, hull):
    xs = [p.coords[0] for p in pts]
    fig = go.Figure()
    if len(xs) >= 2:
        fig.add_trace(go.Scatter(
            x=[min(xs), max(xs)], y=[0, 0], mode='lines',
            line=dict(color=color, width=1), hoverinfo='skip', showlegend=False,
        ))
    fig.add_trace(go.Scatter(
        x=xs, y=[0] * len(xs), mode='markers+text',
        marker=dict(size=marker_size, color=color, line=dict(width=1, color='white')),
        text=labels, textposition='top center',
        hovertext=[_hover(l, p) for l, p in zip(labels, pts)],
        hoverinfo='text', name='hull',
    ))
    fig.update_layout(
        title=title or f'{type(hull).__name__} ({len(pts)} pts)',
        xaxis_title='x', yaxis=dict(visible=False, range=[-1, 1]),
        height=200, showlegend=False,
    )
    return fig


def _plot_2d(go, pts, labels, title, color, marker_size, opacity, hull):
    # `pts` is canonical CCW; close the polygon for the fill trace.
    xs = [p.coords[0] for p in pts]
    ys = [p.coords[1] for p in pts]
    fig = go.Figure()
    if len(pts) >= 3:
        fig.add_trace(go.Scatter(
            x=xs + [xs[0]], y=ys + [ys[0]], mode='lines', fill='toself',
            line=dict(color=color, width=1), fillcolor=color, opacity=opacity,
            hoverinfo='skip', showlegend=False,
        ))
    elif len(pts) == 2:
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode='lines',
            line=dict(color=color, width=1), hoverinfo='skip', showlegend=False,
        ))
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='markers+text',
        marker=dict(size=marker_size, color=color, line=dict(width=1, color='white')),
        text=labels, textposition='top right',
        hovertext=[_hover(l, p) for l, p in zip(labels, pts)],
        hoverinfo='text', name='hull',
    ))
    fig.update_layout(
        title=title or f'{type(hull).__name__} ({len(pts)} pts)',
        xaxis_title='x', yaxis_title='y',
        height=420, width=480, showlegend=False,
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


def _plot_3d(go, pts, labels, title, color, marker_size, opacity, hull):
    coords = np.array([p.coords for p in pts], dtype=float)
    xs, ys, zs = coords[:, 0], coords[:, 1], coords[:, 2]
    fig = go.Figure()
    if len(pts) >= 4:
        try:
            sc = scipy.spatial.ConvexHull(coords)
            i, j, k = sc.simplices.T
            fig.add_trace(go.Mesh3d(
                x=xs, y=ys, z=zs, i=i, j=j, k=k,
                color=color, opacity=opacity, flatshading=True,
                hoverinfo='skip', showlegend=False,
            ))
        except scipy.spatial.QhullError:
            pass
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode='markers+text',
        marker=dict(size=max(3, marker_size // 2), color=color,
                    line=dict(width=1, color='white')),
        text=labels,
        hovertext=[_hover(l, p) for l, p in zip(labels, pts)],
        hoverinfo='text', name='hull',
    ))
    fig.update_layout(
        title=title or f'{type(hull).__name__} ({len(pts)} pts)',
        scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='z'),
        height=520, width=560, showlegend=False,
    )
    return fig


_zero = ConvexHull((), _vetted=True)
_one = ConvexHull((Point(0, 0),), _vetted=True)

ConvexHull.zero = _zero
ConvexHull.one = _one
