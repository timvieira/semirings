"""MERT semiring — the 2D convex-hull semiring (Dyer, 2013)."""
import numpy as np
import pylab as pl
import scipy.spatial

from semirings.base import Semiring


# TODO: there is nothing specific about the convexhull semiring and being
# two-dimensional.
#
# TODO: the union of convex hulls can be computed in linear time (no log factor).
#
# TODO: convex hull commutes with Minkowski sum (no second call to conv needed).


class ConvexHull(Semiring):
    """Convex-hull semiring (two-dimensional).

    Each element is a convex hull of Points. Based on Chris Dyer's 2013 arxiv
    paper (http://arxiv.org/pdf/1307.3675.pdf). Coded for clarity, not efficiency.
    """

    def __init__(self, points):
        self.points = conv(points)

    def __iter__(self):
        return iter(self.points)

    def __eq__(self, other):
        return set(self.points) == set(other.points)

    def __hash__(self):
        return hash(frozenset(self.points))

    def __repr__(self):
        return f'{self.__class__.__name__}({self.points})'

    def __add__(self, other):
        if self is _zero: return other
        if other is _zero: return self
        assert isinstance(other, ConvexHull)
        return ConvexHull(list(self) + list(other))

    def __mul__(self, other):
        if other is _one: return self
        if self is _one: return other
        if other is _zero: return _zero
        if self is _zero: return _zero
        # Minkowski addition: http://en.wikipedia.org/wiki/Minkowski_addition
        assert isinstance(other, ConvexHull)
        return ConvexHull([a + b for a in self for b in other])

    @classmethod
    def multiplicity(cls, v, m):
        return v if m > 0 else cls.zero

    def star(self):
        # Star diverges for nontrivial hulls.
        raise NotImplementedError()

    def draw(self, label_fmt=None, c='k'):   # pragma: no cover
        "Visualize points with interactive scatter plot browser."
        if len(self.points) <= 2:
            print('[warn] ConvexHull has too few points.')
            return
        points = np.array([(p.x, p.y) for p in self.points])
        hull = scipy.spatial.ConvexHull(points, qhull_options='QJ Pp')
        for simplex in hull.simplices:
            pl.plot(points[simplex, 0], points[simplex, 1], zorder=-1, alpha=0.5, lw=.5, color=c)
        for p in self:
            pl.scatter(p.x, p.y, alpha=0.5, zorder=-1, c=c)
        pl.box(False)
        pl.axis('equal')
        pl.xticks([p.x for p in self], rotation='vertical')
        pl.yticks([p.y for p in self])
        for p in self:
            pl.text(x=p.x, y=p.y, s=str(p.d) if label_fmt is None else label_fmt(p))


def conv(points, jiggle=True):
    "Indices of the convex hull."
    if len(points) <= 2:
        return points
    c = scipy.spatial.ConvexHull(
        np.array([(x.x, x.y) for x in points]),
        qhull_options=('QJ Pp' if jiggle else None),
    )
    return [points[i] for i in c.vertices]


class Point:
    """2D point with backpointer for derivation reconstruction."""

    def __init__(self, x, y, d=None):
        self.x = x
        self.y = y
        self.d = d

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, (self, other))

    def __repr__(self):
        return str(self.d)


_zero = ConvexHull([])
_one = ConvexHull([Point(0, 0)])

ConvexHull.zero = _zero
ConvexHull.one = _one
