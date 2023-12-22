"""
MERT semiring
"""
import numpy as np
import pylab as pl
import scipy.spatial
from semirings import Semiring

# TODO: there is nothing specific about the convexhull semiring and being
# two-dimensional.

# TODO: I read somewhere that the union of convex hulls can be computed in
# linear time (no log factor).

# TODO: convex hull operation commutes with minkowski sum (no need for a second
# call to convex hull)

class ConvexHull(Semiring):
    """Convex hull semiring (in two-dimensions).

    Each element of this semiring is a convex hull of `Points`.

    Based closely on Chris Dyer's 2013 arxiv paper
    (http://arxiv.org/pdf/1307.3675.pdf)

    Coded for clarity, not efficiency.

    """

    def __init__(self, points):
        self.points = conv(points)

    def __iter__(self):
        return iter(self.points)

    def __eq__(self, other):
        return set(self.points) == set(other.points)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.points})'

    def __hash__(self):
        return hash(frozenset(self.points))

    def __add__(self, other):
        assert isinstance(other, ConvexHull)
        return ConvexHull(list(self) + list(other))

    def __mul__(self, other):
        # http://en.wikipedia.org/wiki/Minkowski_addition
        assert isinstance(other, ConvexHull)
        return ConvexHull([a + b
                           for a in self
                           for b in other])

    @classmethod
    def multiplicity(cls, v, m):
        return v if m > 0 else cls.zero

    # The star operation diverges
    def star(self):
        #return self.one + self
        raise NotImplementedError()

    def draw(self, label_fmt=None, c='k'):   # pragma: no cover
        "Visualize points with interactive scatter plot browser."
        if len(self.points) <= 2:
            print('[warn] ConvexHull has too few points empty.')
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
    else:
        c = scipy.spatial.ConvexHull(
            np.array([(x.x, x.y) for x in points]),
            qhull_options=('QJ Pp' if jiggle else None)
        )
        return [points[i] for i in c.vertices]


class Point:
    """
    Two-dimensional point with backpointers so that we can reconstruct the
    derivation.
    """

    def __init__(self, x, y, d=None):
        self.x = x
        self.y = y
        self.d = d

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other):
        return Point(self.x + other.x,
                     self.y + other.y,
                     (self, other))

    def __repr__(self):
        #d = '' if self.d is None else f', {self.d}'
        #d = t._pformat_flat(nodesep='', parens='()', quotes=False)
        #return f'Point({self.x}, {self.y}{d})'
        return str(self.d)


zero = ConvexHull([])
one = ConvexHull([Point(0,0)])

ConvexHull.zero = zero
ConvexHull.one = one
