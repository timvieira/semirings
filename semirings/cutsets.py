from . import base


class CutSets(base.Semiring):
    """Set of cut sets

    Implements: Martelli (1976)
    "A Gaussian Elimination Algorithm for the Enumeration of Cut Sets in a Graph"
    https://dl.acm.org/doi/abs/10.1145/321921.321928

    """
    def __init__(self, *sets):
        assert all(isinstance(x, (frozenset, set)) for x in sets)
        self.sets = frozenset(min_sets(map(frozenset, sets)))  # project onto minimal set
    def __repr__(self):
        return '{%s}' % (', '.join(repr(set(x) or {}) for x in self))
    def __hash__(self):
        return hash(self.sets)
    def __eq__(self, other):
        return self.sets == other.sets
    def __mul__(self, other):
        return CutSets(*(self.sets | other.sets))
    def __add__(self, other):
        return CutSets(*{
            (x | y)
            for x in self
            for y in other
        })
    def __iter__(self):
        return iter(self.sets)
    def star(self):
        return self.one #+ self
    @classmethod
    def lift(cls, e): return cls({e})


def min_sets(xs):
    """
    Add r to the collection.
    Returns True if r was subsumed by a more general rule in the program
    Returns False otherwise.
    """
    ys = []
    for x in xs:
        add(ys, x)
    return ys


def add(xs, r):
    # This implementation is pretty inefficient as it is based on linear scan
    rm = []
    for i, s in enumerate(xs):
        # branch r is already subsumed by branch s
        if s < r: return True
        # new branch subsumes existing branch, will be deleted in favor of the
        # more general branch.
        if r < s: rm.append(i)
    for i in reversed(sorted(rm)):
        del xs[i]
    xs.append(r)
    return False


CutSets.zero = CutSets(set())
CutSets.one = CutSets()
