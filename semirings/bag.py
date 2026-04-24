from collections import defaultdict

from semirings.base import Semiring


class Bag(Semiring):
    """Bag semiring: multiset of elements with integer multiplicities.

    Addition merges bags (sum multiplicities). Multiplication takes the
    Cartesian product, tagging each resulting element with the pair of its
    sources. Useful for provenance tracking where derivations are ordered
    tuples and multiplicities matter.
    """

    def __init__(self, bag=None):
        self.bag = defaultdict(int)
        if bag is not None:
            self.bag.update(bag)

    def __add__(self, other):
        new = Bag()
        for k, v in self:
            new.bag[k] += v
        for k, v in other:
            new.bag[k] += v
        return new

    def __sub__(self, other):
        new = Bag()
        for k, v in self:
            new.bag[k] += v
        for k, v in other:
            new.bag[k] -= v
        return new

    def __mul__(self, other):
        return Bag({(x, y): xm * ym for x, xm in self for y, ym in other})

    def __iter__(self):
        return iter(self.bag.items())

    def __eq__(self, other):
        return isinstance(other, Bag) and dict(self.bag) == dict(other.bag)

    def __hash__(self):
        return hash(frozenset(self.bag.items()))

    def __repr__(self):
        return f'Bag({dict(self.bag)})'

    @classmethod
    def lift(cls, w, d):
        return Bag({d: w})


Bag.zero = Bag()
Bag.one = Bag({(): 1})
