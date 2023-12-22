import numpy as np
from arsenal.iterextras import sorted_union, sorted_product, take
from . import base

# TODO: create a LazySort builder which is parametric in the underlying ops,
# just like expectation semiring.  This will reuse code for max-sum, min-sum,
# sum-prod, and even lexicographic semirings.

class _LazySort(base.Semiring):
    def __add__(self, other):
        if self is zero: return other
        if other is zero: return self
        return Sum(self, other)
    def __mul__(self, other):
        if self is one: return other
        if other is one: return self
        if self is zero: return zero
        if other is zero: return zero
        return Prod(self, other)
    def star(self):
        if self is zero: return one
        return Star(self)
    @classmethod
    def multiplicity(cls, v, m):
        return cls.zero if m == 0 else Multiplicity(v, m)
    def take(self, K):
        return take(K, self)


class Multiplicity(_LazySort):
    def __init__(self, a, m):
        assert isinstance(a, _LazySort), a
        self.a = a
        self.m = m
        super().__init__()
    def __iter__(self):
        for v in self.a:
            i = 0
            while i < self.m:
                yield v
                i += 1
    def __repr__(self):
        return f'multiplicity({self.a},{self.m})'


# Notice that in this multiplication is neither associative, nor commutative
# because of the backpointers. This is how we get out a derivation tree (with
# parentheses).  If we ignore backpointers, then we are are associative.
class LazySort(_LazySort):
    def __init__(self, score, data):
        #assert isinstance(score, (int,float)), score
        self.score = score
        self.data = data
        super().__init__()

    # TODO: it is a strange choice to have a different product operation than
    # _LazySort - can we change that?  I think the issue is actually that there
    # is a "base" product that is passed to sorted_product and all that.  And,
    # we've sort of conflated things when we should have kept them separate.
    # The objects that we iterate over should be instances of an ordered base
    # semiring.  The details of backpointers (like we see below) would be
    # handled there too.
    def __mul__(self, other):
        if self is one: return other
        if other is one: return self
        if self is zero: return zero
        if other is zero: return zero
        if isinstance(other, LazySort):
            return LazySort(self.score * other.score, [self.data, other.data])
        else:
            return super().__mul__(other)
    def __lt__(self, other):
        return other.score < self.score    # Warning: this is backwards!
    def __iter__(self):
        yield self
    def __eq__(self, other):
        return isinstance(other, LazySort) and self.score == other.score
    def __repr__(self):
        return repr((self.score, self.data))
    def flat_data_list(self):
        return flatten(self.data)


class Sum(_LazySort):
    def __init__(self, a, b):
        assert isinstance(a, _LazySort) and isinstance(b, _LazySort), [a,b]
        self.a = a
        self.b = b
        super().__init__()
    def __iter__(self):
        yield from sorted_union(self.a, self.b)
    def __repr__(self):
#        return f'({self.a} + {self.b})'
        return 'sum(...)'


class Prod(_LazySort):
    def __init__(self, a, b):
        assert isinstance(a, _LazySort) and isinstance(b, _LazySort), [a,b]
        self.a = a
        self.b = b
        super().__init__()
    def __iter__(self):
        yield from sorted_product(np.prod, self.a, self.b)
    def __repr__(self):
        return f'({self.a} * {self.b})'


class Star(_LazySort):
    def __init__(self, a):
        assert isinstance(a, _LazySort), a
        self.a = a
        super().__init__()
    def __iter__(self):
        if self.a is zero:
            yield one
            return
        v = one
        while True:
            print('star...', self.a)
            yield from v
            v *= self.a
    def __repr__(self):
        return f'star({self.a})'


class Zero(_LazySort):
    def __init__(self):
        self.score = 0
        super().__init__()
    def __iter__(self):
        if 0: yield
    def __repr__(self):
        return 'NULL'


zero = Zero()
one = LazySort(1, ())

LazySort.zero = zero
LazySort.one = one


def flatten(S):
    if not isinstance(S, list):
        return [S]
    else:
        tmp = []
        for x in S:
            tmp.extend(flatten(x))
        return tmp
