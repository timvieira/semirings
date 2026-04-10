"""
Lazy k-best semiring: lazily enumerates the K-best derivations in sorted
order without fixing K in advance.

Use `make_lazysort_semiring` to create variants with different base operations
and orderings (e.g., max-times, min-plus, lexicographic).
"""

import operator
from functools import reduce
from arsenal.iterextras import sorted_union, sorted_product, take
from . import base


def make_lazysort_semiring(name, base_times, base_one, base_lt):
    """
    Build a lazy k-best semiring over a totally ordered monoid.

    Lazily enumerates derivations in sorted order without fixing K in advance.
    Uses lazy sorted merge and sorted cross-product over streams (Huang and
    Chiang, 2005).

    Parameters:
        name: name for the returned element class
        base_times: (a, b) -> a*b, the monoid operation combining scores
        base_one: identity element for base_times
        base_lt: (a, b) -> bool, True if a should be enumerated before b.
                 For the star operation to be well-defined, base_one must be
                 the first element under this ordering (i.e., base_lt(base_one, x)
                 for all reachable x != base_one).

    Returns:
        A Semiring subclass whose instances are lazy sorted streams.
        Construct elements with Cls(score, data).

    Examples::

        # Max-times (Viterbi k-best): probabilities, largest first
        LazyMaxTimes = make_lazysort_semiring('LazyMaxTimes',
            base_times=operator.mul, base_one=1, base_lt=lambda a, b: a > b)

        # Min-plus (tropical k-best): costs, smallest first
        LazyMinPlus = make_lazysort_semiring('LazyMinPlus',
            base_times=operator.add, base_one=0, base_lt=lambda a, b: a < b)

        # Shortlex: string concatenation, shortest-first then lexicographic
        LazyShortLex = make_lazysort_semiring('LazyShortLex',
            base_times=operator.add, base_one='',
            base_lt=lambda a, b: (len(a), a) < (len(b), b))

        a = LazyShortLex('a', 'a')
        b = LazyShortLex('b', 'b')
        list(take(8, (a + b).star()))
        # => ['', 'a', 'b', 'aa', 'ab', 'ba', 'bb', 'aaa']
    """

    def _prod_fn(items):
        return reduce(operator.mul, items)

    class _Base(base.Semiring):
        def __add__(self, other):
            if self is _zero: return other
            if other is _zero: return self
            return _Sum(self, other)
        def __mul__(self, other):
            if self is _one: return other
            if other is _one: return self
            if self is _zero: return _zero
            if other is _zero: return _zero
            return _Prod(self, other)
        def star(self):
            if self is _zero: return _one
            return _Star(self)
        @classmethod
        def multiplicity(cls, v, m):
            return cls.zero if m == 0 else _Multiplicity(v, m)
        def take(self, K):
            return take(K, self)

    # Notice that multiplication is neither associative nor commutative
    # because of the backpointers.  This is how we get out a derivation tree
    # (with parentheses).  If we ignore backpointers, then we are associative.
    class _Elem(_Base):
        def __init__(self, score, data):
            self.score = score
            self.data = data
            super().__init__()
        def __mul__(self, other):
            if self is _one: return other
            if other is _one: return self
            if self is _zero: return _zero
            if other is _zero: return _zero
            if isinstance(other, _Elem):
                return _Elem(base_times(self.score, other.score),
                             [self.data, other.data])
            return super().__mul__(other)
        def __lt__(self, other):
            return base_lt(self.score, other.score)
        def __eq__(self, other):
            return isinstance(other, _Elem) and self.score == other.score
        def __hash__(self):
            return hash(self.score)
        def __iter__(self):
            yield self
        def __repr__(self):
            return repr((self.score, self.data))
        def flat_data_list(self):
            return flatten(self.data)

    class _Multiplicity(_Base):
        def __init__(self, a, m):
            assert isinstance(a, _Base), a
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

    class _Sum(_Base):
        def __init__(self, a, b):
            assert isinstance(a, _Base) and isinstance(b, _Base), [a, b]
            self.a = a
            self.b = b
            super().__init__()
        def __iter__(self):
            yield from sorted_union(self.a, self.b)
        def __repr__(self):
            return 'sum(...)'

    class _Prod(_Base):
        def __init__(self, a, b):
            assert isinstance(a, _Base) and isinstance(b, _Base), [a, b]
            self.a = a
            self.b = b
            super().__init__()
        def __iter__(self):
            yield from sorted_product(_prod_fn, self.a, self.b)
        def __repr__(self):
            return f'({self.a} * {self.b})'

    class _Star(_Base):
        def __init__(self, a):
            assert isinstance(a, _Base), a
            self.a = a
            super().__init__()
        def __iter__(self):
            if self.a is _zero:
                yield _one
                return
            # a* = 1 ⊕ a ⊗ a*: yield one first to break the recursion,
            # then lazily merge the remaining powers via sorted_product.
            # When sorted_product wraps self in buf_iter and reads self[0],
            # it gets the already-yielded one from the buffer.
            #
            # Correctness requires that one is the first element in the
            # enumeration order (i.e., base_lt(base_one, x) for all
            # reachable x != base_one).
            yield _one
            yield from sorted_product(_prod_fn, self.a, self)
        def __repr__(self):
            return f'star({self.a})'

    class _Zero(_Base):
        def __init__(self):
            super().__init__()
        def __iter__(self):
            return iter([])
        def __repr__(self):
            return 'NULL'

    _zero = _Zero()
    _one = _Elem(base_one, ())

    _Elem.__name__ = name
    _Elem.__qualname__ = name
    _Elem.zero = _zero
    _Elem.one = _one

    return _Elem


def flatten(S):
    if not isinstance(S, list):
        return [S]
    else:
        tmp = []
        for x in S:
            tmp.extend(flatten(x))
        return tmp


# Default: max-times (Viterbi k-best), enumerates highest scores first
LazySort = make_lazysort_semiring(
    name='LazySort',
    base_times=operator.mul,
    base_one=1,
    base_lt=lambda a, b: a > b,
)

zero = LazySort.zero
one = LazySort.one
