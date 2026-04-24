"""Free-semiring machinery.

`FreeExpr` is a magma representation of a sum-product-star expression tree:
binary `+` and `*` at internal nodes, unary `star`, generators at leaves.
Equality is syntactic (tree equality), so `FreeExpr` does NOT satisfy the
semiring axioms — `(a+b)+c` is not equal to `a+(b+c)` as trees, etc. It is
useful as a trace of a computation and as a container that other interpreters
(lazysort, sample, weight, maxtimes, backprop) can walk.

For a proper free (closed) semiring with decidable equality, use a WFSA-backed
realization — see DESIGN.md Part I.
"""
from collections import defaultdict

import numpy as np
from arsenal.cache import memoize
from arsenal.iterextras import sorted_union, sorted_product

from semirings.base import Semiring


class FreeExpr(Semiring):
    def __init__(self, *args):
        self.args = args

    def __lt__(self, other):
        if isinstance(other, FreeExpr):
            return self.args < other.args
        else:
            car = self.args[0] if len(self.args) > 1 else None
            return car < other

    def __add__(self, other):
        if self is _zero: return other
        if other is _zero: return self
        return Sum(self, other)

    def __mul__(self, other):
        if self is _one: return other
        if other is _one: return self
        if self is _zero: return _zero
        if other is _zero: return _zero
        return Prod(self, other)

    def star(self):
        return Star(self)

    @staticmethod
    def lift(w, d):
        return FreeExpr(w, d)

    def __repr__(self):
        return f'{self.__class__.__name__}({", ".join(map(repr, self.args))})'


FreeExpr.zero = _zero = FreeExpr()
FreeExpr.one = _one = FreeExpr()


class Sum(FreeExpr):
    def __repr__(self):
        y, z = self.args
        return f'({y} + {z})'


class Prod(FreeExpr):
    def __repr__(self):
        y, z = self.args
        return f'({y} * {z})'


class Star(FreeExpr):
    def __repr__(self):
        [y] = self.args
        return f'star({y})'


def lazysort(x):
    if isinstance(x, Sum):
        y, z = x.args
        yield from sorted_union(lazysort(y), lazysort(z))
    elif isinstance(x, Prod):
        y, z = x.args
        yield from sorted_product(np.prod, lazysort(y), lazysort(z))
    elif isinstance(x, Star):
        [y] = x.args
        v = _one
        while True:
            yield from lazysort(v)
            v *= y
    elif isinstance(x, FreeExpr):
        yield x
    else:
        assert False


def sample(x):
    """Sample a derivation from FreeExpr `x` weighted by `weight()`."""
    if isinstance(x, Sum):
        y, z = x.args
        if np.random.uniform(0, 1) * weight(x) <= weight(y):
            return sample(y)
        else:
            return sample(z)
    elif isinstance(x, Prod):
        y, z = x.args
        return sample(y) * sample(z)
    elif isinstance(x, Star):
        # Geometric-distribution-style repetition sampling.
        zs = _one
        [y] = x.args
        while True:
            if np.random.uniform() <= weight(y):   # XXX: Need to normalize this.
                zs = zs * sample(y)
            else:
                return zs
    elif isinstance(x, FreeExpr):
        return x
    else:
        assert False


@memoize
def weight(x):
    if isinstance(x, Sum):
        y, z = x.args
        return weight(y) + weight(z)
    elif isinstance(x, Prod):
        y, z = x.args
        return weight(y) * weight(z)
    elif isinstance(x, Star):
        [y] = x.args
        return 1 / (1 - weight(y))
    elif isinstance(x, FreeExpr):
        # Weights must be specified as base cases. By convention they should
        # be the first element of the expression's args.
        return x.args[0]
    else:
        assert False


@memoize
def maxtimes(x):
    if isinstance(x, Sum):
        y, z = x.args
        yy = maxtimes(y)
        zz = maxtimes(z)
        if weight(yy) > weight(zz):
            return yy
        else:
            return zz
    elif isinstance(x, Prod):
        y, z = x.args
        return maxtimes(y) * maxtimes(z)
    elif isinstance(x, Star):
        [y] = x.args
        return maxtimes(y)
    elif isinstance(x, FreeExpr):
        return x
    else:
        assert False


def backprop(expr):
    """Transform a sum-product graph into its adjoint graph."""
    adj = defaultdict(lambda: _zero)
    adj[expr] = _one
    for x in reversed(list(toposort(expr))):
        v = adj[x]
        if isinstance(x, Sum):
            X, Y = x.args
            adj[X] += v
            adj[Y] += v
        elif isinstance(x, Prod):
            X, Y = x.args
            adj[X] += v * Y
            adj[Y] += X * v
    return adj


def toposort(root):
    visited = set()

    def t(v):
        if v not in visited:
            visited.add(v)
            if isinstance(v, FreeExpr):
                for y in v.args:
                    yield from t(y)
            yield v

    yield from t(root)
