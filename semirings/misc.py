from semirings.base import Semiring, Wrapped


class Why(Wrapped):
    """
    Set of combinations of tuples (annotations) needed for a tuple to exist in
    the answer to a query (or derive a given item).
    """
    def __add__(self, other): return Why(self.x | other.x)
    def __mul__(self, other): return Why(frozenset(x | y
                                                   for x in self.x
                                                   for y in other.x))
    @classmethod
    def lift(cls, x): return  cls(frozenset([frozenset([x])]))

Why.zero = Why(frozenset())
Why.one = Why(frozenset([frozenset([])]))


class Lineage(Wrapped):
    "Set of reachable annotations"

    def __add__(self, other):
        if self == self.zero: return other
        if other == self.zero: return self
        return Lineage(self.x | other.x)
    def __mul__(self, other):
        if self == self.zero: return self.zero
        if other == self.zero: return self.zero
        return Lineage(self.x | other.x)
    @classmethod
    def lift(cls, x): return  cls(frozenset([x]))

Lineage.zero = Lineage(None)
Lineage.one = Lineage(frozenset())


class Bridge(Semiring):
    "Edge bridge"
    def __init__(self, xs):  # set of edges
        self.xs = frozenset(xs) if xs is not None else xs
    def __add__(self, other):
        if self == self.zero: return other
        if other == self.zero: return self
        return Bridge(self.xs & other.xs)
    def __mul__(self, other):
        if self == self.zero: return self.zero
        if other == self.zero: return self.zero
        return Bridge(self.xs | other.xs)
    def __eq__(self, other):
        return self.xs == other.xs
    def __hash__(self):
        return hash(self.xs)
    @classmethod
    def lift(cls, x):
        return cls({x})
    def __repr__(self):
        return repr(self.xs)
Bridge.zero = Bridge(None)
Bridge.one = Bridge(set())


class VBridge(Bridge):
    "Vertex bridge"
    @classmethod
    def lift(cls, x):
        (i,j) = x
        return cls({i,j})


class Lukasiewicz(Wrapped):
    """The Łukasiewicz semiring: the closed interval [0,1] with addition given by
    max(a,b) and multiplication given by max(0, a + b - 1) appears in
    multi-valued logic.
    """
    def __add__(self, other): return max(self, other)
    def __mul__(self, other): return Lukasiewicz(max(0, self.x + other.x - 1))
Lukasiewicz.zero = Lukasiewicz(0)
Lukasiewicz.one = Lukasiewicz(1)


def make_set(universe, emptyset):
    class Set(Wrapped):
        """
        Semiring of sets
        """
        def __add__(self, other): return Set(self.x | other.x)
        def __mul__(self, other): return Set(self.x & other.x)
    Set.zero = Set(emptyset)
    Set.one = Set(universe)
    return Set


from math import gcd
lcm = lambda x,y: abs(x * y) // gcd(x,y) if x * y != 0 else 0
class Division(Wrapped):
    def __add__(self, other): return Division(gcd(self.x, other.x))
    def __mul__(self, other): return Division(lcm(self.x, other.x))
Division.zero = Division(0)
Division.one = Division(1)


#class WrappedBackpointers(Semiring):
#    def __init__(self, x, d=None):
#        self.x = x
#        self.d = d   # ignored in equality test
#    def __lt__(self, other): return self.x < other.x
#    def __eq__(self, other): return self.x == other.x
#    def __hash__(self):      return hash(self.x)
#    def __repr__(self):      return f'{self.__class__.__name__}({self.x})'
#    @classmethod
#    def lift(cls, *args):    return cls(*args)


def minmax(zero, one):
    class minmax(Wrapped):
        def __add__(self, other):
            if self == zero: return other
            if other == zero: return self
            return min(self, other)
        def __mul__(self, other):
            if self == zero: return zero
            if other == zero: return zero
            if self == one: return other
            if other == one: return self
            return max(self, other)
        def multiplicity(self, m):
            return self if m > 0 else self.zero
    minmax.zero = zero = minmax(zero)
    minmax.one = one = minmax(one)
    return minmax


def maxmin(zero, one):

    class maxmin(Wrapped):
        def __add__(self, other):
            if self == zero: return other
            if other == zero: return self
            return max(self, other)
        def __mul__(self, other):
            if self == zero: return zero
            if other == zero: return zero
            if self == one: return other
            if other == one: return self
            return min(self, other)
        def multiplicity(self, m): return self if m > 0 else self.zero

    maxmin.zero = zero = maxmin(zero)
    maxmin.one = one = maxmin(one)
    return maxmin


Bottleneck = maxmin(float('-inf'), float('inf'))


def dual(base):
    class Dual:
        def __init__(self, p, r):
            self.p = p
            self.r = r
        def __eq__(self, y):  return self.p == y.p and self.r == y.r
        def __hash__(self):   return hash((self.p, self.r))
        def __iter__(self):   return iter((self.p, self.r))
        def __repr__(self):   return f'⟨{self.p}, {self.r}⟩'
        def __add__(self, y): return self.__class__(self.p + y.p, self.r + y.r)
        def star(self):
            ps = base.star(self.p)
            return self.__class__(ps, ps * self.r * ps)
        def __mul__(self, y):
            p1,r1 = self.p, self.r
            p2,r2 = y.p, y.r
            return self.__class__(p1*p2, p1*r2 + r1*p2)
        def multiplicity(self, m):
            return self.__class__(base.multiplicity(self.p, m),
                                  base.multiplicity(self.r, m))

    Dual.zero = Dual(base.zero, base.zero)
    Dual.one = Dual(base.one, base.zero)
    return Dual


class String(Wrapped):

    def star(self): return String.one

    def __add__(self, other):
        if other is self.zero: return self
        if self is self.zero: return other
        return String(longest_common_prefix(self.x, other.x))

    def __mul__(self, other):
        if other is self.one:  return self
        if self is self.one:   return other
        if other is self.zero: return self.zero
        if self is self.zero:  return self.zero
        return String(self.x + other.x)

    def __truediv__(self, other):
        prefix = longest_common_prefix(self.x, other.x)
        return String(self.x[len(prefix):])


def longest_common_prefix(xs, ys):
    "computes the longest common prefix"
    position = -1
    for n in range(min(len(xs), len(ys))):
        if xs[n] == ys[n]:
            position = n
        else:
            break
    return xs[:position+1]


String.zero = String(None)   # unique "infinity" string
String.one = String("")


class ThreeValuedLogic(Wrapped):
    def __add__(self, other): return max(self, other)
    def __mul__(self, other): return min(self, other)
ThreeValuedLogic.zero = ThreeValuedLogic(-1)
ThreeValuedLogic.one = ThreeValuedLogic(+1)
ThreeValuedLogic.unk = ThreeValuedLogic(0)
