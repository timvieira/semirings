from semirings import Semiring
from arsenal.iterextras import merge_roundrobin, fair_product


from semirings.fsa import FSA
from functools import cached_property


class RegularLanguage(Semiring):

    @cached_property
    def fsa(self):
        raise NotImplementedError()

    def __eq__(self, other):
        return isinstance(other, RegularLanguage) and self.fsa.equal(other.fsa)

    def __hash__(self):
        return 0    # terrible hash function :-/

    def star(self):
        if self is NULL: return EPSILON
        if self is EPSILON: return EPSILON
        return Star(self)

    def __iter__(self):
        raise NotImplementedError()

    def __add__(self, other):
        if other is NULL: return self
        if self is NULL: return other
        return Union(self, other)

    def __mul__(self, other):
        if not isinstance(other, RegularLanguage): return other.__rmul__(self)
        # multiplication by zero
        if self is NULL or other is NULL: return NULL
        # multiplication by one.
        if self is EPSILON:  return other
        if other is EPSILON: return self
        return Concat(self, other)


class Union(RegularLanguage):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        super().__init__()
    @cached_property
    def fsa(self):
        return self.x.fsa + self.y.fsa
    def __repr__(self): return f'{self.x} + {self.y}'
    def __iter__(self):
        yield from merge_roundrobin(self.x, self.y)


class Star(RegularLanguage):
    def __init__(self, x):
        self.x = x
        super().__init__()
    @cached_property
    def fsa(self):      return self.x.fsa.star()
    def __repr__(self): return f'({self.x})*'
    def __iter__(self):
        p = EPSILON
        while True:
            yield from p
            p = p * self.x


class Concat(RegularLanguage):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        super().__init__()
    @cached_property
    def fsa(self):      return self.x.fsa * self.y.fsa
    def __repr__(self):
        x = self.x; y = self.y
        if isinstance(x, Union): x = f'({x})'
        if isinstance(y, Union): y = f'({y})'
        return f'{x}⋅{y}'
    def __iter__(self):
        for x,y in fair_product(self.x, self.y):
            yield x * y


class Symbol(RegularLanguage):
    def __init__(self, x):
        self.x = x
        super().__init__()
    @cached_property
    def fsa(self):      return FSA.lift(self.x)
    def __repr__(self): return f'{self.x}'
    def __iter__(self):
        yield self
#    def __call__(self, *args):
#        # check for NULL in args.
#        #if any(x == NULL for x in args): return NULL
#        return LabeledProduct(self, *args)


#class LabeledProduct(RegularLanguage):
#    def __init__(self, fn, *args):
#        self.fn = fn; self.args = args
#        super().__init__()
#    def __repr__(self):
#        if len(self.args) == 0:
#            return f'{self.fn}'
#        elif len(self.args) == 1:
#            #return f'{self.fn}⋅{self.args[0]}'  # to render as a string
#            return f'{self.fn}⟨{self.args[0]}⟩'
#        else:
#            return f'{self.fn}⟨{",".join(str(x) for x in self.args)}⟩'
#    def __iter__(self):
#        for a in fair_product(*self.args):
#            yield LabeledProduct(self.fn, *a)


class Zero(RegularLanguage):
    @cached_property
    def fsa(self):      return FSA.zero
    def __repr__(self): return '∅'
    def __iter__(self):
        if 0: yield
        return


class One(RegularLanguage):
    @cached_property
    def fsa(self):      return FSA.one
    def __repr__(self): return 'ε'
    def __iter__(self):
        yield self


EPSILON = RegularLanguage.one = One()
NULL = RegularLanguage.zero = Zero()
RegularLanguage.lift = Symbol


def evaluate(x, subst, semiring):
    def f(x):
        if   isinstance(x, Symbol): return subst[x]
        elif isinstance(x, Union):  return f(x.x) + f(x.y)
        elif isinstance(x, Concat): return f(x.x) * f(x.y)
        elif isinstance(x, Star):   return semiring.star(f(x.x))
        elif isinstance(x, Zero):   return semiring.zero
        elif isinstance(x, One):    return semiring.one
        else:                       raise ValueError(x)
    return f(x)


def test():
#    from arsenal.iterextras import take
#    from collections import Counter

#    def check(a, b, T=None):
#        assert a == b
#        A = list(take(T, a))
#        B = list(take(T, b))
#        assert Counter(A) == Counter(B), f'\n\n{A}\n{B}\n'

    a,b,c,d = map(Symbol, 'abcd')

    assert set(NULL) == set()
    assert set(EPSILON) == {EPSILON}
    assert set(NULL.star()) == {EPSILON}
    assert set(EPSILON.star()) == {EPSILON}

    assert a != b
    assert a == a

    assert (a*a)*a == a*(a*a)

    assert a*c != a*d, [a*c, a*d]

    ab = a*c + b*c + a*d + b*d
    assert len(set(ab)) == 4, list(ab)

    ab = (a + b) * (c + d)
    assert len(set(a + b)) == 2
    assert len(set(c + d)) == 2
    assert len(set(ab)) == 4, set(ab)

    assert a.star() == EPSILON + a * a.star()
    assert (a * b).star() == EPSILON + a * (b * a).star() * b
    assert (a + b).star() == a.star() * (a.star() * b).star() * a.star()

    assert hash(a.star()) == hash(a.star())
#    assert hash(a.star()) != hash(a)
#    assert hash(a + b) != hash(a * b)

#    for x in take(25, (a.star() + b.star()) * (c.star() + d.star())):
#        print(x)

    # Tests for labeled products
#    f = Symbol('f')
#    g = Symbol('g')
#    check(f(NULL), NULL)
#
#    check(f(a + b), f(a) + f(b))
#
#    check(g(f(a + b), c + d),
#          g(f(a), c + d) + g(f(b), c + d))
#
#    check(g(f(NULL + a + b), c + d),
#          g(f(a), c) + g(f(b), c) + g(f(a), d) + g(f(b), d))

    #check(EPSILON(a), a)

    print('test: pass')


if __name__ == '__main__':
    test()
