class Semiring:
    @classmethod
    def chart(cls):
        return Chart(cls.zero)
    zero = None
    one = None
#    def __lt__(self, other):
#        # Note: this arbitrary ordering is used by `Program.signature`
#        return hash(self) < hash(other)
    def __add__(self, other):
        raise NotImplementedError
    def __mul__(self, other):
        raise NotImplementedError
    @classmethod
    def lift(cls, *args):
        return cls(*args)
    @classmethod
    def multiplicity(cls, v, m):
        if m == 0: return cls.zero
        old = cls.zero
        while True:   # m might be infinite, in that case, run until we saturate
            new = v + old
            m -= 1
            if m == 0 or new == old: return new
            old = new
    @classmethod
    def multiple(cls, m):
        return cls.multiplicity(cls.one, m)
    def __pow__(self, n):
        assert n >= 0
        x = self
        if n == 0:       return self.one
        elif n % 2 == 0: return (x * x) ** (n // 2)
        else:            return x * (x ** (n - 1))
    def star_approx(self, T):
        # star(x) = 1 + x star(x)
        v = self.zero
        for _ in range(T):
            v = self.one + self * v
        return v
    def star_fixpoint(self):
        # star(x) = 1 + x star(x)
        prev = self.zero
        while True:
            curr = self.one + self * prev
            if prev == curr: return curr
            prev = curr
    star = star_fixpoint
    @classmethod
    def approx_zero(cls, x):
        return x == cls.zero
    @classmethod
    def sum(cls, xs):
        y = cls.zero
        for x in xs:
            y += x
        return y
    @classmethod
    def product(cls, xs):
        y = cls.one
        for x in xs:
            y *= x
        return y
    def __round__(self, precision):
        return self
    def metric(self, other):
        return (self != other)


class Chart(dict):

    def __init__(self, zero):
        self.zero = zero
        super(Chart, self).__init__()

    def __missing__(self, key):
        return self.zero

    def copy(self):
        x = Chart(self.zero)
        x.update(self)
        return x

#    def __add__(self, other):
#        assert isinstance(other, Chart), other
#        x = Chart(self.zero)
#        for k,v in self.items():
#            x[k] += v
#        for k,v in other.items():
#            x[k] += v
#        return x
#
#    def __sub__(self, other):
#        assert isinstance(other, Chart), other
#        x = Chart(self.zero)
#        for k,v in self.items():
#            x[k] += v
#        for k,v in other.items():
#            x[k] -= v
#        return x
#
#    def round(self, precision):
#        x = Chart(self.zero)
#        for k,v in self.items():
#            x[k] += round(v, precision)
#        return x


class Wrapped(Semiring):
    def __init__(self, x):
        self.x = x
    def __lt__(self, other): return self.x < other.x
    def __eq__(self, other): return self.x == other.x
    def __hash__(self):      return hash(self.x)
    def __repr__(self):      return f'{self.__class__.__name__}({self.x})'
    @classmethod
    def lift(cls, x):    return cls(x)


def make_semiring(name, plus, times, zero, one, star=None, pp=None, hash=hash, multiplicity=None):
    class SemiringWrapper(Semiring):
        def __init__(self, x, *args):
            self.x = x
        def __add__(self, other): return SemiringWrapper(plus(self.x, other.x))
        def __mul__(self, other):
            if isinstance(other, tuple):
                other = (SemiringWrapper.one if len(other) == 0 else other[0] * other[1:])
            return SemiringWrapper(times(self.x, other.x))
        def __lt__(self, other): return self.x < other.x
        def __eq__(self, other): return self.x == other.x
        def __hash__(self):      return hash(self.x)
        def __repr__(self):
            if pp is None:
                return f'{name}({self.x})'
            else:
                return pp(self.x)
        @classmethod
        def lift(cls, *args):    return cls(*args)

    if star is not None:
        SemiringWrapper.star = lambda self: SemiringWrapper(star(self.x))

    if multiplicity is not None:
        SemiringWrapper.multiplicity = lambda x,m: SemiringWrapper(multiplicity(x.x, m))

    SemiringWrapper.zero = SemiringWrapper(zero)
    SemiringWrapper.one = SemiringWrapper(one)
    SemiringWrapper.__name__ = name
    return SemiringWrapper
