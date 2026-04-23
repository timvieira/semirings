class Semiring:
    """Base class for semirings.

    Subclasses must define ``zero``, ``one``, ``__add__``, and ``__mul__``.

    A ``metric(self, other)`` method provides a distance function between
    semiring elements.  It must satisfy the metric axioms:

    1. ``d(x, y) >= 0``  (non-negativity)
    2. ``d(x, y) == 0`` iff ``x == y``  (identity of indiscernibles)
    3. ``d(x, y) == d(y, x)``  (symmetry)
    4. ``d(x, z) <= d(x, y) + d(y, z)``  (triangle inequality)

    The default implementation returns the discrete metric (0 if equal, 1
    otherwise).  Numeric semirings should override with a meaningful distance
    (e.g., ``abs(self.score - other.score)``).
    """
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
    def star_doubling(self):
        # star(x) = (1 + x) star(x²); doubles the number of series terms per step.
        p = self.one
        x = self
        while True:
            curr = p + p * x
            if curr == p: return curr
            p = curr
            x = x * x
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
    @classmethod
    def assert_equal(cls, x, *ys, tol=1e-10):
        if not all(x == y or cls.metric(x, y) <= tol for y in ys):
            assert False, f'{cls.__name__}: {x} != {ys}'


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
    def metric(self, other):
        if self.x == other.x: return 0
        if isinstance(self.x, (int, float)) and isinstance(other.x, (int, float)):
            d = abs(self.x - other.x)
            s = max(1, abs(self.x), abs(other.x))
            if d == float('inf'): return float('inf')
            return d / s
        return self != other
    @classmethod
    def lift(cls, x):        return cls(x)


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
        def __eq__(self, other): return isinstance(other, SemiringWrapper) and self.x == other.x
        def __hash__(self):      return hash(self.x)
        def __repr__(self):
            if pp is None:
                return f'{name}({self.x})'
            else:
                return pp(self.x)
        def metric(self, other):
            if self.x == other.x: return 0
            if isinstance(self.x, (int, float)) and isinstance(other.x, (int, float)):
                d = abs(self.x - other.x)
                s = max(1, abs(self.x), abs(other.x))
                if d == float('inf'): return float('inf')
                return d / s
            return self != other
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



def _check_axioms_unary(S, A, hash_=True, check_multiplicity=True, star=True,
                        assoc=True, **_):
    assert A == A
    if hash_:
        assert hash(A) == hash(A)
    S.assert_equal(A + S.zero, A, S.zero + A)
    if assoc:
        S.assert_equal(A * S.one, A, S.one * A)
    S.assert_equal(A * S.zero, S.zero, S.zero * A)
    if check_multiplicity:
        S.assert_equal(S.multiplicity(A, 0), S.zero)
        S.assert_equal(S.multiplicity(A, 1), A)
        S.assert_equal(S.multiplicity(A, 2), A + A)
        S.assert_equal(S.multiplicity(A, 3), A + A + A)
    if star:
        S.assert_equal(S.star(A),
                       S.one + S.star(A) * A,
                       S.one + A * S.star(A))


def _check_axioms_binary(S, A, B, **_):
    S.assert_equal(A + B, B + A)


def _check_axioms_ternary(S, A, B, C, left_distrib=True, right_distrib=True,
                         assoc=True, **_):
    S.assert_equal((A + B) + C, A + (B + C))
    if assoc:
        S.assert_equal((A * B) * C, A * (B * C))
    if left_distrib:
        S.assert_equal(A * (B + C), A * B + A * C)
    if right_distrib:
        S.assert_equal((B + C) * A, B * A + C * A)


def check_axioms_samples(S, samples, **kwargs):
    # Hoist unary/binary axioms out of the triple loop: O(n) + O(n²) + O(n³)
    # calls instead of O(n³) with redundant re-checks of the cheaper axioms.
    try:
        for A in samples:
            _check_axioms_unary(S, A, **kwargs)
        for A in samples:
            for B in samples:
                _check_axioms_binary(S, A, B, **kwargs)
        for A in samples:
            for B in samples:
                for C in samples:
                    _check_axioms_ternary(S, A, B, C, **kwargs)
    except AssertionError:
        print()
        print(S.__name__, 'samples:', samples)
        raise


def check_axioms(S, A, B, C, **kwargs):
    try:
        _check_axioms_unary(S, A, **kwargs)
        _check_axioms_binary(S, A, B, **kwargs)
        _check_axioms_ternary(S, A, B, C, **kwargs)
    except AssertionError:
        print()
        print(S.__name__, A, B, C)
        raise


def check_metric_axioms(S, samples):
    for a in samples:
        # d(a,a) == 0
        d_aa = S.metric(a, a)
        assert d_aa == 0 or d_aa != d_aa, f'metric({a}, {a}) = {d_aa} (expected 0 or nan)'
        for b in samples:
            d_ab = S.metric(a, b)
            if a == b or d_ab != d_ab:
                pass  # equal or nan — skip
            else:
                # Non-negativity
                assert d_ab >= 0, f'metric({a}, {b}) = {d_ab} < 0'
                # Symmetry
                assert d_ab == S.metric(b, a), f'metric({a}, {b}) = {d_ab} != metric({b}, {a}) = {S.metric(b, a)}'
            # Triangle inequality (skip if any distance is nan)
            for c in samples:
                d_ac = S.metric(a, c)
                d_bc = S.metric(b, c)
                if d_ac != d_ac or d_ab != d_ab or d_bc != d_bc:
                    continue  # nan — metric is undefined here
                assert d_ac <= d_ab + d_bc + 1e-15, (
                    f'triangle inequality: metric({a}, {c}) = {d_ac}'
                    f' > metric({a}, {b}) + metric({b}, {c}) = {d_ab} + {d_bc} = {d_ab + d_bc}'
                )
