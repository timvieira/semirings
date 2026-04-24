"""First- and second-order expectation semirings (Li & Eisner, 2009).

The first-order Expectation semiring carries pairs (p, r) with addition
coefficientwise and multiplication (p1, r1)*(p2, r2) = (p1*p2, p1*r2 + p2*r1).
This is the truncated dual-number ring P[ε]/(ε²) with ε-coefficient in R.

The second-order version carries 4-tuples (p, r, s, t) with
(p1*p2, p1*r2+p2*r1, p1*s2+p2*s1, p1*t2+p2*t1+r1*s2+r2*s1) — equivalent
to order-2 truncation.

These are candidates for retrofitting as `Jet(P, R, order=k)` aliases in
DESIGN.md Stage 3.
"""
from semirings.base import Semiring
from semirings.logval import LogVal, LogValVector


def make_expectation(P, R):
    class Expectation(Semiring):
        """First-order expectation semiring (p, r) over types (P, R)."""

        def __init__(self, p, r):
            self.p = p
            self.r = r

        def __repr__(self):
            return repr((self.p, self.r))

        def __add__(self, y):
            return Expectation(self.p + y.p, self.r + y.r)

        def __mul__(self, y):
            p1, r1 = self.p, self.r
            p2, r2 = y.p, y.r
            return Expectation(p1 * p2, p1 * r2 + p2 * r1)

        def __eq__(self, other):
            if not isinstance(other, Expectation):
                return NotImplemented
            return self.p == other.p and self.r == other.r

        __hash__ = None  # structural __eq__ with float-ish values → not hashable

        def metric(self, other):
            mp = self.p.metric(other.p) if hasattr(self.p, 'metric') else (0.0 if self.p == other.p else 1.0)
            mr = self.r.metric(other.r) if hasattr(self.r, 'metric') else (0.0 if self.r == other.r else 1.0)
            return max(mp, mr)

    Expectation.zero = Expectation(P.zero, R.zero)
    Expectation.one = Expectation(P.one, R.zero)
    return Expectation


Expectation = make_expectation(LogVal, LogVal)


class SecondOrderExpectation(Semiring):
    """Second-order expectation semiring (p, r, s, t).

    Carries a LogVal scalar `p`, a LogVal `r`, and two LogValVectors `s`, `t`.
    Note: for most computations the inside-outside speed-up is more efficient
    and less memory-intensive than this direct implementation.
    """

    def __init__(self, p, r, s, t):
        self.p = p
        self.r = r
        self.s = s
        self.t = t

    def __repr__(self):
        return repr((self.p, self.r, self.s, self.t))

    def __add__(self, y):
        return SecondOrderExpectation(
            self.p + y.p,
            self.r + y.r,
            self.s + y.s,
            self.t + y.t,
        )

    def __mul__(self, y):
        p1, r1, s1, t1 = self.p, self.r, self.s, self.t
        p2, r2, s2, t2 = y.p, y.r, y.s, y.t
        return SecondOrderExpectation(
            p1 * p2,
            p1 * r2 + p2 * r1,
            p1 * s2 + p2 * s1,
            p1 * t2 + p2 * t1 + r1 * s2 + r2 * s1,
        )

    def __eq__(self, other):
        if not isinstance(other, SecondOrderExpectation):
            return NotImplemented
        return (self.p == other.p and self.r == other.r
                and self.s == other.s and self.t == other.t)

    __hash__ = None  # structural __eq__ with float-ish values → not hashable

    def metric(self, other):
        return max(
            self.p.metric(other.p),
            self.r.metric(other.r),
            self.s.metric(other.s),
            self.t.metric(other.t),
        )


SecondOrderExpectation.zero = SecondOrderExpectation(
    LogVal.zero, LogVal.zero, LogValVector(), LogValVector(),
)
SecondOrderExpectation.one = SecondOrderExpectation(
    LogVal.one, LogVal.zero, LogValVector(), LogValVector(),
)
