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
            # Jet product: (p1 + r1·ε)(p2 + r2·ε) = p1·p2 + (p1·r2 + r1·p2)·ε.
            # Order matters when P is non-commutative (e.g., RegularLanguage).
            p1, r1 = self.p, self.r
            p2, r2 = y.p, y.r
            return Expectation(p1 * p2, p1 * r2 + r1 * p2)

        def star(self):
            # First-order jet star, non-commutative.
            # In S[ε]/(ε²): star(p + rε) = (p*, p* r p*) since
            # x^n = (p^n, sum_{i=0}^{n-1} p^i r p^{n-1-i}) and summing over n
            # telescopes to p* r p*.
            ps = self.p.star()
            return Expectation(ps, ps * self.r * ps)

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
        # Order-2 jet product with independent commuting infinitesimals ε, δ:
        # (p + rε + sδ + tεδ)(p' + r'ε + s'δ + t'εδ)
        #   = pp' + (pr' + rp')ε + (ps' + sp')δ + (pt' + tp' + rs' + sr')εδ
        # Order matters when P is non-commutative.
        p1, r1, s1, t1 = self.p, self.r, self.s, self.t
        p2, r2, s2, t2 = y.p, y.r, y.s, y.t
        return SecondOrderExpectation(
            p1 * p2,
            p1 * r2 + r1 * p2,
            p1 * s2 + s1 * p2,
            p1 * t2 + t1 * p2 + r1 * s2 + s1 * r2,
        )

    def star(self):
        # Order-2 jet star over two independent infinitesimals ε, δ (both nilpotent
        # with ε² = δ² = 0), with εδ as the mixed coefficient. Non-commutative in P.
        #
        # For x = p + rε + sδ + tεδ:
        #   p-component: p*
        #   r-component: p* r p*
        #   s-component: p* s p*
        #   t-component: p* r p* s p* + p* s p* r p* + p* t p*
        #
        # The two cross-terms arise from interleaving rε and sδ in different orders
        # within x^n; they coincide only if P (and the P-action on the s/t vectors)
        # is commutative.
        ps = self.p.star()
        ps_r_ps = ps * self.r * ps
        ps_s_ps = ps * self.s * ps
        t_sum = (ps_r_ps * self.s * ps) + (ps_s_ps * self.r * ps) + (ps * self.t * ps)
        return SecondOrderExpectation(ps, ps_r_ps, ps_s_ps, t_sum)

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
