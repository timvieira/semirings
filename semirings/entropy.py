import numpy as np
from . import base


class Entropy(base.Semiring):
    """Entropy semiring for computing the entropy of a distribution.

    Elements are pairs (p, r) where p is a probability and r = sum of p*log(p)
    terms.  The entropy H = log(p) - r/p can be read off after summation.

    Addition: (p1, r1) + (p2, r2) = (p1+p2, r1+r2)
    Multiplication: (p1, r1) * (p2, r2) = (p1*p2, p1*r2 + r1*p2)
    """

    def __init__(self, p, r):
        self.p = p
        self.r = r

    def __repr__(self):
        return f'Entropy({self.p}, {self.r})'

    def __eq__(self, other):
        return isinstance(other, Entropy) and self.p == other.p and self.r == other.r

    def __hash__(self):
        return hash((self.p, self.r))

    def __add__(self, other):
        if other is zero: return self
        if self is zero: return other
        return Entropy(self.p + other.p, self.r + other.r)

    def __mul__(self, other):
        if other is one: return self
        if self is one: return other
        if other is zero: return zero
        if self is zero: return zero
        return Entropy(
            self.p * other.p,
            self.p * other.r + self.r * other.p,
        )

    def star(self):
        if self.p == 1: return Entropy(float('inf'), float('inf'))
        ps = 1 / (1 - self.p)
        return Entropy(ps, ps * ps * self.r)

    @classmethod
    def lift(cls, x):
        xlogx = x * np.log2(x) if x != 0 else 0
        return cls(x, xlogx)

    @property
    def H(self):
        if self.p == 0: return 0.0    # lim_{p→0+} p·log(p) = 0
        return np.log2(self.p) - self.r / self.p

    def metric(self, other):
        "Technically, a pseudo metric"
        return abs(self.H - other.H)


Entropy.zero = zero = Entropy(0.0, 0.0)
Entropy.one = one = Entropy(1.0, 0.0)
