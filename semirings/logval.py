import numpy as np
from numpy import log as _log, exp, isnan, log1p as _log1p, expm1
from arsenal import colors
from . import base


def log(x):
    if x <= 0:
        return float('-inf')
    return _log(x)


def log1p(x):
    if x <= -1:
        return float('-inf')
    return _log1p(x)


def log1pexp(x):
    """
    Numerically stable implementation of log(1+exp(x)) aka softmax(0,x).

    -log1pexp(-x) is log(sigmoid(x))

    Source:
    http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    """
    if x <= -37:
        return exp(x)
    elif -37 <= x <= 18:
        return log1p(exp(x))
    elif 18 < x <= 33.3:
        return x + exp(-x)
    else:
        return x


def log1mexp(x):
    """
    Numerically stable implementation of log(1-exp(x))

    Source:
    http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    """
    assert x <= 0
    a = abs(x)
    if 0 < a <= 0.693:
        return log(-expm1(-a))
    else:
        return log1p(-exp(-a))


class LogVal(base.Semiring):

    def __init__(self, pos, ell):
        self.pos = pos
        self.ell = ell

    @classmethod
    def lift(cls, x):
        return cls(x >= 0, log(abs(x)))

    def __eq__(self, other):
        return np.allclose(float(self), float(other))

    def __hash__(self):
        return hash(float(self))

    def __lt__(self, other):
        return float(self) < float(other)

    def is_zero(self):
        return self.ell <= float('-inf')

    def __float__(self):
#        if self.is_zero():
#            return 0.0
        if self.pos:
            return +float(exp(self.ell))
        else:
            return -float(exp(self.ell))

    def __mul__(self, b):
        c = LogVal.lift(0.0)
#        if self.is_zero() or b.is_zero():
#            return c
        c.pos = self.pos == b.pos
        c.ell = self.ell + b.ell
        return c

    def __truediv__(self, b):
        c = LogVal.lift(0.0)
        if self.is_zero():
            return c
#        if b.is_zero():
#            return c         # divide by zero: Should probably return NaN.
        c.pos = self.pos == b.pos
        c.ell = self.ell - b.ell
        return c

    __div__ = __truediv__

    def __sub__(self, b):
        c = LogVal.lift(0.0)
        c.ell = b.ell
        c.pos = not b.pos
        return self + c

    def __add__(self, b):
        c = LogVal.lift(0.0)
        a = self
        if a.ell < b.ell:
            a, b = b, a
        if a.is_zero():
            c.pos = b.pos
            c.ell = b.ell
            return c
        if b.is_zero():
            c.pos = a.pos
            c.ell = a.ell
            return c
        x = b.ell - a.ell
        assert not isnan(x)
#        assert x <= 0         # unnecessary assertion.
        if a.pos == 1 and b.pos == 1:
            c.pos = 1
            c.ell = a.ell + log1pexp(x)   # log(1+exp(x))
        elif a.pos == 1 and b.pos == 0:
            c.pos = 1
            c.ell = a.ell + log1mexp(x)   # log(1-exp(x))
        elif a.pos == 0 and b.pos == 1:
            c.pos = 0
            c.ell = a.ell + log1mexp(x)
        else:
            c.pos = 0
            c.ell = a.ell + log1pexp(x)
        return c

    def star(self):
        return LogVal.one / (LogVal.one - self)

    def multiplicity(self, m):
        return self.lift(m) * self

    def __repr__(self):
        x = float(self)
        sign = '+' if self.pos else colors.cyan % '-'
        return f'LogVal({x:.2g}={sign}exp({self.ell:.2g}))'

LogVal.zero = LogVal.lift(0.0)
LogVal.one = LogVal.lift(1.0)
