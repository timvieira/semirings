import math

import numpy as np
from . import base


class Float(base.Semiring):
    def __init__(self):
        assert False, 'should never be called'

    zero = 0.0
    one = 1.0

    def approx_zero(x):
        return abs(x) < 1e-7   # XXX: YUCKY MAGIC NUMBER

    @classmethod
    def lift(cls, x, _):
        return x

    @staticmethod
    def star(x):
        if x == 1: return np.inf
        return 1/(1-x)

    @classmethod
    def multiplicity(cls,x,m):
        return x*m

    def metric(self, other):
        # Equal values are distance 0, including ±inf. Otherwise, the standard
        # formula `|a - b| / max(1, |a|, |b|)` produces nan when either
        # argument is infinite (∞ - ∞ = nan), which fails downstream
        # convergence checks like `metric <= 1e-8`. Special-case the infinite
        # branch to return 1.0 — non-zero, comparable, and consistent with
        # the formula's normalized range.
        if self == other: return 0.0
        if math.isinf(self) or math.isinf(other): return 1.0
        return abs(self - other) / max(1, abs(self), abs(other))
