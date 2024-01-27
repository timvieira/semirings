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
        return abs(self - other)
