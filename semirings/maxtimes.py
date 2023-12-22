from . import base


class MaxTimes(base.Semiring):

    def __init__(self, score, d=None):
        self.score = score
        self.d = d

    def __repr__(self):
        return f'MaxTimes({self.score}, {self.d})'

    def __eq__(self, other):
        return (isinstance(other, MaxTimes)
                and (self.score == other.score   # for infinity
                     or abs(self.score - other.score) < 1e-10))

    def __hash__(self):
        return hash(self.score)

    def __lt__(self, other):
        return isinstance(other, MaxTimes) and self.score < other.score

    def __add__(self, other):
        return max(self, other)

    def __mul__(self, other):
        if other is one: return self
        if self is one: return other
        if self is zero: return zero
        if other is zero: return zero
        return MaxTimes(self.score * other.score, (self.d, other.d))

    @classmethod
    def multiplicity(cls,x,m):
        if m > 0:
            return x
        else:
            return cls.zero

    def star(self):
        return self.one + self if self.score <= 1 else MaxTimes.inf


MaxTimes.zero = zero = MaxTimes(0, None)
MaxTimes.one = one = MaxTimes(1, ())
MaxTimes.inf = MaxTimes(float('+inf'), None)
