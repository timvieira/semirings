from . import base


class MaxPlus(base.Semiring):

    def __init__(self, score, d=None):
        self.score = score
        self.d = d

    def __repr__(self):
        return f'MaxPlus({self.score}, {self.d})'

    def __eq__(self, other):
        return (isinstance(other, MaxPlus)
                and (self.score == other.score
                     or abs(self.score - other.score) < 1e-10))

    def __hash__(self):
        return hash(self.score)

    def __lt__(self, other):
        return isinstance(other, MaxPlus) and self.score < other.score

    def __add__(self, other):
        return max(self, other)

    def __mul__(self, other):
        return MaxPlus(self.score + other.score, [self.d, other.d])

    @classmethod
    def multiplicity(cls,x,m):
        if m > 0:
            return x
        else:
            return cls.zero

    def star(self):
        # a positive reward cycle goes to infinity,
        # everything else stays put x^* = max(0, x, x+x, x+x+x, ...)
        return (self.one + self) if self.score <= 0 else MaxPlus.inf


MaxPlus.zero = zero = MaxPlus(float('-inf'), None)
MaxPlus.one = one = MaxPlus(0.0, ())
MaxPlus.inf = MaxPlus(float('+inf'))
