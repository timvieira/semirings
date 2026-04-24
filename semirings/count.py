from semirings.base import Semiring


class Count(Semiring):
    """Natural-number counting semiring (ℕ, +, ·, 0, 1).

    Counts how many ways something can happen. Not closed under star —
    star(0) = 1 but star(n) diverges for n ≥ 1.
    """

    def __init__(self, x):
        assert isinstance(x, int)
        self.x = x

    def __add__(self, other): return Count(self.x + other.x)
    def __mul__(self, other): return Count(self.x * other.x)
    def __sub__(self, other): return Count(self.x - other.x)
    def __eq__(self, other):  return isinstance(other, Count) and self.x == other.x
    def __hash__(self):       return hash(self.x)
    def __lt__(self, other):  return self.x < other.x
    def __repr__(self):       return f'Count({self.x})'

    def lower(self):  return self.x
    __float__ = lower

    def metric(self, other):
        d = abs(self.x - other.x)
        s = max(1, abs(self.x), abs(other.x))
        return d / s


Count.zero = Count(0)
Count.one = Count(1)
