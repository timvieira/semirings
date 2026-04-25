import math

from semirings.base import Semiring


class Count(Semiring):
    """Counting semiring (ℕ ∪ {∞}, +, ·, 0, 1).

    Counts how many ways something can happen. Closed under star:
    star(0) = 1, star(n) = ∞ for n ≥ 1. The saturating element is
    `math.inf`, with the convention 0·∞ = 0.
    """

    def __init__(self, x):
        assert isinstance(x, int) or x == math.inf, x
        self.x = x

    def __add__(self, other): return Count(self.x + other.x)
    def __mul__(self, other):
        # 0·∞ = 0 by convention; float multiplication would give nan
        if self.x == 0 or other.x == 0: return Count.zero
        return Count(self.x * other.x)
    def __sub__(self, other): return Count(self.x - other.x)
    def __eq__(self, other):  return isinstance(other, Count) and self.x == other.x
    def __hash__(self):       return hash(self.x)
    def __lt__(self, other):  return self.x < other.x
    def __repr__(self):       return f'Count({self.x})'

    def star(self):
        return Count.one if self.x == 0 else Count(math.inf)

    def lower(self):    return self.x
    def __float__(self): return float(self.x)

    def metric(self, other):
        if self.x == other.x: return 0.0
        if self.x == math.inf or other.x == math.inf: return 1.0
        d = abs(self.x - other.x)
        s = max(1, abs(self.x), abs(other.x))
        return d / s


Count.zero = Count(0)
Count.one = Count(1)


def make_count_k(K):
    """k-truncated counting semiring: ℕ saturated at K. Closed under star
    (star(0) = 1, star(n>0) = K), unlike plain Count."""
    from semirings.base import make_semiring
    return make_semiring(
        f'Count[{K}]',
        lambda a, b: min(K, a + b),
        lambda a, b: min(K, a * b),
        0,
        1,
        star=lambda x: 1 if x == 0 else K,
    )
