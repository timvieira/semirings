from . import base


# WARNING: Intervals are not a semiring.
def make_interval(*, bot, zero, top):

    class Interval(base.Semiring):
        def __init__(self, lo, hi):
            self.lo = lo
            self.hi = hi
            self.empty = lo > hi

        def __repr__(self):
            if self.lo > self.hi: return '∅'
            if self.lo == self.hi: return f'[{self.lo}]'
            return f'[{self.lo},{self.hi}]'

        def __hash__(self):
            return hash((self.lo, self.hi))

        def __eq__(self, other):
            return self.empty == other.empty and self.lo == other.lo and self.hi == other.hi

#        def __lt__(self, other):
#            "definitely less than"
#            # [lo, hi] < [lo', hi']
#            #return self.hi < other.lo
#            raise NotImplementedError()

        def __add__(self, other):
            if self.empty: return empty
            if other.empty: return empty
            return Interval(self.lo + other.lo,
                            self.hi + other.hi)
#            a = self.lo + other.lo
#            b = self.lo + other.hi
#            c = self.hi + other.lo
#            d = self.hi + other.hi
#            return Interval(min([a,b,c,d]),
#                            max([a,b,c,d]))

        def __mul__(self, other):
            if self.empty: return empty
            if other.empty: return empty
            a = self.lo * other.lo
            b = self.lo * other.hi
            c = self.hi * other.lo
            d = self.hi * other.hi
            return Interval(min([a,b,c,d]),
                            max([a,b,c,d]))

        def inv(self):
            if self.lo > self.hi:                  return self   # no-op on an empty interval
            elif self.hi == 0:                     return Interval(bot, 1/self.lo)
            elif self.lo == 0:                     return Interval(1/self.hi, top)
            elif (self.lo <= zero <= self.hi):     return Interval.loose
            else:                                  return Interval(1/self.hi, 1/self.lo)  # if zero is not in interval

        def __truediv__(self, other):
            return self * other.inv()

        def __neg__(self):
            if self.empty: return empty
            # -self swaps lo and hi because it *reflects* across 0.
            return Interval(-self.hi, -self.lo)

        def __sub__(self, other):
            return self + -other

        def __or__(self, other):
            "(lossy) union of intervals"
            return Interval(min(self.lo, other.lo), max(self.hi, other.hi))

        def __and__(self, other):
            "intersection of intervals"
            return Interval(max(self.lo, other.lo), min(self.hi, other.hi))

        def __contains__(self, x):
            return not self.empty and self.lo <= x.lo <= x.hi <= self.hi

        @classmethod
        def sum(cls, xs):
            y = cls.zero
            for x in xs:
                y += x
            return y

        @classmethod
        def tightest(cls, xs):
            xs = list(xs)
            if len(xs) == 0: return cls.loose
            return cls(max(x.lo for x in xs), min(x.hi for x in xs))

        @classmethod
        def union(cls, xs):
            xs = list(xs)
            if len(xs) == 0: return cls.loose
            return cls(min(x.lo for x in xs), max(x.hi for x in xs))

        def star(self):
            return Interval(1,1)/(Interval(1,1) - self)

        def multiplicity(self, m):
            return Interval(m, m) * self

    empty = Interval(top, bot)
    Interval.zero = Interval(zero, zero)
    Interval.loose = Interval(bot, top)
    return Interval


Interval = make_interval(bot=float('-inf'), zero=0, top=float('+inf'))
Interval.one = Interval(1,1)
