"""
Uncertainty-set meta-semiring (work-in-progress; was test_uncertainty_sets).

`UncertaintySet(base)` lifts a base semiring pointwise to finite sets of its
values under Minkowski addition and multiplication:

    X + Y = {x + y : x in X, y in Y}    (using base's +)
    X * Y = {x * y : x in X, y in Y}    (using base's *)
    zero  = {base.zero}
    one   = {base.one}

Sub-distributive rather than distributive:

    (X + Y) * Z ⊆ X*Z + Y*Z

    strict example: X = {1}, Y = {1}, Z = {0, 1}
      (X + Y) * Z = {0, 2}
      X*Z + Y*Z   = {0, 1, 2}

The `1` on the right comes from picking z=0 in the first summand and z=1 in
the second — independence the left side can't match. Same structural reason
Interval is sub-distributive.

Kept in experimental/ rather than promoted to the library because it fails
distributivity (the axiom suite rejects it), matching the convention we use
for semilinear/funky.
"""

from semirings.base import Semiring


def UncertaintySet(base):
    """Meta-semiring: lift `base` pointwise to finite sets under Minkowski ops."""

    class U(Semiring):
        def __init__(self, xs):
            self.xs = frozenset(xs)

        def __add__(self, other):
            return U({x + y for x in self.xs for y in other.xs})

        def __mul__(self, other):
            return U({x * y for x in self.xs for y in other.xs})

        def __eq__(self, other):
            return isinstance(other, U) and self.xs == other.xs

        def __hash__(self):
            return hash(self.xs)

        def __repr__(self):
            return f'U({set(self.xs)})'

    U.zero = U({base.zero})
    U.one = U({base.one})
    U.__name__ = f'U[{base.__name__}]'
    return U


if __name__ == '__main__':
    from semirings import Float
    U = UncertaintySet(Float)
    x, y, z = U({1}), U({1}), U({0, 1})
    print(f'(x+y)*z = {(x + y) * z}')
    print(f'x*z+y*z = {x * z + y * z}')
    assert ((x + y) * z).xs <= (x * z + y * z).xs
    assert (x + y) * z != x * z + y * z
