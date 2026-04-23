"""
Endomorphism semiring sketch (work-in-progress; was test_endomorphism).

Elements are endomorphisms f: M → M of some monoid (M, +), with:

    (f + g)(x) = f(x) + g(x)      (pointwise sum in M)
    (f * g)(x) = f(g(x))          (composition)
    zero       = x ↦ 0
    one        = x ↦ x

For this to be a semiring the elements must actually be *endomorphisms* —
they have to preserve the monoid structure:

    f(x + y) = f(x) + f(y)

If M is (ℝ, +), those are exactly the linear (not affine!) functions
f(x) = cx. The sketch in the test file used members like `lambda x: 2*x + 1`
and constant functions `lambda _: 1` — which are affine, not linear, so they
are NOT endomorphisms of (ℝ, +). That's why the sketch only tested two
hand-picked distributivity instances and left the full axiom suite commented
out.

TODO before promoting to a real semiring:
  - Restrict members to honest endomorphisms (linear maps, or something with
    a checkable monoid-homomorphism property).
  - Equality and hashing: functions don't have extensional equality in Python;
    need a finite canonical representation (matrix? symbolic form?).
  - Decide what `check_multiplicity` and `star` mean here.
"""

from semirings import make_semiring


def build():
    return make_semiring(
        'Endo',
        lambda f, g: lambda x: f(x) + g(x),
        lambda f, g: lambda x: f(g(x)),
        lambda _: 0,
        lambda x: x,
    )


if __name__ == '__main__':
    S = build()
    # NOTE: these are affine, not endomorphisms of (ℝ, +). Distributivity
    # happens to hold for the specific triples checked below, not in general.
    members = list(map(S, [
        lambda x: 2*x + 1,
        lambda x: x * 2,
        lambda x: x * 3,
        lambda x: 4 * x,
        lambda x: 1 * x,
        lambda _: 1,
        lambda _: 3,
    ]))
    a, b, c, d, e, f, g = members

    def eval_(f):
        return [f.x(x) for x in range(-5, 5)]

    assert eval_((a + b) * c) == eval_(a * c + b * c)
    assert eval_(c * (a + b)) == eval_(c * a + c * b)
    print('OK (affine sketch; not a proper endomorphism semiring)')
