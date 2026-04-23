"""
Semilinear-set semiring sketch (work-in-progress; was test_semilinear in tests/test_semirings.py).

The values are Counters keyed by count-vectors in N^k (one coordinate per alphabet
symbol). Addition is bag-union; multiplication is the Minkowski sum on keys, with
multiplicities multiplied (see the `# or should it be sum?` comment below for an
unresolved modeling choice). A one-hot lift embeds a single symbol as the basis
vector e_i with multiplicity 1.

This is roughly the commutative image / Parikh-semiring construction: it tracks
how many times each terminal appears, collapsing non-commutative structure. Useful
as a base for Parikh's theorem and semilinear-set reasoning in formal language theory.

Why this is not a test:

1. `star` is undefined. For this semiring, star(x) would be the infinite bag
   {0*x, 1*x, 2*x, ...} — an unbounded support, so the naive `star_fixpoint`
   never converges and `star_approx(x, T)` only ever sees the first T terms.
   A principled star would need a finite symbolic representation of semilinear
   sets (i.e. unions of linear sets a + Nb_1 + ... + Nb_m), not a Counter.

2. The former test file had `assert approx == analytical, [approx, analytical]`
   commented out, alongside `analytical = S.star(x)` — both referencing a star
   that does not exist on this construction. The `if approx.score > 1000:
   approx = S.inf` guard was also cargo-culted from scalar star-approx tests;
   Counter values have no `.score`.

3. The `# or should it be sum?` comment in `minkowski` flags a real modeling
   question: for Parikh-style counting, repeated factorizations should *add*
   (count multiplicities), but the current code *multiplies*. Resolving this
   is prerequisite to writing meaningful tests.

4. `one` is constructed as `Counter({emptyvector: 0})` — multiplicity **zero**.
   Under `minkowski`, that makes `one * x` evaluate to a Counter of keys-from-x
   all with multiplicity 0, i.e. it annihilates. The probable intent was
   `Counter({emptyvector: 1})`. This bug silently propagated through the
   commented-out test because nothing ever asserted `one * x == x`.

TODO before promoting to a real semiring:
  - Fix `one` (finding 4) — currently annihilates everything.
  - Decide multiplication semantics on repeated keys (add vs. multiply).
  - Represent semilinear sets symbolically so that star has a closed form.
  - Add an equality/metric that respects that representation.
"""

from collections import Counter
from semirings import make_semiring
from arsenal import Integerizer


def bag_union(x, y):
    z = Counter()
    for k, v in x.items(): z[k] += v
    for k, v in y.items(): z[k] += v
    return z


def minkowski(x, y):
    z = Counter()
    for k1, v1 in x.items():
        for k2, v2 in y.items():
            assert len(k1) == len(k2)
            k3 = tuple(a + b for a, b in zip(k1, k2))
            z[k3] += v1 * v2   # or should it be sum?
    return z


def build(symbols='abcd'):
    alphabet = Integerizer(list(symbols))
    emptyvector = (0,) * len(symbols)

    S = make_semiring(
        'Semilinear',
        bag_union,
        minkowski,
        Counter(),                          # zero: empty bag
        Counter({emptyvector: 0}),          # one: ??? — see notes above
    )

    def onehot(symbol):
        x = [0] * len(symbols)
        x[alphabet(symbol)] = 1
        return S({tuple(x): 1})

    return S, onehot


if __name__ == '__main__':
    S, onehot = build('abcd')
    a, b, c, d = map(onehot, 'abcd')

    print(a, b, c, d)
    print(a + a*a + a*a*a + a*a*a + a*b*b + b*a*b + b*b*a)

    # star_approx only sees the first 10 geometric-series terms; the full star
    # would be the infinite semilinear set {n*a : n >= 0}.
    print('star_approx(a, 10) =', S.star_approx(a, 10))
