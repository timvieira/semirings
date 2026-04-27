# TODO

## Misc
- retire `approx_zero`

## Audit metrics
- check that each semiring has decent metric implementation

## Reduce hard dependencies
- Inline `arsenal` utilities (`log1pexp`/`log1mexp` in `logval.py`, `sorted_union`/`sorted_product` in `lazysort.py`)—they're small enough to vendor
- Make `regex.py` optional (lazy import) so the core package doesn't require `fsa`/`wfsa`
- Make `mert.py` optional so `scipy` isn't required at install time

## Clean up `__init__.py`
- Replace `from semirings.misc import *` with explicit imports
- Add `__all__` to `__init__.py` and all submodules

## Make `Dual` inherit from `Semiring`
- The `dual()` factory in `misc.py` creates a class that doesn't subclass `Semiring`, so it lacks `chart()`, `sum()`, `product()`, `__pow__`, `star_approx`, etc.

## Remove commented-out code
- `base.py`: commented-out `__lt__`, `Chart.__add__`, `Chart.__sub__`, `Chart.round`
- `misc.py`: commented-out `WrappedBackpointers` class

## Extract shared backpointer logic
- MaxPlus, MinPlus, MaxTimes, MinTimes all duplicate nearly identical backpointer tracking (`d` argument). Extract into a reusable base class.

## Fix `multiplicity` shadowing in `minmax`/`maxmin`
- `minmax` and `maxmin` in `misc.py` define `multiplicity` as an instance method, shadowing the classmethod inherited from `Semiring`.

## Log-space Entropy semiring
- The current `Entropy` semiring operates in probability space `(p, r)`, which suffers from underflow for small probabilities. Implement a version that stores `(log p, log r)` (or similar) and does all arithmetic in log space, analogous to how `LogVal` relates to `Float`.

## Other semirings sitting around in other projects
- Centralize my collection of semirings that are currently scattered across many projects.

## Fill out the GKT 2007 provenance hierarchy

We have `Boolean`, `Count`, `Why`, `Lineage`, and `Bag`, but the universal element of the
[Green-Karvounarakis-Tannen 2007](https://scholar.google.com/scholar?q=Green+Karvounarakis+Tannen+2007+Provenance+Semirings)
hierarchy — **`N[X]`**, the polynomial provenance over tuple-IDs — is missing. `N[X]` is
the *initial* commutative semiring over `X`: every other provenance in the hierarchy
(`B`, `N`, `Why`, `Lin`, `PosBool`, `B[X]`, `Trio`, security/trust lattice) is a
homomorphic image of it. Adding it unlocks the "parse once, interpret many" demo.

**Concrete gaps:**

- **`N[X]`** — ℕ-coefficient polynomials over tuple-IDs with monomial collection.
  This is the headline addition.
- **`Bag` mismatch.** `bag.py:36-37` tags products as pairs `(x, y)` rather than
  collecting monomials, so `{a:1} * {a:1} = {(a,a):1}` instead of `{a²:1}`. It's a
  free-bag-of-traces, not GKT's `N[X]`. Decide: fix it to be `N[X]`, rename it
  (e.g. `Trace` / `FreeBag`) and add `N[X]` separately, or retire it once `N[X]` lands.
- **`PosBool[X]`** — positive Boolean expressions over tuple-IDs (idempotent
  ⊕ and ⊗ with absorption). Natural target for "query without multiplicity."
- **Named security/trust lattice.** The `minmax` factory already produces it, but
  the trust/clearance use case is invisible from the README — wrap it as a named
  class (e.g. `SecurityLattice(levels)` / `Trust(levels)`).
- **`B[X]`** and **`Trio(X)`** — lower priority; add only if a concrete demo wants them.

**Demo target (the EXPLORE.md provenance entry's "first question"):** one notebook
that builds a derivation in `N[X]` and projects it homomorphically to `B`, `N`,
`Why`, `Lin`, `PosBool`, and a trust lattice. Doing the demo will also force the
`Bag`-vs-`N[X]` decision above.

## Deferred: `forget_multiplicity` / support functor

The homomorphism `FreeClosedSemiring → RegularLanguageSet` (drop ℕ-counts to
{0,1}, recover the underlying regular language). Drafted and scrutinized in
`experimental/free_semirings.py`; see the long `DEFERRED:` comment block there
for the full discussion. Two issues blocked shipping:

1. **ε-arcs** introduced by wfsa concatenation/star require explicit handling;
   structural mirror leaks them as literal symbols.
2. **Path cancellation** — the structural mirror is sound only over zerosumfree,
   zero-divisor-free weight semirings (a "positive cone"). Current `Count` is
   really ℤ-Count (it has `__sub__`), not ℕ-Count, so even our motivating case
   isn't safe.

Pickup options (in order of likely value): add a `CountN` carrier and ship the
structural mirror against it; implement via Hankel-rank zero-test for general
`R`; or drop it (users can use `RegularLanguageSet` / `fsa.FSA` directly for
language-level work). Connects to "Audit carrier sets and star conventions"
below — `Count`-vs-`CountN` is the same flavor of carrier-set ambiguity.

## Algebra cherry-picks and redirects (DESIGN.md Stage 4 — deferred)

The `~/projects/algebra/` repo still has its own implementations of every concrete semiring and some unique types. Stage 0 added `semirings` as a dep in `algebra/setup.py`, but no redirects have happened.

**Scope:**

- Duplicate semirings in `algebra/algebra/number/`: MinPlus, MaxPlus, Boolean, LogVal, Bottleneck (in `fuzzy.py`), Regex, Interval. Replace each with a re-export from `semirings`.
- Unique-to-algebra types worth porting: `dual.py` (redundant once Stage 3 Jet lands — becomes `Jet(S, S, order=1)` alias), `poly.py` / `Monomial` (port if needed), `baggr.py` / `BagRing` (port if the DB-provenance use case matters; otherwise skip).
- **Stays in `algebra/`:** `Ring`, `Field`, `ClosedSemiring`, and all the matrix algorithms (`kleene.py`, `lu.py`, `gaussjordan.py`, `bareiss.py`, `eig.py`, `iterative.py`, `wfsa.py`, `graph.py`). These are `algebra/`'s purpose — they're parameterized over semirings.

**Per-item protocol:** same as Stages 1–2. Port to `semirings/` (or confirm it's already there), add AXIOM_CASES entry if a real semiring, replace the algebra-side file with a shim or re-export, run algebra's test suite, verify no regressions beyond baseline.

**Baseline (recorded during Stage 0):** `/tmp/baseline-algebra.txt` — 3 passed, 5 failed, 9 errors. Most failures are pre-existing numpy-version drift (`numpy.product` removed) or missing external deps (`spanning_tree`, `algebra.mtt.bf`). When evaluating the redirect, diff against this baseline — don't expect a green board.

**Ordering tip:** most value is in the concrete-semiring redirects (one source of truth). `dual.py` becomes trivial after Stage 3. `poly.py`/`baggr.py` are lowest priority.

## Pairing framework (DESIGN.md Stage 3 — deferred)

Full design is in `DESIGN.md` Part I. Summary so that future sessions can pick up:

**What:** two semiring-construction combinators plus a free-semiring realization.

- `Product(S, D)` — categorical product. Componentwise +, *, 0, 1. Requires D to be a semiring.
- `Jet(S, M, order=k)` — truncated polynomial S[ε]/(ε^(k+1)). First-order: `(p, r)` with `(p1, r1)·(p2, r2) = (p1·p2, p1·r2 + r1·p2)` (note the non-commutative order — see below). M defaults to S and must be an S-module otherwise. k=1 and k=2 are the cases we need.
- `SelectedDerivation(X)` — chain-of-generators monoid promoted to a semiring via sum-selection. Paired with a selective S (MaxPlus/MinPlus) gives classical backpointers.
- `FreeClosedSemiring(X)` — free closed semiring over X, backed by the existing `wfsa` package with ℕ weights. Decidable equality via WFSA minimization. **Open question: does `wfsa` already support ℕ weights? If not, upstream PR first.**

**Why:** collapses three ad-hoc classes into Jet aliases.

- `Dual(S)` (misc.py)            = `Jet(S, S, order=1)`
- `Expectation` (expectation.py) = `Jet(LogVal, LogVal, order=1)`
- `SecondOrderExpectation`       = `Jet(LogVal, LogValVector, order=2)` (two independent infinitesimals — may need a two-variable variant)
- Plus `Product(S, FreeClosedSemiring)` and `Product(S, SelectedDerivation)` as the general-purpose derivation trackers, replacing hypergraphs' per-class `.d` field pattern.

**Design decisions already settled (from session discussion):**

1. The jet multiplication uses `p1·r2 + r1·p2`, not `p1·r2 + p2·r1`. The session's Stage 2 commit (`a3ad760`) fixed this latent bug in the existing `Expectation` / `SecondOrderExpectation`, and an `ExpectationRL = make_expectation(RegularLanguage, RegularLanguage)` entry was added to AXIOM_CASES to guard against regression.
2. First-order star: `star((p, r)) = (p*, p* · r · p*)`. Non-commutative.
3. Second-order star: `star((p, r, s, t)) = (p*, p* r p*, p* s p*, p* r p* s p* + p* s p* r p* + p* t p*)`. The two cross-terms only collapse if the base is commutative.
4. LogVal's `__mul__` now honors absorbing zero — required so `0 · ∞ = 0` after star introduces infinities (was previously NaN).

**Retrofitting plan (Stage 3 execution):**

- Implement `Product`, `Jet`, `SelectedDerivation`, `FreeClosedSemiring` in `semirings/pairing.py` and `semirings/free.py` (the latter already exists with `FreeExpr` — add `FreeClosedSemiring` there too).
- Replace `Dual`, `Expectation`, `SecondOrderExpectation` implementations with `Jet(...)` aliases. Keep public names.
- AXIOM_CASES: no changes needed — existing entries continue to exercise the same semantics via the aliases.
- Verification: all existing tests must still pass.

**Effort estimate:** 1–2 days, per DESIGN.md Part III Stage 3.

## Audit carrier sets and star conventions per semiring

Several semirings currently have star behavior that is technically correct but
undocumented about *which* algebraic structure is in play. The carrier set
determines whether star is the Conway-axiom solution of `x* = 1 + x·x*` or
the ω-continuous supremum `sup_n (1 + x + x² + … + x^n)`, and these disagree
for many values.

Concrete examples to untangle:

- **`Float.star(2)` returns `-1`.** That's the Conway answer treating `Float`
  as the field ℝ — `-1 = 1 + 2·(-1)` solves the fixed-point equation. The
  ω-continuous answer over `[0, +∞]` would be `+∞` (sup of `1, 3, 7, 15, …`).
  Different structure, different value; both are valid in their own carrier.
- **`Float.star(-1)` returns `0.5`.** Conway on ℝ: `0.5 = 1 + (-1)·0.5`, fine.
  ω-continuous: the partial sums `1, 0, 1, 0, …` oscillate — no supremum, not
  closed. So `Float` as currently written is implicitly the Conway / field-ℝ
  semiring, not the ω-continuous nonneg-reals one.
- **`Float.star(1)` returns `+∞`** — a special case that jumps to the
  ω-continuous answer since Conway has no ℝ-solution here. Inconsistent with
  the rest of `Float.star`.
- **`MinPlus` over extended reals `ℝ ∪ {+∞, -∞}`** is closed in both senses
  and they agree: `star(x) = 0` for `x ≥ 0`, `-∞` for `x < 0`. No ambiguity,
  but the carrier should be stated.
- **`ConvexHull.star` is `NotImplementedError`** because general hulls aren't
  closed under either interpretation. A cones-only variant or a Pareto-frontier
  variant (`experimental/pareto.py`) would be closed — which one we offer
  should be a deliberate choice, not the absence of one.

What to do:

1. For each semiring, state the intended carrier in the class docstring
   (e.g., "Float over ℝ as a Conway semiring" vs "Float over `[0, +∞]` as an
   ω-continuous semiring"), and pick one. If both are useful, expose them as
   separate classes (`Float`, `FloatNonneg`) rather than overloading one name.
2. Name the star convention explicitly: `star_conway`, `star_omega`, or
   whatever pairs well with the existing `star_approx` / `star_fixpoint` /
   `star_doubling` helpers. Let `.star()` dispatch to the class's chosen one.
3. Resolve the `Float.star(1) = inf` special case: it's the ω-continuous
   answer embedded in an otherwise-Conway implementation. Either return `nan`
   (Conway has no ℝ-solution), or make `Float` explicitly ω-continuous on
   `[0, +∞]` and use the Conway reals as a separate class.
4. AXIOM_CASES should exercise star at the boundary values (`x = 1`, `x = -1`,
   `x = ±∞`, `x = 2` for Float; zero / one / negative / `+∞` for MinPlus) so
   the chosen convention is pinned by tests.
5. **Catalog ω-continuity per semiring.** Once the carriers are pinned down,
   mark each class with its structural properties — ω-continuous, Conway,
   idempotent, commutative, complete — as attributes (`is_omega_continuous`,
   `is_conway`, etc.) or a single `structural_properties` set. This is
   pedagogically useful (README / compendium can surface which semirings sum
   infinite series vs solve fixpoint equations) and also enables tests that
   check, e.g., "if a semiring claims ω-continuity, then `star_approx(x, T)`
   converges monotonically to `star(x)` as T grows." Likely labels to track:
   `ω-continuous` (MinPlus, MaxPlus, Boolean, Count, MinTimes on `[0, 1]`,
   MaxTimes on `[0, 1]`; **note:** `LogVal` is the *signed* log-semiring, so
   it's Conway-style on ℝ rather than ω-continuous on `[0, +∞]`),
   `Conway-only` (Float over ℝ, LogVal over ℝ,
   anything we decide is a field rather than a nonneg-semiring), `neither`
   (ConvexHull on general hulls, Interval), `idempotent` (MinPlus, MaxPlus,
   Boolean, CutSets, Why, Lineage, Bottleneck), `commutative` (most), and so
   on.

Related: "Clarify the `star` convention" in the code-quality audit below
touches the same area but from the algorithm-selection angle (fixpoint vs
doubling vs closed form). This item is about the mathematical structure; pick
that first, then the algorithm follows.

## Axiom-test coverage for FreeExpr and sampling semirings
- `FreeExpr` is a magma — syntactic tree equality fails commutativity/associativity
  of `+` and `*` — but it is equal-up-to-semiring-axioms under the evaluation
  homomorphism `weight(·) : FreeExpr → ℝ`. Add a `metric` to `FreeExpr` using
  `abs(weight(self) - weight(other))` so axiom-test `assert_equal` succeeds via
  the tolerance path. First fix `weight()` to handle `FreeExpr.zero` / `FreeExpr.one`
  (their `.args` is empty) and to return `+inf` when `Star(x)` with `weight(x) >= 1`.
- Sampling semirings (`Expon`, lazy `Sample`, SWOR `Sample`) are `FreeExpr + an
  evaluation pass`. They already carry the accumulated weight in `self.w`. Add
  `metric` using `abs(self.w - other.w)` and add to AXIOM_CASES with appropriate
  kwargs (eager/SWOR: `star=False`; lazy: full). Use `w` in `[0, 1)` for samples
  so `Star` doesn't diverge.

## Fix `Float.__init__`
- `Float.__init__()` takes no arguments, so `Float(3)` raises `TypeError`. The README examples don't work.

## Better implementations of chart
- I think I have a better implementation of a chart

## Interval is not a semiring
- `Interval` inherits from `Semiring` but violates distributivity (it's sub-distributive: `a * (b + c) ⊆ a*b + a*c`, not equal). The test (`todo_interval`) is disabled for this reason.
- Options:
  1. **Keep as-is, document clearly.** Interval arithmetic gives sound over-approximations in semiring-style algorithms. Just improve the docs and re-enable the test with sub-distributivity-aware assertions.
  2. **Stop inheriting from `Semiring`.** `Interval` already overrides most methods; it only inherits `chart()`, `product()`, `__pow__`, and `star_approx`/`star_fixpoint` (and overrides `star`).
  3. **Introduce a lighter base class** (e.g., `AlgebraBase`) that provides the shared interface (`zero`, `one`, `+`, `*`, `chart()`, etc.) without the semiring contract. Both `Semiring` and `Interval` inherit from it.
  4. **Move `Interval` out of the package** if the package should be strictly semirings.

## Weighted formal languages (WFSA/WCFG) semiring
- Add a semiring of weighted formal languages, where values are functions from strings to semiring weights (i.e., weighted languages over a base semiring).
- Addition is pointwise sum, multiplication is weighted concatenation (convolution), zero maps everything to the base zero, one maps the empty string to the base one.
- Could be backed by WFSAs for compact/closed representation, connecting to the existing `wfsa` dependency.
- Similarly, consider a WCFG-backed variant for context-free weighted languages.
- See dissertation S3.3 "Formal Languages" / "Weighted formal languages" for the definition.

## Code-quality audit (session review)

Rough edges identified while reviewing the current state. Not urgent individually, but worth tracking.

### Organization

- **Split `misc.py`** — currently a catch-all of seven unrelated semirings. Suggested partition: `semirings/provenance.py` (Why, Lineage, Bridge, VBridge), `semirings/fuzzy.py` (Lukasiewicz, ThreeValuedLogic), `semirings/string.py` (String, longest_common_prefix), `semirings/division.py` (Division). `make_set`, `minmax`, `maxmin`, `dual` are factories — put them with other factories (see next).
- **Decide the fate of `sample.py`** — a one-line shim re-exporting `sampling.lazy2.Sample` for backwards-compat with `experimental/sample_powerseries.py`. Either update experimental to import from `semirings.sampling.lazy2` and delete `sample.py`, or document the shim's purpose in a comment.
- **`experimental/` state** — the files `powerseries.py` and `sample_powerseries.py` sit untracked in git. Either track them properly (commit) or move to a separate working area. Right now they're in neither state.

### Unity

- **Unify the factory-construction pattern.** There are currently eight factories following slightly different conventions: `make_semiring` (base.py), `make_set` / `minmax` / `maxmin` / `dual` (misc.py), `make_vector` (vector.py), `make_lazysort_semiring` (lazysort.py), `make_expectation` (expectation.py). They all do "build a class parametrically" but with different call conventions, different naming, and different ways of setting `zero`/`one`. Consider a single `make_semiring(name, **kwargs)` pattern or a `@semiring` decorator.
- **Clarify the `star` convention.** Base class offers three algorithms (`star_approx`, `star_fixpoint`, `star_doubling`) as free helpers; concrete classes pick one, or roll their own closed form (MaxPlus handles `score > 0` with a sentinel; MinPlus returns `-inf` for negative cost; ConvexHull raises `NotImplementedError`; Expectation now has a closed form; Lukasiewicz uses the inherited fixpoint). Document in `base.py` which algorithm is the default and when each override applies; consider naming a `star_closed` convention for classes that provide it.

### Elegance

- **`zero = None; one = None` in the base class** is semantically clear (subclasses must override) but type-ugly. Pyright flags every override. Options: sentinel values, explicit `@abstractproperty`, or just a class variable with a strong docstring.
- **`__hash__ = None`** on Vector, Expectation, SecondOrderExpectation — the semantic is right (structural `__eq__` over float-ish values → no sane hash) but proliferated. Consider a small mixin `class StructuralEqNoHash: __hash__ = None` for explicitness.
- **Two equality conventions cohabitate.** LogVal uses `np.allclose`; CutSets/String use structural equality; Float uses value equality; some types inherit identity equality and rely on the `metric()` tolerance fallback in `assert_equal`. Works, but the reader has to juggle three models. At minimum, document per class which convention applies.

### Cleanliness

- **No type annotations anywhere.** Pyright diagnostics are noisy because of it (`lift` override mismatches, `zero: None` assignments, etc.). Consider a type-annotation pass on `base.py` + the concrete semirings at minimum. Trade-off: generic types over `Semiring[T]` imply some design decisions about what `T` is.
- **Python-2 vestiges:** `__div__ = __truediv__` in `vector.py` no longer serves any purpose (Python 3 only). Safe to delete.
- **`make_semiring`'s class body** (base.py:158-194) is a 40-line class nested inside a function. Hard to read. Either extract the inner class logic to a module-level base + add attributes via the factory, or add a section-comment to make the structure visible.

### Tests

- **Add section separators to `tests/test_semirings.py`** — 1252 lines, navigable but getting long. A handful of `# ===== <theme> =====` comments would help (e.g., "Per-semiring direct tests", "Metric helpers", "AXIOM_CASES", "Matrix-semiring tests", "WeightedGraph tests").
- **The top-of-file import wall** is one `from semirings import (...)` with 25+ names. Group them (axiomatic: semirings; utility: check_*, make_set; graph: WeightedGraph, scc_decomposition) to make the intent scannable.
- **Module-level test fixtures** `_MinMax`, `_MaxMin`, `_Set` look like private module constants but are test fixtures. Add a comment like `# Test fixtures used by AXIOM_CASES below.`
- **Edge-case coverage gaps.** Examples worth adding: `LogVal` near float underflow, `ConvexHull` on degenerate colinear points, `RegularLanguage` on alphabet-sized regexes, `WeightedGraph` on disconnected components, `MinPlus`/`MaxPlus` behavior at ±inf inside sums.

## Add `metric()` to semiring interface
- In progress: adding `metric(self, other)` methods to semiring classes for approximate equality testing.
- Currently added to `Wrapped`, `make_semiring`, `LogVal`, `MaxPlus`, `MaxTimes`, `MinPlus`, `MinTimes`, and `Dual`.
- Test helper `assert_equal` updated to use metric-based tolerance instead of exact equality.
- Still needed: add to `Semiring` base class, `Float`, `Interval`, `Boolean`, matrix semiring, etc.
