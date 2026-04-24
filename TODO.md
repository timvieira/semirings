# TODO

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
- Missing first- and second-order expectation semirings

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

## Add `metric()` to semiring interface
- In progress: adding `metric(self, other)` methods to semiring classes for approximate equality testing.
- Currently added to `Wrapped`, `make_semiring`, `LogVal`, `MaxPlus`, `MaxTimes`, `MinPlus`, `MinTimes`, and `Dual`.
- Test helper `assert_equal` updated to use metric-based tolerance instead of exact equality.
- Still needed: add to `Semiring` base class, `Float`, `Interval`, `Boolean`, matrix semiring, etc.
