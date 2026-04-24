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
