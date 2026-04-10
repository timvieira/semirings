# TODO

## ~~1. Modernize packaging~~ DONE
- ~~Switch from `setup.py` to `pyproject.toml` (PEP 517/518/621)~~
- ~~Fix `project_url` (wrong key; should be `url` or `project_urls`)~~
- ~~Add `content_type='text/markdown'` for long description~~
- ~~Add a `LICENSE` file (MIT)~~
- ~~Move `pytest` and `graphviz` from `install_requires` to `[project.optional-dependencies] dev`~~
- ~~Delete old `setup.py`~~

## 2. Reduce hard dependencies
- Inline `arsenal` utilities (`log1pexp`/`log1mexp` in `logval.py`, `sorted_union`/`sorted_product` in `lazysort.py`) — they're small enough to vendor
- Make `regex.py` optional (lazy import) so the core package doesn't require `fsa`/`wfsa`
- Make `mert.py` optional so `scipy` isn't required at install time

## 3. Clean up `__init__.py`
- Replace `from semirings.misc import *` with explicit imports
- Add `__all__` to `__init__.py` and all submodules

## 4. Fix `make_semiring` `__eq__` bug
- `base.py:130`: `__eq__` uses `isinstance(other, self.x.__class__)` which checks against the wrapped value's type, not the semiring wrapper type. Should be `isinstance(other, SemiringWrapper)`.

## 5. Make `Dual` inherit from `Semiring`
- The `dual()` factory in `misc.py` creates a class that doesn't subclass `Semiring`, so it lacks `chart()`, `sum()`, `product()`, `__pow__`, `star_approx`, etc.

## 6. Remove commented-out code
- `base.py`: commented-out `__lt__`, `Chart.__add__`, `Chart.__sub__`, `Chart.round`
- `misc.py`: commented-out `WrappedBackpointers` class

## 7. Extract shared backpointer logic
- MaxPlus, MinPlus, MaxTimes, MinTimes all duplicate nearly identical backpointer tracking (`d` argument). Extract into a reusable base class.

## 8. Add type annotations
- Add type hints to the `Semiring` base class at minimum
- Consider making `Semiring` generic: `Semiring[T]`

## ~~9. Add CI~~ DONE
- ~~Add GitHub Actions workflow (`.github/workflows/ci.yml`) — pytest on Python 3.9–3.12~~
- Note: tests currently fail to collect due to pre-existing `ImportError: cannot import name 'EPSILON' from 'wfsa'`

## 10. Fix `multiplicity` shadowing in `minmax`/`maxmin`
- `minmax` and `maxmin` in `misc.py` define `multiplicity` as an instance method, shadowing the classmethod inherited from `Semiring`.

## 11. Other semirings sitting around in other projects
- centralize my collection of semirings that are currently scatter across many projects.
- missing first- and second-order expectation semirings

## 12. Fix broken `wfsa` import in tests
- `tests/test_semirings.py:12`: `from wfsa import WFSA, EPSILON` fails — `EPSILON` doesn't exist in the installed `wfsa` package. Tests can't even collect.

## 13. Fix `Float.__init__`
- `Float.__init__()` takes no arguments, so `Float(3)` raises `TypeError`. The README examples don't work.

## 14. Better implementations of chart
- I think I have a better implementation of a chart

## 15. Interval is not a semiring
- `Interval` inherits from `Semiring` but violates distributivity (it's sub-distributive: `a * (b + c) ⊆ a*b + a*c`, not equal). The test (`todo_interval`) is disabled for this reason.
- Options:
  1. **Keep as-is, document clearly.** Interval arithmetic gives sound over-approximations in semiring-style algorithms. Just improve the docs and re-enable the test with sub-distributivity-aware assertions.
  2. **Stop inheriting from `Semiring`.** `Interval` already overrides most methods; it only inherits `chart()`, `product()`, `__pow__`, and `star_approx`/`star_fixpoint` (and overrides `star`).
  3. **Introduce a lighter base class** (e.g., `AlgebraBase`) that provides the shared interface (`zero`, `one`, `+`, `*`, `chart()`, etc.) without the semiring contract. Both `Semiring` and `Interval` inherit from it.
  4. **Move `Interval` out of the package** if the package should be strictly semirings.

## 16. Weighted formal languages (WFSA/WCFG) semiring
- Add a semiring of weighted formal languages, where values are functions from strings to semiring weights (i.e., weighted languages over a base semiring).
- Addition is pointwise sum, multiplication is weighted concatenation (convolution), zero maps everything to the base zero, one maps the empty string to the base one.
- Could be backed by WFSAs for compact/closed representation, connecting to the existing `wfsa` dependency.
- Similarly, consider a WCFG-backed variant for context-free weighted languages.
- See dissertation §3.3 "Formal Languages" / "Weighted formal languages" for the definition.
