# Unifying semirings: design & staged plan

## Context

Semiring implementations are currently spread across three repos:

- `~/projects/semirings/` (canonical — this repo)
- `~/projects/hypergraphs/hypergraph/semirings/` (~1,800 LOC, derivation-tracking style)
- `~/projects/algebra/algebra/number/` (richer algebraic hierarchy — Ring/Field/ClosedSemiring/Module, plus matrix algorithms built on top)

Goal: centralize concrete semiring types in `semirings/`, eliminate duplicate implementations, absorb the genuinely-unique contributions from each repo (sampling semirings, Expr, Bag, Count, dual numbers, expectation semirings), and do it in stages so easy wins land first and harder design work happens once with a clear framework.

The harder design work is derivation tracking. Most of this document is about that framework — the mechanical merges are secondary.

**Migration discipline:** each port is immediately paired with a downstream redirect, so no type lives in two places for longer than the item that moves it. See Part III for the per-item protocol.

---

# Part I — The derivation-tracking framework

## The key idea: pairing as a first-class operation

Hypergraphs tracks derivations by hanging a `.d` field on every semiring instance. That's ad-hoc — it bakes one specific choice of derivation-tracker into every concrete type.

The right abstraction is **pairing**: given a semiring $S$ and a compatible companion structure $D$, construct a new semiring whose values are pairs $(s, d)$. Derivation tracking, automatic differentiation, expectation semirings, and several provenance patterns are all instances of pairing — they just differ in what $D$ is and how the pair's multiplication is defined.

Two kinds of pairing cover essentially all the cases in the three repos:

### 1. `Product(S, D)` — componentwise pairing

- Elements: $(s, d)$ with $s \in S$, $d \in D$
- Addition: $(s_1, d_1) + (s_2, d_2) = (s_1 + s_2,\; d_1 + d_2)$
- Multiplication: $(s_1, d_1) \cdot (s_2, d_2) = (s_1 s_2,\; d_1 d_2)$
- $D$ must itself be a semiring
- This is the categorical product in **Semiring**

Use: attach a parallel "provenance" or "derivation" semiring whose operations are independent of $S$'s values.

### 2. `Jet(S, M, order=k)` — truncated-polynomial / twisted pairing

- Elements: truncated polynomials $s_0 + s_1 \varepsilon + \ldots + s_k \varepsilon^k$, represented as $(k+1)$-tuples of coefficients in appropriate types
- Ring $S[\varepsilon]/(\varepsilon^{k+1})$ with $\varepsilon$ formally nilpotent of order $k+1$
- Addition is coefficientwise; multiplication is polynomial convolution truncated at degree $k$
- $M$ (the coefficient type for the higher-order terms) must be an $S$-module — i.e., admit left and right $S$-actions compatible with multiplication
- First-order ($k=1$): $(s_1, d_1) \cdot (s_2, d_2) = (s_1 s_2,\; s_1 d_2 + d_1 s_2)$

Use: differentiation, sensitivity, expectations — any setting where the "extra" component depends linearly on the base value.

### Why these two suffice

Every pairing-flavored construction in the three repos falls into one of these buckets:

| Existing thing | Is really |
|---|---|
| `semirings.Dual` | `Jet(S, S, order=1)` with a specific interpretation |
| `algebra.number.Dual` | `Jet(Real, Real, order=1)` for AD |
| `algebra.number.Expectation` | `Jet(Real, Real, order=1)` for E-step |
| `hypergraph.semirings.SecondOrderExpectation` | `Jet(Real, Real, order=2)` |
| `hypergraph.semirings.Expectation` | `Jet(Real, Real, order=1)` |
| `hypergraph.semirings.*` with `.d` field | `Product(S, FreeExpr)` (or `Product(S, FreeClosedSemiring)`) |

Three ad-hoc classes collapse to `Jet` aliases. The per-class `.d` field collapses to `Product(S, D)` for whichever derivation tracker $D$ you actually want.

Pairings compose: `Product(Product(S, D₁), D₂)` and `Jet(Product(S, D), M, k)` are all well-defined.

---

## What goes in the `D` slot: the free-semiring zoo

For derivation tracking specifically, the interesting question is what to use for $D$. The answer: "a representation of the free (closed) semiring over a set of generators $X$," with the realization depending on what properties you want and what equality you can afford.

Four realizations, three tiers of rigor:

### Tier 0 — magma: `FreeExpr(X)`

- Binary trees with `+`, `*`, `0`, `1` at internal nodes and generators at leaves
- Equality is syntactic (tree equality)
- Not a semiring — associativity, distributivity, and identity laws fail as tree equalities
- Port: hypergraphs' `Expr` class (≈230 LOC), renamed and documented as magma
- Use: tracing, debugging, pretty-printing derivations. Semantic equality is the user's problem.

### Tier 1 — free semiring: WFSAs over $\mathbb{N}$

The free semiring over $X$ is $\mathbb{N}\langle X\rangle$ — non-commutative polynomials over $\mathbb{N}$ — and the free *closed* semiring is $\mathbb{N}$-rational series $\mathbb{N}\langle\langle X\rangle\rangle_{\mathrm{rat}}$. Both have a canonical computational realization as WFSAs with weights in $\mathbb{N}$:

- Finite polynomials = acyclic WFSAs
- Rational series (star supported) = general WFSAs

Equality is decidable in polynomial time for $\mathbb{N}$-weighted WFSAs (Schützenberger). Minimization gives a canonical form.

**Implementation: a thin wrapper on the existing `wfsa` dependency.** `FreeClosedSemiring(X)` exposes a `WFSA[X, ℕ]` with semiring operations on top of it. No new algorithms — the hard parts (minimization, equivalence) live in `wfsa`.

### Tier 1' — free idempotent closed semiring: `RegularLanguage(X)`

This is the existing `semirings/regex.py`. Same as Tier 1 but with idempotence ($a + a = a$), so weights collapse from $\mathbb{N}$ to Boolean, and WFSAs reduce to NFAs. Equality via FSA minimization. Already implemented.

### Tier 2 — commutative provenance: `Why`, `Lineage`, `Bag`

For cases where the order of generators doesn't matter (a provenance annotation is a set or multiset of sources, not a sequence). These are the free *commutative* closed semiring over $X$ at two granularities — already implemented in `semirings/misc.py`.

### The picture

```
                        commutative ?                   can afford equivalence ?
                          /       \                         /         \
                     no                yes              yes              no
                      |                 |                |                |
           (non-comm free)      (commutative free)  (decidable)       (magma)
                      |                 |                |                |
            WFSAs over ℕ           Why / Lineage     RegularLanguage   FreeExpr
           (or Boolean for       (boolean / multiset)   (NFA / FSA)      (trees)
            idempotent: NFAs)
```

All four fit into a `Product(S, D)` with $S$ = whatever base semiring you care about (MaxPlus, LogVal, Real, ...).

---

## Selective semirings and backpointers

For selective semirings like MaxPlus, "+" picks a winner, which means the derivation tracker can project away the losing branch at every Add node. The result is a single product chain — the classical backpointer list.

This is a homomorphism
$$\mathrm{FreeExpr}(X) \twoheadrightarrow \mathrm{SelectedDerivation}(X)$$
parameterized by a choice function on sums. So `Product(MaxPlus, SelectedDerivation)` recovers the classical "MaxPlus with backpointers" semiring, and `Product(MaxPlus, FreeExpr)` gives you the full parse/proof forest. Same construction, different $D$.

`SelectedDerivation(X)` is trivially a monoid (concatenation of generators) promoted to a semiring via the sum-selection homomorphism; ~30 LOC.

---

# Part II — API sketch

## Constructions

```python
# semirings/pairing.py

def Product(S, D):
    """Categorical product of semirings. Componentwise +, *, 0, 1."""
    ...

def Jet(S, M=None, order=1):
    """Truncated-polynomial pairing S[ε]/(ε^(order+1)).
    M defaults to S. Must be an S-module when provided."""
    ...
```

## Free-semiring realizations

```python
# semirings/free.py

class FreeExpr:
    """Magma. Syntactic equality. Use for tracing only."""
    ...

class FreeClosedSemiring:
    """Free closed semiring over X. Backed by WFSA with ℕ weights.
    Decidable equality via WFSA minimization."""
    ...

class SelectedDerivation:
    """Chain of generators. Image of FreeExpr under +-selection.
    Use as D in Product(S, SelectedDerivation) for backpointers on
    selective semirings."""
    ...
```

## Retrofits

```python
# Before:
class Expectation(...):
    def __mul__(self, other):
        return Expectation(self.p * other.p,
                           self.p * other.r + self.r * other.p)

# After:
Expectation = Jet(Real, Real, order=1)
SecondOrderExpectation = Jet(Real, Real, order=2)
Dual = Jet(Real, Real, order=1)  # same math, different convention
```

## Usage pattern

```python
from semirings import MaxPlus, Product, FreeClosedSemiring, SelectedDerivation

# MaxPlus with full parse forest
TracedMaxPlus = Product(MaxPlus, FreeClosedSemiring("abc"))

# MaxPlus with Viterbi-style single backpointer
ViterbiMaxPlus = Product(MaxPlus, SelectedDerivation("abc"))

# AD over MaxPlus (uncommon but well-defined)
DiffMaxPlus = Jet(MaxPlus, MaxPlus, order=1)
```

---

# Part III — Staged plan

## Per-item protocol

Every port, merge, and retrofit in Stages 1–4 uses the same five-step protocol. Each item is independently revertible.

1. **Port / merge** into `semirings/`.
2. **Verify in-tree:** `ruff check`, `pytest`, `check_axioms` on the new / modified type.
3. **Redirect downstream:** replace the corresponding implementation in `hypergraph/semirings/` and/or `algebra/number/` with a re-export shim (`from semirings import X`).
4. **Verify downstream:** run the downstream repo's full test suite; run `pip install -e .` to confirm packaging is intact.
5. **Commit** in the downstream repo only if steps 2–4 pass. On failure, revert step 3 and investigate — the `semirings/` port still stands and subsequent items can proceed.

Cross-repo prerequisite (done once, in Stage 0): `hypergraphs` and `algebra` must declare `semirings` as a dependency in their packaging.

## Stage 0 — Prep

- Unification branch in each of the three repos
- Per-file ruff ignores in `semirings/` (done)
- Add `semirings` as a dependency in:
  - `hypergraphs/setup.py` (or `pyproject.toml`)
  - `algebra/setup.py`
- Verify dev installs: `pip install -e .` in each of the three repos succeeds
- Record baseline: run each repo's test suite and note current pass/fail state, so later regressions are attributable

Effort: ~1 hour. Risk: trivial.

## Stage 1 — Port + redirect: types unique to hypergraphs

Each row uses the per-item protocol. Order is low-risk to higher-risk.

| Item | `semirings/` destination | `hypergraph/semirings/` redirect |
|---|---|---|
| `Count` | `semirings/count.py` | replace `count.py` with shim |
| `Bag` | `semirings/bag.py` | replace `bag.py` with shim |
| Derivation-tree pretty-printer | `semirings/util.py` | replace `util.py` with shim |
| `Vector` (module factory) | `semirings/vector.py` | replace `vector.py` with shim |
| `Expr` → `FreeExpr` | `semirings/free.py` | replace `expr.py` with shim |
| Sampling trio (`Expon`, lazy `Sample`, SWOR `Sample`) | `semirings/sampling/` | replace `sampling/` subpackage with shims |
| Expectation / SecondOrderExpectation | `semirings/expectation.py` (interim; becomes `Jet(...)` aliases in Stage 3) | replace `expectation.py` with shim |

Effort: ~1 day total. Risk: low per item. Bag, Count, util, Vector are trivial; Expr, sampling, Expectation carry slightly more test surface.

## Stage 2 — Reconcile + redirect: near-duplicates with hypergraphs

Per pair: merge winner into `semirings/`, then replace hypergraphs' copy with a shim. For the merges flagged as "True merge needed" or where semantics shift, run hypergraphs' full test suite before committing step 5 — this is where most risk sits.

| Pair | Outcome |
|---|---|
| Boolean | `semirings/` wins; port `samples()` classmethod from hypergraphs |
| MaxPlus | same |
| MinPlus | same |
| MaxTimes | same |
| MinTimes | same |
| Float | same |
| LogVal | True merge: keep `semirings/`'s `metric()`; pull `__neg__()`, `to_real` alias, arsenal imports from hypergraphs |
| Bottleneck | `semirings/` wins wholesale |
| Lukasiewicz | `semirings/` wins wholesale (hypergraphs' version has a line-19 bug) |
| LazySort | `semirings/` wins wholesale (factory pattern is superior) |
| ConvexHull/Point (MERT) | adopt hypergraphs' version (identity guards, explicit Point, no Semiring inheritance) |
| Regex | keep current FSA-backed `RegularLanguage`; hypergraphs' symbolic classes were already absorbed by `FreeExpr` / `FreeClosedSemiring` in Stage 1/3 |

After Stage 2, `hypergraph/semirings/` should be entirely shims.

Effort: ~1 day. Risk: low per pair; `check_axioms` catches regressions.

## Stage 3 — Pairing framework + retrofit in place

1. Implement in `semirings/`:
   - `Product(S, D)` (~40 LOC)
   - `Jet(S, M, order=k)` (~80 LOC)
   - `SelectedDerivation` (~30 LOC)
   - `FreeClosedSemiring` wrapper on `wfsa` (~50 LOC; see open question on $\mathbb{N}$ weights)
2. Retrofit in-tree: replace `Dual`, `Expectation`, `SecondOrderExpectation` implementations with `Jet(...)` aliases; keep their existing public names.
3. Redirect downstream: `hypergraph/semirings/expectation.py` shim now points to the `Jet`-backed version. Same for `algebra/number/dual.py` and `algebra/number/expectation.py`.
4. Verify `check_axioms` and existing tests across all three repos.

Effort: 1–2 days. Risk: medium. Most risk is in `Jet` algebra correctness and the `wfsa`-over-ℕ integration.

**Open question for this stage:** does `wfsa` support arbitrary weight semirings including $\mathbb{N}$? If only floats / max-plus / log, a small upstream PR to `wfsa` is a prerequisite.

## Stage 4 — Cherry-pick + redirect from `algebra/number/`

Types unique to `algebra/number/` that fit `semirings/`. Per-item protocol as above.

- `dual.py` — already redirected in Stage 3 (subsumed by `Jet`); nothing to do here beyond deleting the stale file
- `poly.py` / `Monomial` → port to `semirings/poly.py` if there's a real use case; otherwise skip
- `BagRing` / `baggr` → port to `semirings/bagring.py` if database-provenance use case matters; otherwise skip
- **Stays in `algebra/`:** `Ring`, `Field`, `ClosedSemiring`, and every matrix algorithm that depends on them (LU, Gauss-Jordan, Bareiss, eigenvalue, Kleene, WFSA, graph). These remain local to `algebra/` because the matrix algorithms are parameterized over them, and that's `algebra/`'s core purpose.

After Stage 4, `algebra/number/`'s concrete semiring types (MaxPlus, MinPlus, Boolean, LogVal, Regex, Interval, etc.) are all shims, and the package's identity becomes clearly "linear algebra over semirings" rather than "yet another semiring zoo."

Effort: few hours to a half day. Risk: low.

---

# Part IV — Verification

Verification is distributed across stages — each item in Stages 1–4 carries its own per-item protocol (steps 2 and 4). That replaces the old "big-bang downstream smoke test at the end" approach.

Global checks, to be run at the end of each stage:

- `ruff check` is clean in all three repos
- `check_axioms` passes for every semiring and every construction (`Product`, `Jet`, `FreeClosedSemiring`, `Product(S, FreeExpr)` with magma caveat noted)
- `check_metric_axioms` passes for types that define `metric()`
- `pytest` passes in all three repos
- `pip install -e .` succeeds in all three repos
- `from semirings import *` yields a superset of the previous public names
- No orphaned files in `hypergraph/semirings/` or `algebra/number/` — every former implementation is either a shim or deleted

---

# Part V — Open questions

1. **`wfsa` over $\mathbb{N}$.** Does the existing `wfsa` package support $\mathbb{N}$ as a weight semiring? If not, small upstream PR needed before Stage 3.

2. **`FreeExpr` as pluggable $D$.** Should we allow `Product(S, FreeExpr)` even though the result isn't strictly a semiring (due to `FreeExpr` being a magma)? Proposal: yes, with a documented caveat. The `S`-projection is still a homomorphism, which is what most users care about.

3. **Commutative generators.** Do any concrete use cases want commutative generators that aren't already served by `Why`/`Lineage`/`Bag`? If not, we can punt on a dedicated "free commutative semiring" type.

4. **Scope of `Jet` module coefficients.** First-order `Jet(S, S)` is universal for AD-over-$S$. Do we need `Jet(S, M)` with $M \ne S$ (e.g., gradients lying in an $S$-module rather than $S$ itself)? Probably yes eventually, but we can start with `Jet(S, S)` and generalize later.

5. **Polymorphism under Kleene star.** Does `Product(S, D).star()` work if only `S` or only `D` is closed? If not, `ClosedProduct(S, D)` requires both to be closed semirings. Same for `Jet`.

---

# Part VI — Critical files & existing utilities to reuse

**Modify in `semirings/`:**
- `semirings/__init__.py` — expand exports
- `semirings/base.py` — possibly add `samples()` classmethod default (used by many ported types)
- `semirings/logval.py` — merge from hypergraphs
- `semirings/mert.py` — replace with hypergraphs version
- `semirings/misc.py` — may shrink as `Dual` etc. migrate to `Jet`

**Create in `semirings/`:**
- `semirings/pairing.py` — `Product`, `Jet`
- `semirings/free.py` — `FreeExpr`, `FreeClosedSemiring`, `SelectedDerivation`
- `semirings/count.py`, `semirings/bag.py`, `semirings/vector.py`, `semirings/util.py`
- `semirings/sampling/` — subpackage
- `semirings/expectation.py` — interim in Stage 1, then `Jet(...)` alias in Stage 3

**Replace with shims (downstream):**
- `hypergraph/semirings/*.py` — nearly every file becomes a one-line `from semirings import X`
- `algebra/number/{dual,expectation,logval,regex,interval,fuzzy}.py` and the relevant classes in `algebra/number/semirings.py` — shims after their Stage

**Packaging:**
- `hypergraphs/setup.py` or `pyproject.toml` — add `semirings` dep
- `algebra/setup.py` — add `semirings` dep

**Reuse (existing, don't reinvent):**
- `semirings.base.check_axioms` — axiom verification (`semirings/base.py:196`)
- `semirings.base.check_metric_axioms` — metric verification (`semirings/base.py:246`)
- `semirings.misc.Why`, `Lineage` — commutative provenance (`semirings/misc.py`)
- `semirings.regex.RegularLanguage` — free idempotent closed semiring (already the right tier)
- External `wfsa` package — backbone for `FreeClosedSemiring`
- `arsenal.maths.log1pexp`, `log1mexp` — used by `LogVal`

---

# Part VII — What I recommend starting with

Stage 0 first, always — the cross-repo prerequisite of declaring `semirings` as a dep and getting a green baseline test run in all three repos makes every subsequent per-item protocol safe to execute.

Within Stage 1, start with the smallest trivial items (`Count`, `Bag`, `util`) to warm up the per-item protocol before tackling the larger ports (`Expr`, sampling, `Expectation`).

Stage 3 (pairing framework) is the real design stage. It benefits from having `FreeExpr` in-tree and a clean base class with no drift (Stages 1–2). Doing it third avoids redesigning twice.
