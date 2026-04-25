# EXPLORE

A working list of [semirings](https://en.wikipedia.org/wiki/Semiring), applications,
and results to dig into together. Each entry includes enough context to recover
*why* it's on the list. Check items off as we implement, prove, or write them up.

## Semirings to add

- [x] **[Expectation semiring](https://scholar.google.com/scholar?q=Eisner+2002+Parameter+Estimation+Probabilistic+Finite-State+Transducers)** (Eisner 2002) — `expectation.py` (`Expectation`).
      Open follow-up: notebook reproducing the forward-backward-as-expectation reduction.

- [x] **[Higher-order expectation semiring](https://scholar.google.com/scholar?q=Li+Eisner+2009+First+Second+Order+Expectation+Semirings)** (Li & Eisner 2009) — `expectation.py` (`SecondOrderExpectation`).
      Open follow-up: how it composes with lazy / truncated evaluation.

- [ ] **[Provenance semiring](https://scholar.google.com/scholar?q=Green+Karvounarakis+Tannen+2007+Provenance+Semirings)** $\mathbb{N}[X]$ (Green, Karvounarakis, Tannen, PODS 2007). See also [data lineage](https://en.wikipedia.org/wiki/Data_lineage).
      Partial: `Why` / `Lineage` (`misc.py`) cover the lattice-projection targets, but the full $\mathbb{N}[X]$ polynomial carrier and the homomorphism-projection demo are still missing. First question: a clean end-to-end demo showing one symbolic computation projected to Bool (why-provenance), $\mathbb{N}$ (counting), tropical (trust), and a security lattice.

- [x] **[Bottleneck / widest-path](https://en.wikipedia.org/wiki/Widest_path_problem)** $(\mathbb{R}, \max, \min)$ — `misc.py` (`Bottleneck = maxmin(-inf, +inf)`).
      Open follow-up: confirm Kleene-star specializes correctly here.

- [ ] **[Gödel](https://en.wikipedia.org/wiki/T-norm)** $([0,1], \max, \min)$.
      Łukasiewicz $([0,1], \max, \cdot)$ is done (`misc.py`, `Lukasiewicz`); Gödel is the remaining t-norm. *Why:* idempotent $\otimes$ changes which fixed-point iterations converge; concrete bridge to [t-norm fuzzy logic](https://en.wikipedia.org/wiki/T-norm_fuzzy_logics). First question: which metric axioms from `metric()` survive on $[0,1]$?

- [x] **Counting semiring** $(\mathbb{N}, +, \times)$ as its own type — `count.py` (`Count`, `make_count_k`).
      Open follow-up: integrate with `cutsets` to count min-cuts exactly.

- [ ] **Viterbi-derivation** — best score paired with its witness derivation, as a proper semiring (not the `lazysort` sorting trick or the per-class `.d` backpointer field).
      Currently only the ad-hoc `.d` pattern exists. The planned replacement is `Product(MaxPlus, SelectedDerivation)` from DESIGN Stage 3 (see TODO.md "Pairing framework"). First question: is it actually a semiring, or does the tie-breaking rule sneak in non-associativity?

- [ ] **Free / polynomial semiring over an alphabet** $\mathbb{N}\langle\langle \Sigma^*\rangle\rangle$ — [formal power series](https://en.wikipedia.org/wiki/Formal_power_series).
      Partial: `FreeExpr` (`free.py`) is the syntactic free term algebra; no truncated-degree power series or rational-expression carrier yet. First question: what's the right finite representation (truncated-degree? rational expressions?) to make this usable as the "semantic target" for `regex` / `kleene`?

- [x] **Min-max / max-min on lattices** — `misc.py` (`minmax`, `maxmin` factories) work over any total order.

- [x] **Signed log-semiring** — `logval.py` (`LogVal`).
      Carries a sign bit (`pos`) and uses `log1pexp` / `log1mexp` for the four sign-cases of $\oplus$, with absorbing-zero in $\otimes$. Open follow-up: write up the numerically-stable $\oplus$ recipe — it's the "small public good" entry — and pin it down with axiom-test edge cases (near-cancellation, $\log 0$, sign flips with $\pm\infty$).

## Applications worth a notebook or writeup

- [ ] **[Interprocedural dataflow (IFDS/IDE)](https://scholar.google.com/scholar?q=Reps+Horwitz+Sagiv+1995+precise+interprocedural+dataflow)** (Reps–Horwitz–Sagiv, POPL 1995).
      *Hook:* distributive dataflow reduces to graph reachability over a semiring. Literally "semiring parsing, but for programs." First question: can we stage a small IFDS demo using the existing `regex`/`kleene` infrastructure?

- [ ] **[Tropical geometry](https://en.wikipedia.org/wiki/Tropical_geometry) / Newton polytopes** (Pachter & Sturmfels, [*Algebraic Statistics for Computational Biology*](https://scholar.google.com/scholar?q=Pachter+Sturmfels+Algebraic+Statistics+Computational+Biology)).
      *Hook:* Minkowski sum + convex hull forms a semiring — the right algebra for *parametric* shortest-path ("how does the solution change as weights vary?"). First question: implement Newton-polytope propagation and replicate a parametric-alignment figure.

- [ ] **[Network calculus](https://en.wikipedia.org/wiki/Network_calculus)** (Le Boudec & Thiran, book [free online](https://scholar.google.com/scholar?q=Le+Boudec+Thiran+Network+Calculus)).
      *Hook:* min-plus convolution gives worst-case delay bounds for queueing networks. A semiring deployed in actual industrial telecom tooling. First question: min-plus convolution vs. tropical matrix product — are they actually the same operator in disguise?

- [ ] **Coalitional game theory via polynomial semirings**.
      *Hook:* Shapley values and related indices can be phrased as evaluations of a polynomial in a commutative semiring. Surprisingly clean, rarely presented this way. First question: find the cleanest reference and see whether it unifies Shapley / Banzhaf / Shapley–Shubik.

- [ ] **[Kleene Algebra with Tests (KAT)](https://scholar.google.com/scholar?q=Kozen+Kleene+algebra+with+tests)** (Kozen).
      *Hook:* semiring framework in which imperative-program equivalence is decidable. The "RE equivalence" story scaled up to real programs. First question: work a small program-equivalence example end-to-end.

- [ ] **Pedagogical unification**: shortest path = Viterbi = CYK = DAG reachability = DAG-SAT.
      *Hook:* same algorithm, different semiring. Aimed at the audience who's used dynamic programming for years and never seen it as *one* theorem. First question: write the blog post / tutorial notebook; this is where your dissertation's framing would shine.

## Theory threads to pull on

- [ ] **[Krob's theorem](https://scholar.google.com/scholar?q=Krob+1994+equality+problem+tropical+semiring+undecidable)** (1994): the equational theory of $(\mathbb{N}\cup\{\infty\}, \min, +)$ is undecidable.
      *Hook:* no complete equivalence procedure for tropical/Viterbi DP expressions — every "optimizer" is necessarily heuristic. Encodes Post correspondence into tropical matrix products. First question: write this up with a small worked instance.

- [ ] **[Simon's factorization theorem](https://scholar.google.com/scholar?q=Imre+Simon+factorization+forest+theorem)**.
      *Hook:* carves out a decidable fragment that sits underneath Krob's undecidability — the reason real tropical tools work at all. First question: find the cleanest modern exposition.

- [ ] **[Hashiguchi's limitedness decidability](https://scholar.google.com/scholar?q=Hashiguchi+limitedness+distance+automata)**.
      *Hook:* what distance-automata minimization actually relies on in practice. First question: connect to any existing tooling we could compare against.

- [ ] **[Schützenberger's theorem](https://scholar.google.com/scholar?q=Schutzenberger+rational+series+weighted+automata)**: rational power series $\iff$ recognizable by a weighted automaton.
      *Hook:* the fundamental equivalence for weighted-automata theory. Ties the "free semiring" entry above to concrete machines.

- [ ] **[Tropical Perron–Frobenius / Karp's mean-cycle theorem](https://scholar.google.com/scholar?q=Karp+1978+minimum+mean+cycle)**.
      *Hook:* tropical eigenvalue = max cycle mean. Gives the "steady-state rate" for a min-plus linear system — the tropical analogue of PF's dominant eigenvalue. First question: implement and visualize on a small weighted graph.

- [ ] **[Eilenberg's variety theorem](https://scholar.google.com/scholar?q=Eilenberg+variety+theorem+recognizable+series)** for recognizable series.
      *Hook:* the algebraic classification that makes "what can a weighted automaton recognize?" a clean question. First question: the weighted generalization specifically — which references are canonical?

- [ ] **$\omega$-continuous semirings and [quantales](https://en.wikipedia.org/wiki/Quantale)**.
      *Hook:* the right setting for Kleene's star / fixed-point theorems to hold in full generality. Cleans up what `kleene.py` currently assumes. First question: audit the `kleene` module against the standard $\omega$-continuous axioms.
