# EXPLORE

A working list of [semirings](https://en.wikipedia.org/wiki/Semiring), applications,
and results to dig into together. Each entry includes enough context to recover
*why* it's on the list. Check items off as we implement, prove, or write them up.

## Semirings to add

- [ ] **[Expectation semiring](https://scholar.google.com/scholar?q=Eisner+2002+Parameter+Estimation+Probabilistic+Finite-State+Transducers)** (Eisner 2002).
      Pairs $(p, r)$ with $(p_1,r_1)\oplus(p_2,r_2)=(p_1{+}p_2,\ r_1{+}r_2)$ and
      $(p_1,r_1)\otimes(p_2,r_2)=(p_1 p_2,\ p_1 r_2 + p_2 r_1)$. One forward pass yields $\mathbb{E}[r]$ over derivations.
      *Why it's here:* the cleanest "algebra computes the gradient" story outside AD. First question: implement it, reproduce the forward-backward-as-expectation reduction.

- [ ] **[Higher-order expectation semiring](https://scholar.google.com/scholar?q=Li+Eisner+2009+First+Second+Order+Expectation+Semirings)** (Li & Eisner 2009).
      Tuples carrying $p, r, s, rs, r^2, \ldots$ for variance, covariance, Hessian-vector products.
      *Why:* variance of a feature over derivations — without sampling — in one pass. First question: how does it compose with lazy / truncated evaluation?

- [ ] **[Provenance semiring](https://scholar.google.com/scholar?q=Green+Karvounarakis+Tannen+2007+Provenance+Semirings)** $\mathbb{N}[X]$ (Green, Karvounarakis, Tannen, PODS 2007). See also [data lineage](https://en.wikipedia.org/wiki/Data_lineage).
      Symbolic polynomial over tuple-IDs; project homomorphically to Bool (why-provenance), $\mathbb{N}$ (counting), tropical (trust), security lattice.
      *Why:* "parse once, interpret many" — a direct analogue to what semiring parsing does for DPs, but for databases. First question: can we get a clean demo showing the homomorphism-projection story end-to-end?

- [ ] **[Bottleneck / widest-path](https://en.wikipedia.org/wiki/Widest_path_problem)** $(\mathbb{R}, \max, \min)$.
      *Why:* dual to shortest-path but *not* isomorphic to min-plus — distributivity holds for different reasons. Good teaching example that the semiring axioms don't pin down the "shape" you expect. First question: does our Kleene-star implementation specialize correctly here?

- [ ] **Fuzzy / [Gödel](https://en.wikipedia.org/wiki/T-norm)** $([0,1], \max, \min)$ and **Łukasiewicz** $([0,1], \max, \cdot)$.
      *Why:* idempotent $\otimes$ in the Gödel case changes which fixed-point iterations converge; concrete bridge to [t-norm fuzzy logic](https://en.wikipedia.org/wiki/T-norm_fuzzy_logics). First question: which metric axioms from `metric()` survive on $[0,1]$?

- [ ] **Counting semiring** $(\mathbb{N}, +, \times)$ as its own type, not aliased to `float`.
      *Why:* derivation counts blow up factorially; exactness is load-bearing for pedagogy (and for catching bugs where someone "counts" in float and doesn't notice). First question: integrate with `cutsets` to count min-cuts exactly.

- [ ] **Viterbi-derivation** — best score paired with its witness derivation, as a proper semiring (not the `lazysort` sorting trick).
      *Why:* lets us separate "the algebra" from "the implementation." First question: is it actually a semiring, or does the tie-breaking rule sneak in non-associativity?

- [ ] **Free / polynomial semiring over an alphabet** $\mathbb{N}\langle\langle \Sigma^*\rangle\rangle$ — [formal power series](https://en.wikipedia.org/wiki/Formal_power_series).
      *Why:* the initial object that every other weighted-automaton computation factors through. Gives a clean semantic target for your `regex` / `kleene` modules. First question: what's the right finite representation (truncated-degree? rational expressions?) to make this usable?

- [ ] **Min-max / max-min on lattices**.
      *Why:* bottleneck assignment, reliability analysis, capacity planning — all reduce here, and the lattice structure is different enough from the real-line semirings to be instructive. First question: generalize the existing ones to work over any total order.

- [ ] **Signed log-semiring** — log-space with a sign bit.
      *Why:* everyone reinvents this (poorly) when doing gradient-like sums in log space with negative contributions. Cataloging it precisely is a small public good. First question: what's the numerically-stable $\oplus$ that handles $\log 0$, sign flips, and near-cancellation together?

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
