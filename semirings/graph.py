"""
Weighted graphs and algebraic path problems over closed semirings.

The core primitive is `WeightedGraph`, a sparse weighted adjacency structure
equipped with SCC-decomposed algorithms for:

  * `closure()` — reflexive, transitive Kleene closure (matrix star).
  * `solve_left(b)`  — solve x = x A + b.
  * `solve_right(b)` — solve x = A x + b.

Both solvers run in block-upper-triangular order over the SCC decomposition, so
tree/DAG-structured graphs take linear time and graphs with a few large SCCs
pay O(|block|³) only within each block.

Ported from genparse's linear.py (same author).
"""

import html
from collections import defaultdict
from functools import cached_property

import numpy as np
from arsenal import Integerizer

from semirings.base import Semiring
from semirings.kleene import kleene


class WeightedGraph(Semiring):
    """Sparse weighted graph that is also a semiring under edgewise addition,
    relation composition, and Kleene star (`closure`)."""

    def __init__(self, WeightType):
        self.WeightType = WeightType
        self.N = set()
        self.outgoing = defaultdict(set)
        self.E = WeightType.chart()

    def __add__(self, other):
        """Edgewise semiring sum over the union of node sets."""
        if not isinstance(other, WeightedGraph):
            return NotImplemented
        assert self.WeightType is other.WeightType, \
            f'WeightType mismatch: {self.WeightType} vs {other.WeightType}'
        g = WeightedGraph(self.WeightType)
        g.N |= self.N | other.N
        for (i, j), w in self.E.items():
            g[i, j] = w
        for (i, j), w in other.E.items():
            g[i, j] = g.E[i, j] + w
        return g

    def __mul__(self, other):
        """Relation composition: `(self * other)[i, j] = Σ_k self[i, k] · other[k, j]`."""
        if not isinstance(other, WeightedGraph):
            return NotImplemented
        assert self.WeightType is other.WeightType, \
            f'WeightType mismatch: {self.WeightType} vs {other.WeightType}'
        S = self.WeightType
        g = WeightedGraph(S)
        g.N |= self.N | other.N
        for (i, k), w1 in self.E.items():
            for j in other.outgoing.get(k, ()):
                g[i, j] = g.E[i, j] + w1 * other.E[k, j]
        return g

    def star(self):
        """Kleene star. Same as `closure()` — reflexive, transitive closure."""
        return self.closure()

    @property
    def zero(self):
        """Empty graph — additive identity under `+` (no edges, no nodes)."""
        return WeightedGraph(self.WeightType)

    @property
    def one(self):
        """Identity graph over `self.N` — multiplicative identity restricted to
        this graph's node universe. Composition preserves `self.N`."""
        g = WeightedGraph(self.WeightType)
        for i in self.N:
            g[i, i] = self.WeightType.one
        return g

    def metric(self, other):
        """Max entrywise distance over the union of edge supports."""
        if not isinstance(other, WeightedGraph):
            return float('inf')
        if self.WeightType is not other.WeightType:
            return float('inf')
        S = self.WeightType
        keys = set(self.E) | set(other.E)
        if not keys:
            return 0.0
        return max(S.metric(self.E.get((i, j), S.zero),
                            other.E.get((i, j), S.zero))
                   for (i, j) in keys)

    def __eq__(self, other):
        if not isinstance(other, WeightedGraph):
            return NotImplemented
        return (self.WeightType is other.WeightType
                and self.N == other.N
                and dict(self.E) == dict(other.E))

    def __hash__(self):
        return hash((self.WeightType,
                     frozenset(self.N),
                     frozenset(self.E.items())))

    def __iter__(self):
        return iter(self.E)

    def __getitem__(self, item):
        i, j = item
        return self.E[i, j]

    def __setitem__(self, item, value):
        i, j = item
        self.N.add(i)
        self.N.add(j)
        if value != self.WeightType.zero:
            self.E[i, j] = value
            self.outgoing[i].add(j)
        return self

    def edge(self, i, j, w=None):
        """Add `w` to the edge (i, j); lift (i, j) as the default weight."""
        self.N.add(i)
        self.N.add(j)
        self.E[i, j] += w if w is not None else self.WeightType.lift((i, j))
        if self.E[i, j] != self.WeightType.zero:
            self.outgoing[i].add(j)

    @cached_property
    def incoming(self):
        """Reverse adjacency, derived from `E` on first access.

        Read-only: mutating the graph after reading `incoming` makes the cache
        stale (same contract as `blocks`, `Blocks`, `buckets`).
        """
        inc = defaultdict(set)
        for (i, j) in self.E:
            inc[j].add(i)
        return inc

    def transpose(self):
        """Return a new WeightedGraph with every edge reversed."""
        T = WeightedGraph(self.WeightType)
        for (i, j), w in self.E.items():
            T[j, i] = w
        return T

    def closure(self):
        """Reflexive, transitive closure as a new WeightedGraph."""
        C = WeightedGraph(self.WeightType)
        K = self.closure_scc_based()
        for i, j in K:
            C[i, j] += K[i, j]
        return C

    def closure_reference(self):
        """Naive Floyd–Warshall closure over all nodes; useful as a correctness oracle."""
        return self._closure(self.E, self.N)

    def closure_scc_based(self):
        K = self.WeightType.chart()
        for i in self.N:
            b = self.WeightType.chart()
            b[i] = self.WeightType.one
            sol = self.solve_left(b)
            for j in sol:
                K[i, j] = sol[j]
        return K

    def solve_left(self, b):
        """Solve `x = x A + b` via block-upper-triangular back-substitution."""
        sol = self.WeightType.chart()
        for block, B in self.Blocks:
            enter = self.WeightType.chart()
            for j in block:
                enter[j] += b[j]
                for i in self.incoming[j]:
                    enter[j] += sol[i] * self.E[i, j]
            for j, k in B:
                sol[k] += enter[j] * B[j, k]
        return sol

    def solve_right(self, b):
        """Solve `x = A x + b` via block-upper-triangular back-substitution."""
        sol = self.WeightType.chart()
        for block, B in reversed(self.Blocks):
            enter = self.WeightType.chart()
            for j in block:
                enter[j] += b[j]
                for k in self.outgoing[j]:
                    enter[j] += self.E[j, k] * sol[k]
            for i, j in B:
                sol[i] += B[i, j] * enter[j]
        return sol

    def _closure(self, A, N):
        """Reflexive, transitive closure of `A` restricted to the node set `N`."""
        if not isinstance(N, list):
            N = list(N)

        S = self.WeightType

        # Short-circuit: one node — closure is just star of its self-loop.
        if len(N) == 1:
            i = N[0]
            return {(i, i): S.star(A[i, i])}

        # Build the N×N restriction as a dense ndarray, delegate to the shared
        # Kleene–Lehmann implementation (`MatrixSemiring.star` uses the same
        # routine), then unpack back into a sparse chart keyed by node names.
        n = len(N)
        M = np.empty((n, n), dtype=object)
        for ii, i in enumerate(N):
            for kk, k in enumerate(N):
                M[ii, kk] = A[i, k]

        M = kleene(M, S, reflexive=True)

        B = S.chart()
        for ii, i in enumerate(N):
            for kk, k in enumerate(N):
                v = M[ii, kk]
                if v != S.zero:
                    B[i, k] = v
        return B

    @cached_property
    def blocks(self):
        """Strongly connected components in reverse topological order."""
        return list(scc_decomposition(self.incoming.__getitem__, self.N))

    @cached_property
    def buckets(self):
        """Map from node to block index."""
        return {x: i for i, block in enumerate(self.blocks) for x in block}

    @cached_property
    def Blocks(self):
        """List of (block, within-block closure) pairs."""
        return [(block, self._closure(self.E, block)) for block in self.blocks]

    def _repr_svg_(self):
        return self.graphviz()._repr_image_svg_xml()

    def graphviz(self, label_format=str, escape=lambda x: html.escape(str(x))):
        from graphviz import Digraph
        name = Integerizer()
        g = Digraph(
            node_attr=dict(
                fontname='Monospace', fontsize='9',
                height='0', width='0', margin='0.055,0.042',
                penwidth='0.15', shape='box', style='rounded',
            ),
            edge_attr=dict(
                penwidth='0.5', arrowhead='vee', arrowsize='0.5',
                fontname='Monospace', fontsize='8',
            ),
        )
        for i, j in self:
            g.edge(str(name(i)), str(name(j)), label=label_format(self.E[i, j]))
        for i in self.N:
            g.node(str(name(i)), label=escape(i))
        return g


def scc_decomposition(successors, roots):
    """
    Tarjan (1972) strongly connected components — O(V + E), O(V) space.
    Yields components in reverse topological order (sinks first).
    """
    lowest = {}
    stack = []
    trail = set()
    t = 0

    def dfs(v):
        nonlocal t
        t += 1
        num = t
        lowest[v] = t
        trail.add(v)
        stack.append(v)

        for w in successors(v):
            if lowest.get(w) is None:
                yield from dfs(w)
                lowest[v] = min(lowest[v], lowest[w])
            elif w in trail:
                lowest[v] = min(lowest[v], lowest[w])

        if lowest[v] == num:
            C = []
            while True:
                w = stack.pop()
                trail.remove(w)
                C.append(w)
                if w == v:
                    break
            yield C

    for v in roots:
        if lowest.get(v) is None:
            yield from dfs(v)
