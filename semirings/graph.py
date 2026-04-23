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

from arsenal import Integerizer


class WeightedGraph:
    """Sparse weighted graph with SCC-decomposed closure and linear solvers."""

    def __init__(self, WeightType):
        self.WeightType = WeightType
        self.N = set()
        self.incoming = defaultdict(set)
        self.outgoing = defaultdict(set)
        self.E = WeightType.chart()

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
            self.incoming[j].add(i)
            self.outgoing[i].add(j)
        return self

    def edge(self, i, j, w=None):
        """Add `w` to the edge (i, j); lift (i, j) as the default weight."""
        self.N.add(i)
        self.N.add(j)
        self.E[i, j] += w if w is not None else self.WeightType.lift((i, j))
        if self.E[i, j] != self.WeightType.zero:
            self.incoming[j].add(i)
            self.outgoing[i].add(j)

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

        if len(N) == 1:
            i = N[0]
            return {(i, i): self.WeightType.star(A[i, i])}

        star = self.WeightType.star
        one = self.WeightType.one
        zero = self.WeightType.zero

        # Start with only the N×N restriction of A so that cross-block edges
        # (present in A but not in N) don't leak into the returned closure.
        B = self.WeightType.chart()
        for i in N:
            for k in N:
                B[i, k] = A[i, k]

        for j in N:
            sjj = star(B[j, j])
            for i in N:
                x = B[i, j] * sjj
                if x == zero:
                    continue
                for k in N:
                    B[i, k] += x * B[j, k]

        for i in N:
            B[i, i] += one

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
