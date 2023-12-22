"""
Implementation of weighted finite state automata over fields.

This implementation of equivalence and minimization follows Kiefer (2020) very closely


References

 - Stefan Kiefer (2020) "Notes on Equivalence and Minimization of Weighted Automata"
   arXiv:2009.01217v1 [cs.FL]


"""
import numpy as np
from arsenal import Integerizer, colors
from collections import defaultdict, Counter, deque
from functools import cached_property

from numpy import linalg
#from scipy import linalg

#cached_property = property


EPSILON = ε = "ε"


class WFSA:
    """
    Weighted finite state automata where weights are a field (e.g., real-valued).
    """
    def __init__(self):
        self.alphabet = set()
        self.states = set()
        self.delta = defaultdict(lambda: defaultdict(Counter))
        self.start = Counter()
        self.stop = Counter()

    def __repr__(self):
        return f"WFSA({self.dim} states)"

    def __hash__(self):
        return hash(self.simple)

    def __str__(self):
        output = []
        output.append('{')
        for p in self.states:
            output.append(f'  {p} \t\t({self.start[p]}, {self.stop[p]})')
            for a, q, w in self.arcs(p):
                output.append(f'    {a}: {q}\t[{w}]')
        output.append('}')
        return "\n".join(output)

    @property
    def dim(self):
        return len(self.states)

    def add_state(self, q):
        self.states.add(q)

    def add_arc(self, i, a, j, w):
        self.add_state(i)
        self.add_state(j)
        self.alphabet.add(a)
        self.delta[i][a][j] += w

    def add_I(self, q, w):
        self.add_state(q)
        self.start[q] += w

    def add_F(self, q, w):
        self.add_state(q)
        self.stop[q] += w

    def __call__(self, xs):
        "Evaluate the sequence `xs`."
        return self.simple(xs)

    @property
    def I(self):
        for q, w in self.start.items():
            if abs(w) > 1e-5:
                yield q, w

    @property
    def F(self):
        for q, w in self.stop.items():
            if abs(w) > 1e-5:
                yield q, w

    def arcs(self, i=None):
        if i is None:
            for i in self.delta:
                for a, T in self.delta[i].items():
                    for j, w in T.items():
                        if abs(w) > 1e-5:
                            yield i, a, j, w
        else:
            for a, T in self.delta[i].items():
                for j, w in T.items():
                    if abs(w) > 1e-5:
                        yield a, j, w

    def rename(self, f):
        "Note: If `f` is not bijective, states may merge."
        m = self.spawn()
        for i, w in self.I:
            m.add_I(f(i), w)
        for i, w in self.F:
            m.add_F(f(i), w)
        for i, a, j, w in self.arcs():
            m.add_arc(f(i), a, f(j), w)
        return m

    def rename_apart(self, other):
        f = Integerizer()
        return (self.rename(lambda i: f((0,i))), other.rename(lambda i: f((1,i))))

    @cached_property
    def renumber(self):
        return self.rename(Integerizer())

    def spawn(self, *, keep_init=False, keep_arcs=False, keep_stop=False):
        "returns a new FSA in the same semiring"
        m = WFSA()
        if keep_init:
            for q, w in self.I:
                m.add_I(q, w)
        if keep_arcs:
            for i,a,j,w in self.arcs():
                m.add_arc(i, a, j, w)
        if keep_stop:
            for q, w in self.F:
                m.add_F(q, w)
        return m

    @cached_property
    def reverse(self):
        "creates a reversed machine"
        # create the new machine
        R = self.spawn()
        # reverse each arc
        for i, a, j, w in self.arcs():
            R.add_arc(j, a, i, w)
        # reverse initial and final states
        for q, w in self.I:
            R.add_F(q, w)
        for q, w in self.F:
            R.add_I(q, w)
        return R

    def __add__(self, other):

        # non-essential optimizations
        if self is zero: return other
        if other is zero: return self

        self, other = self.rename_apart(other)
        U = self.spawn(keep_init=True, keep_arcs=True, keep_stop=True)
        # add arcs, initial and final states from argument
        for i, a, j, w in other.arcs():
            U.add_arc(i, a, j, w)
        for q, w in other.I:
            U.add_I(q, w)
        for q, w in other.F:
            U.add_F(q, w)
        return U

    def __sub__(self, other):
        self, other = self.rename_apart(other)
        U = self.spawn(keep_init=True, keep_arcs=True, keep_stop=True)
        # add arcs, initial and final states from argument
        for q, w in other.I:            U.add_I(q, -w)
        for i, a, j, w in other.arcs(): U.add_arc(i, a, j, w)
        for q, w in other.F:            U.add_F(q, w)
        return U

    def __mul__(self, other):
        if not isinstance(other, WFSA): return other.__rmul__(self)

        # non-essential optimizations
        if self is one: return other
        if other is one: return self
        if other is zero: return zero
        if self is zero: return zero

        self, other = self.rename_apart(other)
        C = self.spawn(keep_init=True, keep_arcs=True)
        # add arcs, initial and final states from argument
        for i, a, j, w in other.arcs():
            C.add_arc(i, a, j, w)
        for q, w in other.F:
            C.add_F(q, w)
        # connect the final states from `self` to initial states from `other`
        for i1, w1 in self.F:
            for i2, w2 in other.I:
                C.add_arc(i1, EPSILON, i2, w1 * w2)
        return C

    def star(self):
        return self.one + self.kleene_plus()

    def kleene_plus(self):
        "self^+"
        m = self.spawn(keep_init=True, keep_arcs=True, keep_stop=True)
        for i, w1 in self.F:
            for j, w2 in self.I:
                m.add_arc(i, EPSILON, j, w1*w2)
        return m

    def _repr_svg_(self):
#        return self.graphviz()._repr_svg_()
        return self.graphviz()._repr_image_svg_xml()

    def graphviz(self, fmt=lambda x: f'{round(x,3):g}', fmt_node=lambda x: ' '):
        from graphviz import Digraph
        g = Digraph(
            graph_attr=dict(rankdir='LR'),
            node_attr=dict(
                fontname='Monospace',
                fontsize='10',
                height='.05', width='.05',
                #margin="0.055,0.042"
                margin="0,0"
            ),
            edge_attr=dict(
                #arrowhead='vee',
                arrowsize='0.3',
                fontname='Monospace',
                fontsize='9'
            ),
        )
        f = Integerizer()
        for i,w in self.I:
            start = f'<start_{i}>'
            g.node(start, label='', shape='point', height='0', width='0')
            g.edge(start, str(f(i)), label=f'{fmt(w)}')
        for i in self.states:
            g.node(str(f(i)), label=str(fmt_node(i)), shape='circle')
        for i,w in self.F:
            stop = f'<stop_{i}>'
            g.node(stop, label='', shape='point', height='0', width='0')
            g.edge(str(f(i)), stop, label=f'{fmt(w)}')
        for i, a, j, w in sorted(self.arcs()):
            g.edge(str(f(i)), str(f(j)), label=f'{a}/{fmt(w)}')
        return g

    @cached_property
    def simple(self):
        self = self.renumber

        S = self.dim
        start = np.zeros(S)
        arcs = {a: np.zeros((S,S)) for a in self.alphabet}
        stop = np.zeros(S)

        for i, w in self.I:
            start[i] += w
        for i, a, j, w in self.arcs():
            arcs[a][i,j] += w
        for i, w in self.F:
            stop[i] += w

        if EPSILON in arcs:
            W = arcs.pop(EPSILON)   # remove it
            E = linalg.inv(np.eye(self.dim) - W)
            for a in arcs:
                arcs[a] = arcs[a] @ E
            start = start @ E

        return Simple(start, arcs, stop)

    def __eq__(self, other):
        return self.simple == other.simple

    def counterexample(self, other):
        return self.simple.counterexample(other.simple)

    @cached_property
    def min(self):
        return self.simple.min.to_wfsa()

    def multiplicity(self, m):
        return WFSA.lift(EPSILON, m) * self

#    def push(self):
#        """
#        Weight pushing
#        """
#        V = self.backward()
#        p = self.spawn()
#        for i in self.states:
#            p.add_I(i, self.start[i] * V[i])
#            p.add_F(i, 1/V[i] * self.stop[i])
#            for a, j, w in self.arcs(i):
#                p.add_arc(i, a, j, 1/V[i] * w * V[j])
#        return p
#
#    def backward(self):
#        b = self.R.chart()
#        K = self.matrix_star()
#        for p in self.states:
#            for q in self.states:
#                b[p] += K[p, q] * self.stop[q]
#        return b

    @classmethod
    def lift(cls, x, w):
        m = cls()
        m.add_I(0, 1)
        m.add_arc(0, x, 1, w)
        m.add_F(1, 1)
        return m

    @classmethod
    def from_string(cls, xs):
        m = cls()
        m.add_I(xs[:0], 1)
        for i in range(len(xs)):
            m.add_arc(xs[:i], xs[i], xs[:i+1], 1)
        m.add_F(xs, 1)
        return m

    @classmethod
    def from_strings(cls, Xs):
        m = cls()
        for xs in Xs:
            m.add_I(xs[:0])
            for i in range(len(xs)):
                m.add_arc(xs[:i], xs[i], xs[:i+1], 1)
            m.add_F(xs)
        return m


WFSA.zero = zero = WFSA()
WFSA.one = one = WFSA.lift(EPSILON, 1)


class Simple:
    def __init__(self, start, arcs, stop):
        assert EPSILON not in arcs
        self.start = start
        self.arcs = arcs
        self.stop = stop
        [self.dim] = start.shape

    def __call__(self, xs):
        forward = self.start.copy()
        for x in xs:
            if x not in self.arcs: return 0
            forward = forward @ self.arcs[x]
        return forward @ self.stop

    def __eq__(self, other):
#        return (self is other) or hash(self) == hash(other) and self.counterexample(other) is None
        return self.counterexample(other) is None

    def __hash__(self):
#        return self._hash
        return 0

    # TODO: not working yet :-(
    @cached_property
    def _hash(self):
        # Assign a random value in to each symbols of the alphabet, then compute
        # the weight of all paths.
        #
        # The idea for this hash function is related to this randomized
        # algorithm for equivalence testing: https://arxiv.org/abs/1302.2818
        #
        # TODO: We might be able to skip the matrix inverse! Rather than summing
        # over all paths, we can basically just sample some paths.  That's more
        # similar to what is done in that paper.
        self = self.min
        # this gives us a weighted graph where (i,a:w,j) ==> (i,a*w,j)
        n = self.dim
        W = np.zeros((n,n))
        for a,M in self.arcs.items():
            W += hash(a) * M
        z = self.start @ linalg.solve(np.eye(n) - W, self.stop)
        return hash(int(np.round(z)))

    def counterexample(self, B):

        alphabet = set(self.arcs) | set(B.arcs)

        va = self.start @ self.stop
        vb = B.start @ B.stop
        if not approx_equal(va, vb):
            return ((), va, vb)

        eta = np.hstack([self.stop, B.stop])
        if approx_equal(eta, 0):
            return  # both empty; thus, equivalent.

        worklist = deque()
        worklist.append(((), self.stop, B.stop))
        basis = [eta]

        while worklist:

            (w, VA, VB) = worklist.pop()

            for a in alphabet:

                ua = self.arcs[a] @ VA if a in self.arcs else 0*VA
                ub = B.arcs[a] @ VB if a in B.arcs else 0*VB

                w1 = (a, w)

                # counterexample?
                va = self.start @ ua
                vb = B.start @ ub
                if not approx_equal(va, vb):
                    return (w1, va, vb)   # yup!

                u = np.hstack([ua, ub])
                q = proj(u, basis)

                # add to the basis if it's not redundant; updates to the basis
                # go on the worklist.
                if not approx_equal(u - q, u):
                    worklist.append((w1, ua, ub))
                    basis.append(q)

    @cached_property
    def min(self):
        # Proposition 3.4: Forward (backward) conjugates are forward (backward) minimal.
        # Proposition 3.5: The backward (forward) conjugate of a forward (backward) minimal automaton is minimal.
        return self.forward_conjugate().backward_conjugate()

    def __repr__(self):
        return f'<Simple states={len(self.start)}, syms={len(self.arcs)}>'

    def forward_conjugate(self):
        # The forward basis F is an n' x n where n is the set of old states, and n' the set of new states
        #
        # we want to solve
        #
        #   old_start = new_start    F
        #       n          n'     (n' x n)
        #
        # There are n >= n' constraints, but they aren't linearly independent so
        # the system has a unique solution.
        #
        #    old_start @ F.T = new_start @ F @ F.T
        #                                 (n' x n) (n x n')
        #        n
        #
        #    old_start @ F.T @ inv(F @ F.T) = new_start
        F = self.forward_basis()

        #P = F.T @ linalg.inv(F @ F.T)
        P = linalg.pinv(F)

        #start = self.start @ P
        #assert np.allclose(self.start, start @ F)

        # We also need to solve for our new transition function
        #                   F M = M' F
        #               F M F.T = M' (F F.T)
        #    F M F.T inv(F F.T) = M' (F F.T) inv(F F.T)
        return Simple(
            start = self.start @ P,
            arcs = {a: F @ M @ P for a, M in self.arcs.items()},   # apply change of basis
            stop = F @ self.stop,
        )

    @cached_property
    def reverse(self):
        return Simple(
            start = self.stop,
            arcs = {a: M.T for a, M in self.arcs.items()},
            stop = self.start,
        )

    def forward_basis(self):
        worklist = [self.start]
        basis = [self.start]
        while worklist:
            V = worklist.pop()
            for a in self.arcs:
                u = V @ self.arcs[a]
                q = proj(u, basis)
                if not approx_equal(u - q, u):
                    worklist.append(u)
                    basis.append(q)
        return np.array(basis)

    def backward_conjugate(self):
        return self.reverse.forward_conjugate()

    def to_wfsa(self):
        m = WFSA()
        for i in range(self.dim):
            m.add_I(i, self.start[i])
        for a in self.arcs:
            for i in range(self.dim):
                for j in range(self.dim):
                    m.add_arc(i,a,j,self.arcs[a][i,j])
        for i in range(self.dim):
            m.add_F(i, self.stop[i])
        return m


def approx_equal(x,y):
    return np.allclose(x,y)


def proj(u, Q):
    # simple implementation of the Gram–Schmidt process
    # https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    for q in Q:
        u = u - ((q @ u) / (q @ q)) * q
    return u


def test_min():

    a = WFSA.lift('a', 1)
    b = WFSA.lift('b', 1)
    c = WFSA.lift('c', 1)

    print(colors.yellow % 'Example 1')
    M = one + a + a * b + a * b * c
    m = M.min
    assert m.dim == 4, m.dim

    print(colors.yellow % 'Example 2')
    z = zero
    for x1 in [a,b]:
        for x2 in [a,b]:
            for x3 in [a,b]:
                z += (x1 * x2 * x3)
    assert z.min.dim == 4

    print(colors.yellow % 'Example 3')
    assert ((a + b)*(a + b)*(a + b)).min.dim == 4

    print(colors.yellow % 'Example 4')
    M = ((one + a * a.star()) * b)
    m = M.min
    #print(m)
    #print(m.to_wfsa().push())

    assert m.dim == 2, m.dim


def test_equivalence():

    a = WFSA.lift('a', 1)
    b = WFSA.lift('b', 1)
    c = WFSA.lift('c', 1)


    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 0')
    A = a.star()
    B = a.star()

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    assert np.allclose(A('a'), 1)
    assert np.allclose(A('a'), 1)

    assert np.allclose(A('b'), 0)
    assert np.allclose(B('b'), 0)

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 1')
    A = a * (b + c)
    B = a * b + a * c

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    assert np.allclose(A('ab'), 1), A('ab')
    assert np.allclose(A('ac'), 1)
    assert np.allclose(A('ab'), 1)
    assert np.allclose(B('ac'), 1)

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 2')
    A = a * (b + c + c)
    B = a * b + a * c

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is not None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 3')
    A = a * (b + c + c)
    B = a * b + a * c + a * c

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 4')
    A = a * b.star() * c.star()
    B = a * (one + b * b.star()) * c.star()

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 5')
    A = a.star()
    B = one + a * a.star()

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 6')

    A = ((one + a * a.star()) * b)
    B = (a.star() * b)
    assert A.counterexample(B) is None

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 7')
    A = (a * b).star()
    B = one + a * (b * a).star() * b

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None

    #___________________________________________________________________________
    #
    print(colors.yellow % 'Example 8')
    a = WFSA.lift('a', .1)
    b = WFSA.lift('b', .2)
    c = WFSA.lift('c', .3)

    A = (a + b).star()
    B = (a.star() * b).star() * a.star()

    w = A.counterexample(B)
    #print('counterexample?', w)
    assert w is None



if __name__ == '__main__':
    from arsenal import testing_framework
    testing_framework(globals())
