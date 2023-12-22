# TODO: fsa should move into utilities since it is not a semiring - it is a data
# structure that supports the semiring for formal languages
from arsenal import Integerizer, colors
from collections import defaultdict
from functools import lru_cache


def dfs(Ps, arcs):
    stack = list(Ps)
    m = FSA()
    for P in Ps: m.add_start(P)
    while stack:
        P = stack.pop()
        for a, Q in arcs(P):
            if Q not in m.nodes:
                stack.append(Q)
                m.nodes.add(P)
            m.add(P, a, Q)
    return m


_frozenset = frozenset
class frozenset(_frozenset):
    def __repr__(self):
        return '{%s}' % (','.join(str(x) for x in sorted(self)))


class FSA:

    def __init__(self):
        self.start = set()
        self.edges = defaultdict(lambda: defaultdict(set))
        self.nodes = set()
        self.stop = set()
        self.syms = set()

    def as_tuple(self):
        return (frozenset(self.nodes),
                frozenset(self.start),
                frozenset(self.stop),
#                frozenset(self.syms),
                frozenset(self.arcs()))

    def __hash__(self):
        return hash(self.as_tuple())

    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()

    def __repr__(self):
#        return f'<{self.__class__.__name__} id={id(self)}>'
        return repr(self.to_regex())

    def __str__(self):
        x = ['{']   # todo: better print; include start/stop
        for s in sorted(self.nodes):
            ss = f'{s}'
            if s in self.start:
                ss = f'^{ss}'
            if s in self.stop:
                ss = f'{ss}$'
            x.append(f'  {ss}:')
            for a, t in sorted(self.arcs(s)):
                x.append(f'    {a} -> {t}')
        x.append('}')
        return '\n'.join(x)

    def _repr_html_(self):
        return self.graphviz()._repr_image_svg_xml()

    def graphviz(self, show_label=True):
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

        # FIXME: make sure this name is actually unqiue
        start = '<start>'
        g.node(start, label='', shape='point', height='0', width='0')
        for i in self.start:
            g.edge(start, str(f(i)), label='')

        for i in sorted(self.nodes):
            shape = 'circle'
            if i in self.stop: shape = 'doublecircle'
            label = str(i) if show_label else ''
            #if i in self.start: label = '*'
            g.node(str(f(i)), label=label, shape=shape)
            for a, j in sorted(self.arcs(i)):
                g.edge(str(f(i)), str(f(j)), label=str(a))

        return g

    def D(self, x):
        "left derivative"
        m = FSA()

        e = self.epsremoval()
        for i,a,j in e.arcs():

            if i in e.start and a == x:
                m.add(i,eps,j)
            else:
                m.add(i,a,j)

        m.start = set(e.start)
        m.stop = set(e.stop)
        return m

    def add(self, i, a, j):
        self.edges[i][a].add(j)
        self.nodes.add(i); self.syms.add(a); self.nodes.add(j)
        return self

    def add_start(self, i):
        self.start.add(i)
        self.nodes.add(i)
        return self

    def add_stop(self, i):
        self.stop.add(i)
        self.nodes.add(i)
        return self

    def arcs(self, i=None, a=None):
        if i is None and a is None:

            for i in self.edges:
                for a in self.edges[i]:
                    for j in self.edges[i][a]:
                        yield (i,a,j)

        elif i is not None and a is None:

            for a in self.edges[i]:
                for j in self.edges[i][a]:
                    yield (a,j)

        elif i is not None and a is not None:

            for j in self.edges[i][a]:
                yield j

        else:
            raise NotImplementedError()

    def reverse(self):
        m = FSA()
        for i in self.start:
            m.add_stop(i)
        for i in self.stop:
            m.add_start(i)
        for i, a, j in self.arcs():
            m.add(j, a, i)     # pylint: disable=W1114
        return m

    def _accessible(self, start):
        return dfs(start, self.arcs).nodes

    def accessible(self):
        return self._accessible(self.start)

    @lru_cache(None)
    def trim(self):
        c = self.accessible() & self.reverse().accessible()
        m = FSA()
        for i in self.start & c:
            m.add_start(i)
        for i in self.stop & c:
            m.add_stop(i)
        for i,a,j in self.arcs():
            if i in c and j in c:
                m.add(i,a,j)
        return m

    def renumber(self):
        return self.rename(Integerizer())

    def rename(self, f):
        "Note: f is not bijective, states may split/merge."
        m = FSA()
        for i in self.start:
            m.add_start(f(i))
        for i in self.stop:
            m.add_stop(f(i))
        for i, a, j in self.arcs():
            m.add(f(i), a, f(j))
        return m

    def rename_apart(self, other):
        f = Integerizer()
        self = self.rename(lambda i: f((0, i)))
        other = other.rename(lambda i: f((1, i)))
        assert self.nodes.isdisjoint(other.nodes)
        return (self, other)

    def __mul__(self, other):
        m = FSA()
        self, other = self.rename_apart(other)
        m.start = self.start
        m.stop = other.stop
        for i in self.stop:
            for j in other.start:
                m.add(i,eps,j)
        for i,a,j in self.arcs():
            m.add(i,a,j)
        for i,a,j in other.arcs():
            m.add(i,a,j)
        return m

    def __add__(self, other):
        m = FSA()
        [self, other] = self.rename_apart(other)
        m.start = self.start | other.start
        m.stop = self.stop | other.stop
        for i,a,j in self.arcs():
            m.add(i,a,j)
        for i,a,j in other.arcs():
            m.add(i,a,j)
        return m

    def p(self):
        "self^+"
        m = FSA()
        m.start = set(self.start)
        m.stop = set(self.stop)
        for i,a,j in self.arcs():
            m.add(i,a,j)
        for i in self.stop:
            m.add_stop(i)
            for j in self.start:
                m.add(i, eps, j)
        return m

    def star(self):
        "self^*"
        return one + self.p()

#    def L(self, s):
#        assert s in self.nodes
#        return dfs({s}, self.arcs)

    @lru_cache(None)
    def epsremoval(self):

        eps_m = FSA()
        for i,a,j in self.arcs():
            if a == eps:
                eps_m.add(i,a,j)

        @lru_cache
        def eps_accessible(i):
            return eps_m._accessible({i})

        m = FSA()

        for i,a,j in self.arcs():
            if a == eps: continue
            m.add(i, a, j)
            for k in eps_accessible(j):
                m.add(i, a, k)

        for i in self.start:
            m.add_start(i)
            for k in eps_accessible(i):
                m.add_start(k)

        for i in self.stop:
            m.add_stop(i)
            for k in eps_accessible(i):
                m.add_stop(k)

        return m

    @lru_cache(None)
    def det(self):

        self = self.epsremoval()

        def powerarcs(Q):
            for a in self.syms:
                yield a, frozenset({j for i in Q for j in self.edges[i][a]})

        m = dfs([frozenset(self.start)], powerarcs)

        for powerstate in m.nodes:
            if powerstate & self.stop:
                m.add_stop(powerstate)

        return m

    def min_brzozowski(self):
        "Brzozowski's minimization algorithm"
        # https://en.wikipedia.org/wiki/DFA_minimization#Brzozowski's_algorithm

        # Proof of correctness:
        #
        # Let M' = M.r.d.r
        # Clearly,  [[M']] = [[M]]
        #
        # In M', there are no two states that can accept the same suffix
        # language because the reverse of M' is deterministic.
        #
        # The determinization of M' then creates powerstates, where every pair
        # of distinct powerstates R and S, there exists by construction at least
        # one state q of M' where q \in R and q \notin S. Such a q contributes
        # at least one word w \in [[q]] to the suffix language of q in [[R]]
        # that is not present in [[S]], since this word is unique to q (i.e., no
        # other state accepts it).  Thus, all pairs of states in M'.d are
        # distinguishable.
        #
        # Thus, after trimming of M'.d, we have a DFA with no indistinguishable
        # or unreachable states, which is must minimal.

        return self.reverse().det().reverse().det().trim()

    def min_fast(self):
        self = self.det().renumber()

        # calculate inverse of transition function (i.e., reverse arcs)
        inv = defaultdict(set)
        for i,a,j in self.arcs():
            inv[j,a].add(i)

        final = self.stop
        nonfinal = self.nodes - final

        P = [final, nonfinal]
        W = [final, nonfinal]

        while W:
            A = W.pop()
            for a in self.syms:
                X = {i for j in A for i in inv[j,a]}
                R = []
                for Y in P:
                    if X.isdisjoint(Y) or X >= Y:
                        R.append(Y)
                    else:
                        YX = Y & X
                        Y_X = Y - X
                        R.append(YX)
                        R.append(Y_X)
                        W.append(YX if len(YX) < len(Y_X) else Y_X)
                P = R

        # create new equivalence classes of states
        minstates = {}
        for i, qs in enumerate(P):
            #minstate = frozenset(qs)
            for q in qs:
                minstates[q] = i #minstate

        return self.rename(lambda i: minstates[i]).trim()

    def min_faster(self):
        self = self.det().renumber()

        # calculate inverse of transition function (i.e., reverse arcs)
        inv = defaultdict(set)
        for i,a,j in self.arcs():
            inv[j,a].add(i)

        final = self.stop
        nonfinal = self.nodes - final

        P = [final, nonfinal]
        W = [final, nonfinal]

        find = {i: block for block, elements in enumerate(P) for i in elements}

        while W:

            A = W.pop()
            for a in self.syms:

                X = {i for j in A for i in inv[j,a]}

                blocks = {find[i] for i in X}

                for block in blocks:
                    Y = P[block]

                    if X >= Y: continue

                    # TODO: use indexing to find nonempty (Y-X).
                    # Some notes:
                    #  - To be nonempty we need to find an element i* that is in Y
                    #    but not in X.  We already have an element i that is in Y
                    #    and X.
                    #  - We know that X and Y overlap (thanks to our indexing
                    #    trick).  Now, we want to filter out cases where X >= Y
                    #    because they do not need to be split.

                    YX = Y & X
                    Y_X = Y - X

                    # we will replace block with the intersection case (no
                    # need to update `find` index for YX elements)
                    P[block] = YX

                    new_block = len(P)
                    for i in Y_X:
                        find[i] = new_block

                    P.append(Y_X)
                    W.append(YX if len(YX) < len(Y_X) else Y_X)

        return self.rename(lambda i: find[i]).trim()

    min = lru_cache(None)(min_faster)

    def equal(self, other):
        return self.min()._dfa_isomorphism(other.min())

    def _dfa_isomorphism(self, other):
        "Find isomorphism between DFAs (if one exists)."

        # Requires that self and other are minimal DFAs

        # Theorem. If `self` and `other` are graphs with out-degree at most 1, then
        # the DFA works to determine whether G and H are isomorphic

        # A deterministic machine has exactly one start state

        # Two minimized DFAs are input
        # If the number of states is differs, these machines cannot be isomorphic
        if len(self.nodes) != len(other.nodes): return False
        if len(self.start) == 0: return len(other.start) == 0

        assert len(self.start) == 1 and len(other.start) == 1

        #self = self.renumber()
        #other = other.renumber()

        [p] = self.start
        [q] = other.start

        stack = [(p, q)]
        iso = {p: q}

        syms = self.syms | other.syms

        done = set()
        while stack:
            (p, q) = stack.pop()
            done.add((p,q))
            for a in syms:

                # presences of the arc has to be the same
                if (a in self.edges[p]) != (a in other.edges[q]):
                    return False

                if a not in self.edges[p]:
                    continue

                # machines are assumed deterministic
                [r] = self.edges[p][a]
                [s] = other.edges[q][a]

                if r in iso and iso[r] != s:
                    return False

                iso[r] = s
                if (r,s) not in done:
                    stack.append((r,s))

        #print('-----')
        #print(iso)
        #print(self)
        #print(self.rename(lambda x: iso.get(x, -1)))
        #print(other)
        #print('-----')

        return self.rename(iso.get) == other

    def to_regex(self):
        import numpy as np
        from semirings.regex import Symbol
        from semirings.kleene import kleene

        n = len(self.nodes)

        A = np.full((n,n), Symbol.zero)
        start = np.full(n, Symbol.zero)
        stop = np.full(n, Symbol.zero)

        ix = Integerizer(list(self.nodes))

        for i in self.nodes:
            for a, j in self.arcs(i):
                if a == eps:
                    A[ix(i),ix(j)] += Symbol.one
                else:
                    A[ix(i),ix(j)] += Symbol(a)

        for i in self.start:
            start[ix(i)] += Symbol.one

        for i in self.stop:
            stop[ix(i)] += Symbol.one

        return start @ kleene(A, Symbol) @ stop

    def __and__(self, other):
        "intersection"

        self = self.epsremoval().renumber()
        other = other.epsremoval().renumber()

        def product_arcs(Q):
            (q1, q2) = Q
            for a, j1 in self.arcs(q1):
                for j2 in other.edges[q2][a]:
                    yield a, (j1,j2)

        m = dfs({(q1, q2) for q1 in self.start for q2 in other.start},
                product_arcs)

        # final states
        for q1 in self.stop:
            for q2 in other.stop:
                m.add_stop((q1, q2))

        return m

    def add_sink(self, syms):
        "constructs a complete FSA"

        syms = set(syms)

        self = self.renumber()

        sink = len(self.nodes)
        for a in syms:
            self.add(sink, a, sink)

        for q in self.nodes:
            if q == sink: continue
            for a in syms - set(self.edges[q]):
                if a == eps: continue  # ignore epsilon
                self.add(q, a, sink)

        return self

    def __sub__(self, other):
        return self & other.invert(self.syms | other.syms)

    __or__ = __add__

    def __xor__(self, other):
        "Symmetric difference"
        return (self | other) - (self & other)

    def invert(self, syms):
        "create the complement of the machine"

        self = self.det().add_sink(syms)

        m = FSA()

        for i in self.nodes:
            for a, j in self.arcs(i):
                m.add(i, a, j)

        for q in self.start:
            m.add_start(q)

        for q in self.nodes - self.stop:
            m.add_stop(q)

        return m

    def __floordiv__(self, other):
        "left quotient self//other ≐ {y | ∃x ∈ other: x⋅y ∈ self}"

        # TODO: support NFA/epsilon arcs?
        self = self.epsremoval()
        other = other.epsremoval()

        # quotient arcs are very similar to product arcs except that the common
        # string is "erased" in the new machine.
        def quotient_arcs(Q):
            (q1, q2) = Q
            for a, j1 in self.arcs(q1):
                for j2 in other.edges[q2][a]:
                    yield eps, (j1, j2)

        m = dfs({(q1, q2) for q1 in self.start for q2 in other.start},
                quotient_arcs)

        # If we have managed to reach a final state of q2 then we can move into
        # the post-prefix set of states
        for (q1,q2) in set(m.nodes):
            if q2 in other.stop:
                m.add((q1, q2), eps, (q1,))

        # business as usual
        for q1 in self.nodes:
            for a, j1 in self.arcs(q1):
                m.add((q1,), a, (j1,))
        for q1 in self.stop:
            m.add_stop((q1,))

        return m

    def __truediv__(self, other):
        "right quotient self/other ≐ {x | ∃y ∈ other: x⋅y ∈ self}"
        return (self.reverse() // other.reverse()).reverse()   # reduce to left quotient on reversed languages

    def __lt__(self, other):
        "self ⊂ other"
        if self.equal(other): return False
        return (self & other).equal(self)

    def __le__(self, other):
        "self ⊆ other"
        return (self & other).equal(self)

    @classmethod
    def lift(cls, x):
        m = cls()
        m.add_start(0); m.add_stop(1); m.add(0,x,1)
        return m

    @classmethod
    def from_string(cls, xs):
        m = cls()
        m.add_start(xs[:0])
        for i in range(len(xs)):
            m.add(xs[:i], xs[i], xs[:i+1])
        m.add_stop(xs)
        return m

    @classmethod
    def from_strings(cls, Xs):
        m = cls()
        for xs in Xs:
            m.add_start(xs[:0])
            for i in range(len(xs)):
                m.add(xs[:i], xs[i], xs[:i+1])
            m.add_stop(xs)
        return m

    def __contains__(self, xs):
        d = self.det()
        [s] = d.start
        for x in xs:
            t = d.edges[s][x]
            if not t: break
            [s] = t
        return (s in d.stop)

    def merge(self, S, name=None):
        "merge states in `S` into a single state."
        if name is None: name = min(S)
        def f(s):
            return name if s in S else s
        m = FSA()
        for x in self.start:
            m.add_start(f(x))
        for x,a,y in self.arcs():
            m.add(f(x),a,f(y))
        for x in self.stop:
            m.add_stop(f(x))
        return m

eps = 'ε'

FSA.one = one = FSA()
one.add_start(0); one.add_stop(0)

FSA.zero = zero = FSA()


def test_visualization():
    a,b = map(FSA.lift, 'ab')
    ((a+b).star())._repr_html_()


def test_intersection():
    a,b,c = map(FSA.lift, 'abc')
    assert a.equal((a + a) & a)
    assert b.equal((a + b) & (b + c))
    assert b.equal((a.star() + b.star()) & (b + c))
    assert b.star().equal((a+b).star() & (b.star() + c.star()))
    assert ((a*b*c) & (a+b+c).star()).equal(a*b*c)


def test_complement():
    a,b,c = map(FSA.lift, 'abc')

    have = a - a
    want = zero
    assert have.equal(want)

    have = (a+b) - b
    want = a
    assert have.equal(want)

    have = (a+b) - (b+c)
    want = a
    assert have.equal(want)

    have = (a + b.star()) - (b+c)
    want = a + b * b * b.star() + one
    assert have.equal(want)

    # 1 b bb bbb bbbb bbbbb bbbbbb
    # 1   bb     bbbb       bbbbbb
    have = b.star() - (b*b).star()
    want = b * (b * b).star()
    assert have.equal(want)


def test_equality():
    a,b,c = map(FSA.lift, 'abc')

    assert a.equal(a)
    assert not a.equal(b)

    x = ((one + a) * a.star() * a.star() * a.star() * b)
    y = x.min()
    assert x.equal(y)
    assert y.equal(x)

    assert not (a * a.star()).equal(a.star())

    assert (one + a.star()).equal(a.star())
    assert (b * one + b * a.star()).equal(b * a.star())

    assert (a * a.star()).equal(a * a.star() + a * a * a.star())

    assert (a + b) != (b + a)
    assert (a + b + c).equal(b + a + c)

    assert ((a+b).star()).equal((a.star()*b.star()).star())

    assert (a+a).equal(a)   # idempotent

    x = a*a + a
    y = (a*a + a)*(a*a + a) + a

    assert not x.equal(y)


def test_min():
    a,b = map(FSA.lift, 'ab')
    z = zero
    for x1 in [a,b]:
        for x2 in [a,b]:
            for x3 in [a,b]:
                z += (x1 * x2 * x3).min()
    assert len(z.min().nodes) == 4
    print(z.min())

    assert len(((a + b)*(a + b)*(a + b)).min().nodes) == 4

    assert len(((one + a) * a.star() * a.star() * a.star() * b).min().nodes) == 2


def test_fsa_to_regex():
    from semirings.regex import Symbol
    a, b = map(FSA.lift, 'ab')
    m = b * a.star()
    m = m.min()

    A, B = map(Symbol, 'ab')
    assert m.to_regex() == B + B * A.star() * A


def test_quotient():
    a, b, c, d = map(FSA.lift, 'abcd')

    q = (a * b) // a
    assert q.equal(b)

    q = (a * b * c + c * d) // a
    assert q.equal(b * c + zero * c * d)

    q = (a * b * c * d + c * a * b * d) // (a * b)
    assert q.equal(c * d)

    q = (a * b) / b
    assert q.equal(a)

    q = ((a * b) / b) / a
    assert q.equal(one)

    q = (a * b) / (a * b * c)
    assert q.equal(zero)

    q = (a.star() * b) // a
    assert q.equal(a.star() * b)

    q = (a.star() * b) // a.star()
    print(q)
    print(q.min().renumber())
    assert q.equal(a.star() * b)   # is this correct?


    # L1\L2 = {y | ∃x ∈ L2: xy ∈ L1}
    # L1/L2 = {x | ∃y ∈ L2: xy ∈ L1}

    L1 = (a.star() * b)
    L2 = a.star()

    def checker(y): return (L2 * y) <= L1

    assert checker(a * b)
    assert checker(a.star() * b)
    assert not checker(c)


if __name__ == '__main__':
    from arsenal import testing_framework
    testing_framework(globals())
