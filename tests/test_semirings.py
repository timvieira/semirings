import numpy as np
from arsenal import colors, assert_throws
from arsenal.iterextras import take
from semirings import (
    MinPlus, MinTimes, MaxPlus, MaxTimes, Float, CutSets, Boolean,
    Bottleneck, minmax, maxmin, LogVal, ConvexHull, Point,
    Lukasiewicz, Interval, LazySort, Dual, Bridge, Division,
    make_set, String, ThreeValuedLogic, VBridge, Wrapped,
    Why, Lineage, make_semiring
)
from semirings.regex import RegularLanguage
from semirings.wfsa import WFSA, EPSILON
from semirings.fsa import FSA


class WeightedGraph:

    def __init__(self, WeightType):
        self.N = set()
        self.WeightType = WeightType
        self.E = WeightType.chart()

    def __iter__(self):
        return iter(self.E)

    def incoming(self, j):  # TODO: use slice notation
        return {I: self.E[I,J] for I,J in self if J == j}

    def outgoing(self, i):  # TODO: use slice notation
        return {J: self.E[I,J] for I,J in self if I == i}

    def edge(self, i, j, w=None):
        self.N.add(i); self.N.add(j)
        self.E[i,j] += w if w is not None else self.WeightType.lift((i,j))

    # TODO: return another WeightedGraph that corresponds to the closure
    def __getitem__(self, item):
        i,j = item
        return self.E[i,j]

    def __setitem__(self, item, value):
        i,j = item
        self.N.add(i); self.N.add(j)
        self.E[i,j] = value
        return self

    def closure(self):
        # initialization
        A = self.E
        old = A.copy()
        for j in self.N:
            new = self.WeightType.chart()
            sjj = self.WeightType.star(old[j,j])
            for i in self.N:
                for k in self.N:
                    new[i,k] = old[i,k] + old[i,j] * sjj * old[j,k]
            old, new = new, old   # swap to repurpose space
        # post processing fix-up: add the identity matrix
        for i in self.N: old[i,i] += self.WeightType.one
        c = WeightedGraph(self.WeightType)
        c.E.update(old)
        for i,j in c.E: c.N.add(i); c.N.add(j)
        return c

    def _repr_svg_(self):
        return self.graphviz()._repr_image_svg_xml()

    def graphviz(self):
        from graphviz import Digraph

        name = Integerizer()

        g = Digraph(
            node_attr=dict(
                fontname='Monospace',
                fontsize='9',
                height='0', width='0',
                margin="0.055,0.042",
                penwidth='0.15',
                shape='box',
                style='rounded',
            ),
            edge_attr=dict(
                penwidth='0.5',
                arrowhead='vee',
                arrowsize='0.5',
                fontname='Monospace',
                fontsize='8'
            ),
        )

        def escape(x):
            if isinstance(x, frozenset): x = set(x) or {}
            if isinstance(x, float): x = f'{x:.2g}'
            x = str(x)
            # remove any ANSI color codes
            x = re.sub(r'\033\[.*?m(.*?)\033\[0m', lambda m: '`%s`' % m.group(1), x)
            return x

        for i,j in self.E:
            if self.E[i,j] == self.WeightType.zero: continue
            g.edge(str(name(i)), str(name(j)), label=escape(self.E[i,j]))

        for i in self.N:
            g.node(str(name(i)), label=escape(i))

        return g


def test_pow():
    assert Float.lift(0, None) ** 0 == 1
    assert Float.lift(1, None) ** 10 == 1
    assert Float.lift(2, None) ** 10 == 2**10
    assert Float.lift(0.5, None) ** 12 == 0.5**12

    assert (MaxPlus.lift(0.5, None) ** 12).score == 0.5 * 12
    assert (MaxPlus.lift(2, None) ** 0) == MaxPlus.one

    assert MaxTimes.lift(2, None) ** 5 == MaxTimes.lift(2 ** 5, None)


#def test_pareto():
#
#    def pareto(Y):
#        return np.array([y1
#                for y1 in Y
#                if not any(np.all(y2 < y1) for y2 in Y if y1 is not y2)
#               ])
#
#    points = np.random.uniform(size=(30,2))
#    P = pareto(points)
#
#    def show(points, c='r'):
#        if isinstance(points, Pareto): points = points.x
#        points = np.array(points)
#        P = pareto(points)
#        pl.scatter(points[:,0], points[:,1], alpha=0.5, c='b')
#        pl.scatter(P[:,0], P[:,1], c=c, alpha=0.5)
#
#    Pareto = make_semiring(
#        'Pareto',
#        lambda X,Y: pareto(list(X) + list(Y)),
#        lambda X,Y: pareto([x + y for x in X for y in Y]),
#        [],
#        [np.zeros(2)],
#    )

def test_uncertainty_sets():
    U = make_semiring(
        'uncertainty',
        plus = lambda X,Y: {x + y for x in X for y in Y},
        times = lambda X,Y: {x * y for x in X for y in Y},
        zero = {0},
        one = {1},
    )

    x = U({1})
    y = U({1})
    z = U({0,1})

    print(x, '+', y, '==', x + y)
    print(y, '*', z, '==', y * z)
    print(x, '*', z, '==', x * z)

    assert ((x + y) * z).x <= (x * z + y * z).x, [(x + y) * z, x * z + y * z]


def todo_wfsa():

    S = WFSA
    S.multiplicity = lambda x, m: x * WFSA.lift(EPSILON, m)

    a = S.lift('a', 1)
    b = S.lift('b', 1)
    c = S.lift('c', 1)

    members = [
        a,
        b,
        c,
        a.star() * b,
#        b.star(),
#        (a*b).star(),
        a * b.star() * c,
        (a + b) * c,
        (a * c + b * c),
#        (a.star() * b).star() * a.star(),
    ]

    # TODO: experimental support for hashing weighted regular languages
    if 0:
        have = set(); want = set()
        for A in members:
            for B in members:
                if hash(A) != hash(B):
                    have.add((A,B))
                if A.counterexample(B) is not None:
                    want.add((A,B))

        print('hash function:')
        precision = len(have & want) / len(have) if len(have) > 0 else 1
        recall = len(have & want) / len(want) if len(want) > 0 else 1
        print('  precision:', precision)
        print('  recall:   ', recall)

        if precision != 1:
            print()
            print(colors.light.red % 'false positives:')
            for a,b in have - want:
                #print(a.simple.to_wfsa())   # remove epsilons
                #print(b.simple.to_wfsa())

                print(a.min)
                print(b.min)

                print('hash:', hash(a), hash(b))
                print('counterexample:', a.counterexample(b))
                print()

        assert precision == 1

    check_axioms_samples(S,members)


def test_funky():
    # The semiring requires can be relaxed in many ways that still admit general
    # algorithms.  For example, we can relax the requirement that * is a monoid.
    # We can allow it to be a multi-arity operator with no structure other than
    # closure \forall x_1, \ldots, x_K \mathcal{X}: f(X_1...X_K) \in
    # \mathcal{X}.  But, we still want a kind of distributive property:
    #
    #   \forall k: \sum_{x_k} f(x_1, ..., x_k, ... x_K) = f(x_1, ..., (\sum_{x_k} x_k), ... x_K)
    #
    # Example from Gilea (2020) efficient outside computation
    # https://aclanthology.org/2020.cl-4.2/

    def F(x1, x2): return

    S = make_semiring(
        'Funky',
        lambda x,y: min(x, y),
        lambda x,y: x + np.exp(y),
        np.inf,
        1,   # ???
#        star = lambda x: 1 if x == 0 else K
#        pp = lambda x: repr(dict(x)),
#        hash = lambda x: hash(frozenset(x.items())),
#        multiplicity = multiplicity
    )

    members = list(map(S, np.linspace(-5,5,11)))
    print(members)

    check_axioms_samples(S,members,hash_=False,assoc=False,star=False)

    #print([A,B,C], A * (B + C), '~~', A * B + A * C)


def test_endomorphism():
    from collections import Counter

    S = make_semiring(
        'Endo',
        lambda f,g: lambda x: f(x) + g(x),
        lambda f,g: lambda x: f(g(x)),
        lambda _: 0,
        lambda x: x,
#        star = lambda x: 1 if x == 0 else K
#        pp = lambda x: repr(dict(x)),
#        hash = lambda x: hash(frozenset(x.items())),
#        multiplicity = multiplicity
    )

    # the tricky thing about the members is that the need to be endomorphisms
    # i.e., they have to preserve the monoid structure f(x) + f(y) = f(x + y).
    #
    # If the monoid is real addition, then f are affine functions.

    members = list(map(S, [
        lambda x: 2*x+1,
        lambda x: x*2,
        lambda x: x*3,
        lambda x: 4*x,
        lambda x: 1*x,
        lambda _: 1,
        lambda _: 3,
    ]))

    a,b,c,d,e,f,g = members

    def eval(f):
        return [f.x(x) for x in range(-5, 5)]

    assert eval((a+b)*c) == eval(a*c+b*c)

    print(eval(c*(a+b)))
    print(eval(c*a+c*b))

    assert eval(c*(a+b)) == eval(c*a+c*b)

#    for A in members:
#        for B in members:
#            for C in members:
#                check_axioms(S,A,B,C,star=False,hash_=False)


# XXX: not the "semiring needed" for the mst to work doesn't check out as a semiring!
def test_max_capacity():

    S = make_semiring(
        'MaxCapacity',
        max,
        min,
        -np.inf,
        np.inf,
        star = lambda x: np.inf
    )

    members = [
        S.zero,
        S.one,
        S(0.5),
        S(0.1),
        S(3),
        S(3.1415),
        S(2),
        S(0),
    ]

    for x in members:
        print(x)
        approx = S.star_approx(x, 100)
        analytical = S.star(x)
        print('  approx:    ', approx)
        print('  analytical:', analytical)
        assert approx == analytical, [approx, analytical]

    check_axioms_samples(S,members)


# XXX: work-in-progress
def test_semilinear():
    from collections import Counter

    def bag_union(x,y):
        z = Counter()
        for k,v in x.items(): z[k] += v
        for k,v in y.items(): z[k] += v
        return z

    def minkowski(x,y):
        z = Counter()
        for k1,v1 in x.items():
            for k2,v2 in y.items():
                assert len(k1) == len(k2)
                k3 = tuple(a + b for a,b in zip(k1,k2))
                z[k3] += v1*v2   # or should it be sum?
        return z


    from arsenal import Integerizer
    symbols = 'abcd'
    alphabet = Integerizer(list(symbols))
    emptyvector = (0,)*len(symbols)

    S = make_semiring(
        'Semilinear',
        bag_union,
        minkowski,
        Counter(),     # empty bag
        Counter({emptyvector: 0}),
#        star = lambda x: 1 if x == 0 else K
    )

    def onehot(symbol):
        x = [0]*len(symbols)
        x[alphabet(symbol)] = 1
        return S({tuple(x): 1})


    a,b,c,d = list(map(onehot, symbols))

    print(a,b,c,d)
    print(a + a*a + a*a*a + a*a*a + a*b*b + b*a*b + b*b*a)
    print(a + a*a + a*a*a + a*a*a + a*b*b + b*a*b + b*b*a)

    x = a
    print(x)
    approx = S.star_approx(x, 10)
    #analytical = S.star(x)
    print('  approx:    ', approx)
    #print('  analytical:', analytical)
    #if approx.score > 1000: approx = S.inf
    #assert approx == analytical, [approx, analytical]


def test_countK():

    # k-collapsed integers
    K = 5
    S = make_semiring(
        f'Count[{K}]',
        lambda a,b: min(K, a+b),
        lambda a,b: min(K, a*b),
        0,
        1,
        star = lambda x: 1 if x == 0 else K
    )

    members = list(map(S, range(K+1)))
    print(members)

    for x in members:
        print(x, S.star(x), S.star_fixpoint(x))

    check_axioms_samples(S,members)


def test_max_times():

    S = MaxTimes

    members = [
        S.zero,
        S.one,
        S(0.5),
        S(0.1),
        S(3),
        S(3.1415),
        S(2),
        S(0),
    ]

    for x in members:
        #print(x)
        approx = S.star_approx(x, 100)
        analytical = S.star(x)
        #print('  approx:    ', approx)
        #print('  analytical:', analytical)
        if approx.score > 1000: approx = S.inf
        assert approx == analytical, [approx, analytical]

    check_axioms_samples(S,members)


def test_ints():

    S = make_semiring(
        'Int',
        lambda a,b: a+b,
        lambda a,b: a*b,
        0,
        1,
        star = lambda x: 1 if x == 0 else np.inf,
        multiplicity = lambda x,m: x*m,
    )

    members = [
        S.zero,
        S.one,
        S(0),
        S(1),
        S(2),
        S(100),
#        S(-5),      # star fails (not well-defined)
#        S(np.inf),
    ]

    for x in members:
        print(x)
        approx = S.star_approx(x, 1000)
        analytical = S.star(x)
        print('  approx:    ', approx, [S.star_approx(x, k) for k in range(10)])
        print('  analytical:', analytical)
        if approx.x >= 1000: approx = S(np.inf)
        assert approx == analytical, [approx, analytical]

    check_axioms_samples(S,members)


def test_max_plus():

    S = MaxPlus

    members = [
        S.zero,
        S.one,
        S(0.5),
        S(0.1),
        S(3),
        S(3.1415),
        S(2),
        S(0),
        S(-3),
        S(+3),
        S(2),
        S(-1),
    ]

    for x in members:
        #print(x)
        approx = S.star_approx(x, 20_000)
        analytical = S.star(x)
        #print('  approx:    ', approx.score)
        #print('  analytical:', analytical.score)
        if approx.score > 1000: approx = S.inf
        assert approx == analytical, [approx.score, analytical.score]


    check_axioms_samples(S,members)


def test_security():

    # k-collapsed integers
    K = {-1,2,3,4}
    bot = min(K)
    top = max(K)
    S = make_semiring(
        'Security',
        min,
        max,
        top,
        bot,
        star = lambda _: bot,
    )

    members = list(map(S, K))

    for x in members:
        assert S.star(x) == S.star_fixpoint(x)

    check_axioms_samples(S,members)


def test_why_semiring():

    members = [
        Why.lift('a'),
        Why.lift('b'),
        Why.lift('c'),
    ]
    check_axioms_samples(Why,members)

    graph = WeightedGraph(Why)
    graph.edge(0,1)
    graph.edge(2,3)
    graph.edge(1,2)
    graph.edge(0,1.5)
    graph.edge(1.5,2)
    graph.edge(3,4)
    graph.edge(3,5)
    graph.edge(4,5)

    if 0:
        import pylab as pl
        graph.draw()
        pl.show()

    why = graph.closure()[0,5]   # we have a bridge at edge 2->3
    assert why == Why(frozenset({
        frozenset({(0, 1), (2, 3), (4, 5), (1, 2), (3, 4)}),
        frozenset({(0, 1.5), (2, 3), (4, 5), (3, 4), (1.5, 2)}),
        frozenset({(0, 1), (1, 2), (2, 3), (3, 5)}),
        frozenset({(0, 1.5), (2, 3), (1.5, 2), (3, 5)})
    }))


def test_lineage_semiring():

    members = [
        Lineage.lift('a'),
        Lineage.lift('b'),
        Lineage.lift('c'),
    ]

    check_axioms_samples(Lineage,members)

    graph = WeightedGraph(Lineage)
    graph.edge(0,1)
    graph.edge(2,3)
    graph.edge(1,2)
    graph.edge(0,1.5)
    graph.edge(1.5,2)
    graph.edge(3,4)
    graph.edge(3,5)
    graph.edge(4,5)

    if 0:
        import pylab as pl
        graph.draw()
        pl.show()

    c = graph.closure()

    # we use all the edges for 0->5 paths
    lineage = c[0,5]
    assert lineage == Lineage(frozenset(graph.E))

    # we only use four edges for 0->2 paths
    lineage = c[0,2]
    assert lineage == Lineage(frozenset({(0, 1), (1, 2), (1.5, 2), (0, 1.5)}))


def test_cut_sets():
    S = CutSets

    # minimality
    assert CutSets( {(2,3)}, {(1,3), (2,3)} ) == CutSets( {(2,3)} )

    a = CutSets( {(2,3)}, {(1,3), (2,4)} )
    b = CutSets( {(1,3), (2,3)}, {(1,3), (2,4)} )
    c = CutSets( {(1,3)} )
    assert a * b == CutSets({(2,3)}, {(1,3), (2,4)})
    assert a + b == CutSets({(1,3), (2,3)}, {(1,3), (2,4)})

    one = CutSets.one; zero = CutSets.zero
    assert a * one == a
    assert a * zero == zero
    assert a * (b + c) == a * b + a * c
    assert a * a == a
    assert a + a == a

    members = [
        CutSets({(2,3)}, {(1,3), (2,4)}),
        CutSets({(1,3), (2,3)}, {(1,3), (2,4)}),
        CutSets({(1,3)}),
        CutSets.zero,
        CutSets.one,
    ]
    check_axioms_samples(S,members)


def test_edge_bridge():
    S = Bridge
    members = [
        Bridge({(2,3), (1,3), (2,4)}),
        Bridge({(1,3), (2,3), (1,3), (2,4)}),
        Bridge({(1,3), (2,3)}),
        Bridge({(1,3), (2,4)}),
        Bridge.lift((1,3)),   # creates a singleton set
        Bridge.zero,
        Bridge.one,
    ]

    check_axioms_samples(S,members)

    graph = WeightedGraph(Bridge)
    graph.edge(0,1)
    graph.edge(2,3)
    graph.edge(1,2)
    graph.edge(0,1.5)
    graph.edge(1.5,2)
    graph.edge(3,4)
    graph.edge(3,5)
    graph.edge(4,5)

    if 0:
        import pylab as pl
        graph.draw()
        pl.show()

    bridges = graph.closure()[0,5]   # we have a bridge at edge 2->3
    assert bridges == Bridge({(2,3)})


def test_vertex_bridge():
    graph = WeightedGraph(VBridge)
    graph.edge(0,1)
    graph.edge(2,3)
    graph.edge(1,2)
    graph.edge(0,1.5)
    graph.edge(1.5,2)
    graph.edge(3,4)
    graph.edge(3,5)
    graph.edge(4,5)

    # clearly, removing 0 or 5 disconnects 0 from 5; but so will 2 or 3 (see picture above)
    b = graph.closure()[0,5]

    if 0:
        import pylab as pl
        graph.draw()
        pl.show()

    assert b == VBridge({0,2,3,5})


def test_division():
    S = Division
    members = [Division(x) for x in range(10)]
    check_axioms_samples(S,members)


def test_three_valued_logic():
    S = ThreeValuedLogic

    [f,u,t] = members = S.zero, S.unk, S.one

    assert f < u < t
    assert_equal(f + f, f)
    assert_equal(f + t, t)
    assert_equal(t + t, t)
    assert_equal(f + u, u)
    assert_equal(t + u, t)
    assert_equal(u + u, u)

    assert_equal(f * f, f)
    assert_equal(f * t, f)
    assert_equal(t * t, t)
    assert_equal(f * u, f)
    assert_equal(t * u, u)
    assert_equal(u * u, u)

    check_axioms_samples(S,members)


def test_sets():
    from arsenal.maths.combinatorics import powerset

    universe = frozenset('abc')
    S = make_set(universe, frozenset())

    members = [S(frozenset(x)) for x in powerset(universe)]

    check_axioms_samples(S,members)


def test_string():
    S = String
    members = [
        S.zero,
        S.one,
        S("a"),
        S("aa"),
        S("aaa"),
        S("b"),
        S("aaab"),
        S("abab"),
    ]

    assert S("aaa") + S("aab") == S("aa")
    assert S("aaa") + S("") == S("")

    check_axioms_samples(S,members,right_distrib=False)


def test_regular_sets():

    # We needed a version of RegularLanguage with | and & defined -- those ops
    # are currently only available on FSAs and FSAs have the wrong equality
    # operation.
    class RegularLanguageSet:
        def __init__(self, fsa):
            assert isinstance(fsa, FSA)
            self.fsa = fsa
        def __or__(self, other): return RegularLanguageSet(self.fsa | other.fsa)
        def __and__(self, other): return RegularLanguageSet(self.fsa & other.fsa)
        @classmethod
        def lift(cls, x): return cls(FSA.lift(x))
        def __eq__(self, other): return self.fsa.equal(other.fsa)
        def __hash__(self): return 0
        def __repr__(self): return f'RegularLanguageSet({self.fsa.min()})'
    RegularLanguageSet.zero = RegularLanguageSet(FSA.zero)
    RegularLanguageSet.one = RegularLanguageSet(FSA.one)

    a,b,c = map(FSA.lift, 'abc')

    S = make_set(RegularLanguageSet((a + b + c).star()), RegularLanguageSet.zero)

    members = [
        S.zero,
        S.one,
        S(RegularLanguageSet(a.star())),
        S(RegularLanguageSet((a+c).star() * b + b * (a+b).star())),
        S(RegularLanguageSet(b.star())),
        S(RegularLanguageSet((a + b).star())),
    ]

    check_axioms_samples(S,members)


def test_lazysort():

    def check(a, b, T=None):
        A = list(a.take(T))
        B = list(take(T, b))
        assert A == B, f'\n\n{A}\n{B}\n'

    a = LazySort(2,'a')
    b = LazySort(3,'b')
    c = LazySort(5,'c')

    z = LazySort.zero

    check(z.star(), LazySort.one)

    check(a * LazySort.one, a)
    check(LazySort.one * a, a)
    check(LazySort.one * LazySort.one, LazySort.one)

    check((a+b) * z, z)
    check(z * (a+b), z)

    check((a+b) * LazySort.one, a+b)
    check(LazySort.one * (a+b), a+b)

    check(a * LazySort.zero, [])
    check(LazySort.zero * a, [])

    check(a + LazySort.zero, a)
    check(LazySort.zero + a, a)

    check(a * b + a * c, a * (b + c))

    check(LazySort.multiplicity(a, 0), [])
    check(LazySort.multiplicity(a, 4), [a,a,a,a])
    check(LazySort.multiplicity(a, np.inf), [a,a,a,a,a], T=5)

    check(a.star(), [LazySort.one, a, a*a, a*a*a, a*a*a*a], T=5)


def test_bottleneck():

    one = Bottleneck.one
    zero = Bottleneck.zero

    a = Bottleneck(1)
    b = Bottleneck(2)
    c = Bottleneck(3)
    d = Bottleneck(4)

    assert a.star_approx(10) == a.star()
    assert a.star_fixpoint() == a.star()

    assert Bottleneck.sum([]) == zero
    assert Bottleneck.sum([a,b,c,d]) == a + b + c + d
    assert Bottleneck.product([]) == one
    assert Bottleneck.product([a,b,c,d]) == a * b * c * d

    assert (a * b + b * c + d * c) == c
    assert (a + b) == b
    assert (a * b) == a
    assert a * one == a
    assert a * zero == zero
    assert a + zero == a
    assert b + one == one
    print('ok')


def test_minmax():
    S = minmax(None, "")
    members = list(map(S, "abcd")) + [S.zero, S.one]
    check_axioms_samples(S,members)


def test_maxmin():
    S = maxmin(None, "")
    members = list(map(S, "abcd"))
    check_axioms_samples(S,members)


def test_convex_hull():
    check_axioms(
        ConvexHull,
        ConvexHull([Point(1,1.5), Point(2,1.75), Point(3,4)]),
        ConvexHull([Point(2,2)]),
        ConvexHull([Point(1,10)]),
        star=False,
    )


def test_star_operations():
    assert Float.star(1) == np.inf
    assert MinTimes(-1).star() == MinTimes(-np.inf)
    assert MinPlus(-1).star() == MinPlus(-np.inf)
    assert MaxPlus(1).star() == MaxPlus(np.inf)


# TODO:Intervals are not an exact semiring because they are sub-distributive
# rather than distributive.  Add specific tests for that!
def todo_interval():

    check_axioms(
        Interval,
        Interval(0.5,0.5),
        Interval(-1.1,+1.1),
        Interval(0.25,0.75),
    )

    with assert_throws(AssertionError):
        # fails distributivity
        check_axioms_samples(
            Interval,
            [Interval(1,1),
             Interval(-1,1),
             Interval(1,2)],
            star=False
        )

    assert Interval(0,1) | Interval(3,4) == Interval(0,4)

    assert Interval.tightest([]) == Interval.loose
    assert Interval.tightest([Interval(0,1)]) == Interval(0,1)
    assert Interval.tightest([Interval(0,1), Interval(-1,0.5)]) == Interval(0,0.5)

    assert Interval.union([]) == Interval.loose
    assert Interval.union([Interval(0,1)]) == Interval(0,1)
    assert Interval.union([Interval(0,1), Interval(-1,0.5)]) == Interval(-1,1)

    # TODO: should this be some special empty set?
    assert Interval(0,1) & Interval(3,4) == Interval(3,1)


def test_logval():

    samples = [
        LogVal.lift(.1),
        LogVal.lift(.2),
        LogVal.lift(.3),
        LogVal.lift(-0.1),
        LogVal.lift(0),
        LogVal.lift(10),
        LogVal.lift(23),
        LogVal.lift(-3.1e-5),
        LogVal.lift(-123.5),
        LogVal.lift(1e-100),
        LogVal.lift(1e+100),
    ]


    from semirings.logval import log1pexp, log1mexp
    xs = np.linspace(-40, 40, 100)
    ys = [log1pexp(x) for x in xs]
    zs = np.log(1+np.exp(xs))
    np.allclose(ys, zs)
    if 0:
        import pylab as pl
        pl.plot(xs, ys, c='b', alpha=0.5, label='log1pexp')
        pl.plot(xs, zs, c='r', alpha=0.5, linestyle=':', label='log(1+exp(x))')
        pl.show()

    xs = np.linspace(-40, 0, 100)
    ys = [log1mexp(x) for x in xs]
    zs = np.log(1-np.exp(xs))
    np.allclose(ys, zs)
    if 0:
        import pylab as pl
        pl.plot(xs, ys, c='b', alpha=0.5, label='log1mexp')
        pl.plot(xs, zs, c='r', alpha=0.5, linestyle=':', label='log(1-exp(x))')
        pl.show()

    assert LogVal.lift(.1) < LogVal.lift(.2)
    assert LogVal.lift(-.2) < LogVal.lift(-.1)
    assert LogVal.zero / LogVal.one == LogVal.zero

    for a in samples:
        for b in samples:
            for c in samples:
                check_axioms(LogVal, a, b, c)


def test_axioms():

    for S,samples in [
            (Lukasiewicz, [
                Lukasiewicz.zero,
                Lukasiewicz.one,
                Lukasiewicz(0.5),
                Lukasiewicz(0.25),
                Lukasiewicz(0.75),
            ]),
            (MinPlus, [
                MinPlus(float('+inf')),
                MinPlus.zero,
                MinPlus.one,
                MinPlus(-3),
                MinPlus(+3),
                MinPlus(2),
                MinPlus(-1)
            ]),
            (MaxPlus, [
                MaxPlus(float('-inf')),
                MaxPlus.zero,
                MaxPlus.one,
                MaxPlus(-3),
                MaxPlus(+3),
                MaxPlus(2),
                MaxPlus(-1),
            ]),
            (MinTimes, [
                MinTimes(float('+inf')),
                MinTimes.zero,
                MinTimes.one,
                MinTimes(+3),
                MinTimes(2),
            ]),
            (Dual, [
                Dual(0,1),
                Dual(-3,2),
                Dual(2,4),
                Dual(-1,-1),
            ]),
            (Float, [
                0, -3, 2,
                1, -1,
            ]),
            (Boolean, [
                Boolean(True),
                Boolean(False)
            ]),
            (RegularLanguage, [
                RegularLanguage.lift('a'),
                RegularLanguage.lift('b'),
                RegularLanguage.lift('c'),
                RegularLanguage.zero,
                RegularLanguage.one,
            ])
            #Symbol,   # equality operation of regex is insufficient
            #LazySort,
    ]:

        check_axioms_samples(S,samples)


def check_axioms_samples(S,samples,**kwargs):
    for A in samples:
        for B in samples:
            for C in samples:
                check_axioms(S,A,B,C,**kwargs)


def check_axioms(S,A,B,C,star=True,left_distrib=True,right_distrib=True,hash_=True,assoc=True,check_multiplicity=True):
    try:

        assert A == A
        if hash_:
            assert hash(A) == hash(A)
#            assert hash((A + B) * C) == hash(A * C + B * C)

        # Addition is a commutative monoid
        assert_equal((A + B) + C,
                     A + (B + C))
        assert_equal((A + B),
                     (B + A))
        assert_equal(A + S.zero,
                     A,
                     S.zero + A)

        # Multiplication is a monoid
        if assoc:
            assert_equal((A * B) * C, A * (B * C))
            assert_equal(A * S.one, A, S.one * A)

        # distributivity from the left and the right
        if left_distrib:  assert_equal(A * (B + C), A * B + A * C)
        if right_distrib: assert_equal((B + C) * A, B * A + C * A)

        # Annihilation
        assert_equal(A * S.zero, S.zero, S.zero * A)

        # Multiplicity operation
        if check_multiplicity:
            assert_equal(S.multiplicity(A,0) == S.zero)
            assert_equal(S.multiplicity(A,1) == A)
            assert_equal(S.multiplicity(A,2) == A + A)
            assert_equal(S.multiplicity(A,3) == A + A + A)
            #assert_equal(S.multiplicity(A,np.inf) == S.zero)

        if star:
            # Kleene star operation
            assert_equal(S.star(A),
                         S.one + S.star(A) * A,
                         S.one + A * S.star(A))

    except AssertionError as e:
        print()
        print(colors.light.cyan % S.__name__,A,B,C)
        # check semiring axioms for a subset of values
        raise e


def assert_equal(x, *ys):  # pragma: no cover
    if all(x == y for y in ys):
        #print(colors.ok, x, *ys)
        pass
    else:
        print(colors.bad, x, *ys)
        assert 0


if __name__ == '__main__':
    from arsenal import testing_framework
    testing_framework(globals())
