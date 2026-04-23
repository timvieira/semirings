import numpy as np
import pytest
from arsenal import assert_throws
from arsenal.iterextras import take
from semirings import (
    MinPlus, MinTimes, MaxPlus, MaxTimes, Float, CutSets, Boolean,
    Bottleneck, minmax, maxmin, LogVal, ConvexHull, Point,
    Lukasiewicz, Interval, LazySort, Dual, Bridge, Division,
    make_set, String, ThreeValuedLogic, VBridge, Why, Lineage, make_semiring, MatrixSemiring, Entropy,
    check_axioms_samples, check_axioms, check_metric_axioms,
    WeightedGraph,
)
from semirings.regex import RegularLanguage
from fsa import FSA


# 8-edge graph used by the provenance / bridge tests (Why, Lineage, Bridge,
# VBridge). Has a bridge at 2→3 and two parallel paths 0→1→2 / 0→1.5→2.
PROVENANCE_EDGES = [
    (0, 1), (2, 3), (1, 2), (0, 1.5), (1.5, 2),
    (3, 4), (3, 5), (4, 5),
]


def provenance_graph(WeightType):
    g = WeightedGraph(WeightType)
    for i, j in PROVENANCE_EDGES:
        g.edge(i, j)
    return g


def assert_star_approx_matches(S, members, T, *, inf=None):
    """Verify S.star_approx(x, T) == S.star(x) for all members; skip divergent
    cases where the analytical star equals the supplied `inf` value."""
    for x in members:
        analytical = S.star(x)
        if inf is not None and analytical == inf:
            continue
        S.assert_equal(S.star_approx(x, T), analytical)


def test_pow():
    assert Float.lift(0, None) ** 0 == 1
    assert Float.lift(1, None) ** 10 == 1
    assert Float.lift(2, None) ** 10 == 2**10
    assert Float.lift(0.5, None) ** 12 == 0.5**12

    assert (MaxPlus.lift(0.5, None) ** 12).score == 0.5 * 12
    MaxPlus.assert_equal(MaxPlus.lift(2, None) ** 0, MaxPlus.one)

    MaxTimes.assert_equal(MaxTimes.lift(2, None) ** 5, MaxTimes.lift(2 ** 5, None))


# Pareto-front semiring sketch moved to experimental/pareto.py.

# `test_uncertainty_sets` moved to experimental/uncertainty.py —
# UncertaintySet(base) is a sub-distributive meta-semiring, not a semiring.


# `test_funky` moved to experimental/funky.py — notes on relaxing * from a
# binary monoid to a K-arity operator (Gilea 2020), not an actual test.


# `test_endomorphism` moved to experimental/endomorphism.py — the sketch used
# affine (not linear) members, which aren't endomorphisms of (ℝ, +), so it
# never formed a proper semiring.


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
        S.assert_equal(S.star_approx(x, 100), S.star(x))

    check_axioms_samples(S,members)


# `test_semilinear` moved to experimental/semilinear.py — the Semilinear/Parikh
# semiring needs a symbolic representation of star before it can be tested.


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

    for x in members:
        S.assert_equal(S.star(x), S.star_fixpoint(x))

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

    assert_star_approx_matches(S, members, 100, inf=S.inf)
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

    assert_star_approx_matches(S, members, 1000, inf=S(np.inf))
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

    assert_star_approx_matches(S, members, 100, inf=S.inf)
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
        S.assert_equal(S.star(x), S.star_fixpoint(x))

    check_axioms_samples(S,members)


def test_why_semiring():

    members = [
        Why.lift('a'),
        Why.lift('b'),
        Why.lift('c'),
    ]
    check_axioms_samples(Why,members)

    graph = provenance_graph(Why)
    why = graph.closure()[0,5]   # we have a bridge at edge 2->3
    Why.assert_equal(why, Why(frozenset({
        frozenset({(0, 1), (2, 3), (4, 5), (1, 2), (3, 4)}),
        frozenset({(0, 1.5), (2, 3), (4, 5), (3, 4), (1.5, 2)}),
        frozenset({(0, 1), (1, 2), (2, 3), (3, 5)}),
        frozenset({(0, 1.5), (2, 3), (1.5, 2), (3, 5)})
    })))


def test_lineage_semiring():

    members = [
        Lineage.lift('a'),
        Lineage.lift('b'),
        Lineage.lift('c'),
    ]

    check_axioms_samples(Lineage,members)

    graph = provenance_graph(Lineage)
    c = graph.closure()

    # we use all the edges for 0->5 paths
    lineage = c[0,5]
    Lineage.assert_equal(lineage, Lineage(frozenset(graph.E)))

    # we only use four edges for 0->2 paths
    lineage = c[0,2]
    Lineage.assert_equal(lineage, Lineage(frozenset({(0, 1), (1, 2), (1.5, 2), (0, 1.5)})))


def test_cut_sets():
    S = CutSets

    # minimality
    CutSets.assert_equal(CutSets({(2,3)}, {(1,3), (2,3)}), CutSets({(2,3)}))

    a = CutSets( {(2,3)}, {(1,3), (2,4)} )
    b = CutSets( {(1,3), (2,3)}, {(1,3), (2,4)} )
    c = CutSets( {(1,3)} )
    CutSets.assert_equal(a * b, CutSets({(2,3)}, {(1,3), (2,4)}))
    CutSets.assert_equal(a + b, CutSets({(1,3), (2,3)}, {(1,3), (2,4)}))

    one = CutSets.one; zero = CutSets.zero
    CutSets.assert_equal(a * one, a)
    CutSets.assert_equal(a * zero, zero)
    CutSets.assert_equal(a * (b + c), a * b + a * c)
    CutSets.assert_equal(a * a, a)
    CutSets.assert_equal(a + a, a)

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

    graph = provenance_graph(Bridge)
    bridges = graph.closure()[0,5]   # we have a bridge at edge 2->3
    Bridge.assert_equal(bridges, Bridge({(2,3)}))


def test_vertex_bridge():
    graph = provenance_graph(VBridge)
    # clearly, removing 0 or 5 disconnects 0 from 5; but so will 2 or 3 (see picture above)
    b = graph.closure()[0,5]

    VBridge.assert_equal(b, VBridge({0,2,3,5}))


def test_three_valued_logic():
    S = ThreeValuedLogic

    [f,u,t] = members = S.zero, S.unk, S.one

    assert f < u < t
    S.assert_equal(f + f, f)
    S.assert_equal(f + t, t)
    S.assert_equal(t + t, t)
    S.assert_equal(f + u, u)
    S.assert_equal(t + u, t)
    S.assert_equal(u + u, u)

    S.assert_equal(f * f, f)
    S.assert_equal(f * t, f)
    S.assert_equal(t * t, t)
    S.assert_equal(f * u, f)
    S.assert_equal(t * u, u)
    S.assert_equal(u * u, u)

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

    S.assert_equal(S("aaa") + S("aab"), S("aa"))
    S.assert_equal(S("aaa") + S(""), S(""))

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


def test_lazysort_sorted_invariant():
    """Verify that LazySort iterators always yield elements in descending score order."""

    def scores(expr, T=20):
        return [x.score for x in take(T, expr)]

    def is_descending(xs):
        return all(a >= b for a, b in zip(xs, xs[1:]))

    a = LazySort(0.5, 'a')
    b = LazySort(0.3, 'b')
    c = LazySort(0.7, 'c')

    # Sum should be sorted
    s = scores(a + b + c)
    assert is_descending(s), f'Sum not sorted: {s}'

    # Prod should be sorted
    s = scores((a + b) * (b + c))
    assert is_descending(s), f'Prod not sorted: {s}'

    # Nested expression
    s = scores((a + b) * (b + c) + a * c)
    assert is_descending(s), f'Nested expr not sorted: {s}'


def test_lazysort_star_multi_element():
    """Star of a multi-element sum must yield in sorted order.

    This is the critical test: Star.__iter__ currently concatenates
    `yield from` on each power rather than merging sorted streams,
    so the output is NOT globally sorted when the argument is a sum.

    Example: star(0.5 + 0.3) yields
      [1, 0.5, 0.3, 0.25, 0.15, 0.15, 0.09, 0.125, ...]
                                              ^^^^^ ^^^^^ BUG: 0.125 > 0.09
    The power-2 terms (0.25, 0.15, 0.15, 0.09) are interleaved correctly
    within that power, but the transition to power-3 (0.125) is not merged
    with the tail of power-2.
    """

    def scores(expr, T=15):
        return [x.score for x in take(T, expr)]

    def is_descending(xs):
        return all(a >= b for a, b in zip(xs, xs[1:]))

    a = LazySort(0.5, 'a')
    b = LazySort(0.3, 'b')
    ab = a + b

    s = scores(ab.star())
    assert is_descending(s), (
        f'Star of multi-element sum is not sorted:\n{s}'
    )


def test_lazysort_star_single_descending():
    """Star of a single element with score in (0,1) should be sorted descending."""

    def scores(expr, T=10):
        return [x.score for x in take(T, expr)]

    def is_descending(xs):
        return all(a >= b for a, b in zip(xs, xs[1:]))

    a = LazySort(0.5, 'a')
    s = scores(a.star())
    # powers: 1, 0.5, 0.25, 0.125, ... — should be descending
    assert is_descending(s), f'Star of single element (score<1) not sorted: {s}'


def test_lazysort_star_score_lt_1():
    """Star of an element with score < 1 yields descending scores."""

    def scores(expr, T=10):
        return [x.score for x in take(T, expr)]

    def is_descending(xs):
        return all(a >= b for a, b in zip(xs, xs[1:]))

    a = LazySort(0.5, 'a')
    s = scores(a.star(), T=5)
    # scores are [1, 0.5, 0.25, 0.125, 0.0625] — descending
    assert is_descending(s), f'Star of element with score<1 not sorted: {s}'


def test_lazysort_star_composed():
    """Star result used inside a larger expression must maintain sort order."""

    def scores(expr, T=20):
        return [x.score for x in take(T, expr)]

    def is_descending(xs):
        return all(a >= b for a, b in zip(xs, xs[1:]))

    a = LazySort(0.5, 'a')
    b = LazySort(0.3, 'b')
    c = LazySort(0.8, 'c')

    # star(a) * c — should still be sorted
    expr = a.star() * c
    s = scores(expr)
    assert is_descending(s), f'star(a) * c not sorted: {s}'

    # c + star(a+b) — should be sorted
    expr = c + (a + b).star()
    s = scores(expr)
    assert is_descending(s), f'c + star(a+b) not sorted: {s}'


def test_lazysort_distributivity():
    """Test distributivity: a*(b+c) == a*b + a*c for various expressions."""

    def check(a, b, T=None):
        A = list(a.take(T))
        B = list(take(T, b))
        assert A == B, f'\n\n{A}\n{B}\n'

    a = LazySort(0.5, 'a')
    b = LazySort(0.3, 'b')
    c = LazySort(0.7, 'c')
    d = LazySort(0.2, 'd')

    # basic distributivity
    check(a * (b + c), a * b + a * c)

    # right distributivity
    check((a + b) * c, a * c + b * c)

    # with more terms
    check(a * (b + c + d), a * b + a * c + a * d)

    # double distribution
    check((a + b) * (c + d), a*c + a*d + b*c + b*d)


def test_lazysort_associativity_of_sum():
    """(a + b) + c should yield same elements as a + (b + c)."""

    def check(a, b, T=None):
        A = list(a.take(T))
        B = list(take(T, b))
        assert A == B, f'\n\n{A}\n{B}\n'

    a = LazySort(0.5, 'a')
    b = LazySort(0.3, 'b')
    c = LazySort(0.7, 'c')

    check((a + b) + c, a + (b + c))


def test_bottleneck():

    one = Bottleneck.one
    zero = Bottleneck.zero

    a = Bottleneck(1)
    b = Bottleneck(2)
    c = Bottleneck(3)
    d = Bottleneck(4)

    Bottleneck.assert_equal(a.star_approx(10), a.star())
    Bottleneck.assert_equal(a.star_fixpoint(), a.star())

    Bottleneck.assert_equal(Bottleneck.sum([]), zero)
    Bottleneck.assert_equal(Bottleneck.sum([a,b,c,d]), a + b + c + d)
    Bottleneck.assert_equal(Bottleneck.product([]), one)
    Bottleneck.assert_equal(Bottleneck.product([a,b,c,d]), a * b * c * d)

    Bottleneck.assert_equal(a * b + b * c + d * c, c)
    Bottleneck.assert_equal(a + b, b)
    Bottleneck.assert_equal(a * b, a)
    Bottleneck.assert_equal(a * one, a)
    Bottleneck.assert_equal(a * zero, zero)
    Bottleneck.assert_equal(a + zero, a)
    Bottleneck.assert_equal(b + one, one)


def test_convex_hull():
    check_axioms(
        ConvexHull,
        ConvexHull([Point(1,1.5), Point(2,1.75), Point(3,4)]),
        ConvexHull([Point(2,2)]),
        ConvexHull([Point(1,10)]),
        star=False,
    )


def test_star_operations():
    Float.assert_equal(Float.star(1), np.inf)
    MinTimes.assert_equal(MinTimes(-1).star(), MinTimes(-np.inf))
    MinPlus.assert_equal(MinPlus(-1).star(), MinPlus(-np.inf))
    MaxPlus.assert_equal(MaxPlus(1).star(), MaxPlus(np.inf))


# Intervals are sub-distributive rather than distributive, so
# check_axioms_samples fails on a generic member set — the test explicitly
# asserts that failure below, alongside positive checks on a curated triple.
def test_interval():

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

    Interval.assert_equal(Interval(0,1) | Interval(3,4), Interval(0,4))

    Interval.assert_equal(Interval.tightest([]), Interval.loose)
    Interval.assert_equal(Interval.tightest([Interval(0,1)]), Interval(0,1))
    Interval.assert_equal(Interval.tightest([Interval(0,1), Interval(-1,0.5)]), Interval(0,0.5))

    Interval.assert_equal(Interval.union([]), Interval.loose)
    Interval.assert_equal(Interval.union([Interval(0,1)]), Interval(0,1))
    Interval.assert_equal(Interval.union([Interval(0,1), Interval(-1,0.5)]), Interval(-1,1))

    # TODO: should this be some special empty set?
    Interval.assert_equal(Interval(0,1) & Interval(3,4), Interval(3,1))


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
    assert np.allclose(ys, zs, equal_nan=True)

    # Exclude the x=0 endpoint: naive log(1-exp(0)) = log(0) = -inf, and the
    # sub-epsilon region near 0 is exactly where naive is catastrophically
    # cancelling — this test only validates agreement in the benign range.
    xs = np.linspace(-40, 0, 100, endpoint=False)
    ys = [log1mexp(x) for x in xs]
    zs = np.log(1-np.exp(xs))
    assert np.allclose(ys, zs)

    assert LogVal.lift(.1) < LogVal.lift(.2)
    assert LogVal.lift(-.2) < LogVal.lift(-.1)
    LogVal.assert_equal(LogVal.zero / LogVal.one, LogVal.zero)

    check_metric_axioms(LogVal, samples)

    check_axioms_samples(LogVal, samples)


def _build_set_semiring():
    from arsenal.maths.combinatorics import powerset
    universe = frozenset('abc')
    Set = make_set(universe, frozenset())
    members = [Set(frozenset(x)) for x in powerset(universe)]
    return Set, members


_MinMax = minmax(None, "")
_MaxMin = maxmin(None, "")
_Set, _set_members = _build_set_semiring()


# (semiring, samples, kwargs) — each entry becomes test_axioms[<id>] under pytest.
AXIOM_CASES = [
    (Lukasiewicz, [
        Lukasiewicz.zero, Lukasiewicz.one,
        Lukasiewicz(0.5), Lukasiewicz(0.25), Lukasiewicz(0.75),
    ], {}),
    (MinPlus, [
        MinPlus(float('+inf')), MinPlus.zero, MinPlus.one,
        MinPlus(-3), MinPlus(+3), MinPlus(2), MinPlus(-1),
    ], {}),
    (MaxPlus, [
        MaxPlus(float('-inf')), MaxPlus.zero, MaxPlus.one,
        MaxPlus(-3), MaxPlus(+3), MaxPlus(2), MaxPlus(-1),
    ], {}),
    (MinTimes, [
        MinTimes(float('+inf')), MinTimes.zero, MinTimes.one,
        MinTimes(+3), MinTimes(2),
    ], {}),
    (Entropy, [
        Entropy.zero, Entropy.one,
        Entropy.lift(0.3), Entropy.lift(0.5), Entropy.lift(0.7),
    ], {}),
    (Dual, [Dual(0,1), Dual(-3,2), Dual(2,4), Dual(-1,-1)], {}),
    (Float, [0, -3, 2, 1, -1], {}),
    (Boolean, [Boolean(True), Boolean(False)], {}),
    (RegularLanguage, [
        RegularLanguage.lift('a'),
        RegularLanguage.lift('b'),
        RegularLanguage.lift('c'),
        RegularLanguage.zero,
        RegularLanguage.one,
    ], {}),
    (Division, [Division(x) for x in range(10)], {}),
    (_Set, _set_members, {}),
    (_MinMax, list(map(_MinMax, "abcd")) + [_MinMax.zero, _MinMax.one], {}),
    (_MaxMin, list(map(_MaxMin, "abcd")), {}),
    (MaxTimes, [
        MaxTimes.zero, MaxTimes.one,
        MaxTimes(0.5), MaxTimes(0.1), MaxTimes(0),
    ], {}),
    (Bottleneck, [
        Bottleneck.zero, Bottleneck.one,
        Bottleneck(1), Bottleneck(2), Bottleneck(3), Bottleneck(4),
    ], {}),
    (Why, [Why.zero, Why.one, Why.lift('a'), Why.lift('b'), Why.lift('c')], {}),
    (Lineage, [
        Lineage.zero, Lineage.one,
        Lineage.lift('a'), Lineage.lift('b'), Lineage.lift('c'),
    ], {}),
    (Bridge, [
        Bridge.zero, Bridge.one,
        Bridge({(2,3), (1,3), (2,4)}),
        Bridge({(1,3), (2,3)}),
        Bridge({(1,3), (2,4)}),
        Bridge.lift((1,3)),
    ], {}),
    (CutSets, [
        CutSets.zero, CutSets.one,
        CutSets({(2,3)}, {(1,3), (2,4)}),
        CutSets({(1,3), (2,3)}, {(1,3), (2,4)}),
        CutSets({(1,3)}),
    ], {}),
    (ThreeValuedLogic, [ThreeValuedLogic.zero, ThreeValuedLogic.unk, ThreeValuedLogic.one], {}),
    (String, [
        String.zero, String.one,
        String("a"), String("aa"), String("aaa"),
        String("b"), String("aaab"), String("abab"),
    ], {'right_distrib': False}),
    # Symbol, LazySort: regex equality is insufficient / lazy.
    # VBridge: same axiom shape as Bridge; tested via test_vertex_bridge.
]


@pytest.mark.parametrize(
    'S,samples,kwargs',
    AXIOM_CASES,
    ids=[case[0].__name__ for case in AXIOM_CASES],
)
def test_axioms(S, samples, kwargs):
    check_metric_axioms(S, samples)
    check_axioms_samples(S, samples, **kwargs)


def test_metric_extended_ranges():
    # Metric should handle values outside the semiring's normal range.
    check_metric_axioms(Float, [0, -3, 2, 1, -1, float('inf'), float('-inf'), 1e100, -1e100])
    check_metric_axioms(MinTimes, [MinTimes(x) for x in [0, 1, 3, float('inf'), 1e100]])
    check_metric_axioms(MaxTimes, [MaxTimes(x) for x in [0, 1, 3, float('inf'), 1e100]])


def test_matrix_semiring_boolean():
    M = MatrixSemiring(Boolean, 'ab')
    T, F = Boolean(True), Boolean(False)

    members = [
        M.zero,
        M.one,
        M({('a','a'): T, ('a','b'): T, ('b','a'): F, ('b','b'): F}),
        M({('a','a'): F, ('a','b'): T, ('b','a'): T, ('b','b'): F}),
    ]

    check_metric_axioms(M, members)
    check_axioms_samples(M, members)


def test_matrix_semiring_minplus():
    M = MatrixSemiring(MinPlus, [0, 1])

    members = [
        M.zero,
        M.one,
        M({(0,0): MinPlus(1), (0,1): MinPlus(2),
           (1,0): MinPlus(3), (1,1): MinPlus(4)}),
        M({(0,0): MinPlus(0), (0,1): MinPlus.zero,
           (1,0): MinPlus.zero, (1,1): MinPlus(0)}),
    ]

    check_metric_axioms(M, members)
    check_axioms_samples(M, members)


def test_matrix_semiring_logval():
    M = MatrixSemiring(LogVal, [0, 1])

    members = [
        M.zero,
        M.one,
        M({(0,0): LogVal.lift(0.3), (0,1): LogVal.lift(0.1),
           (1,0): LogVal.lift(0.2), (1,1): LogVal.lift(0.4)}),
        M({(0,0): LogVal.lift(0.5), (0,1): LogVal.lift(0.0),
           (1,0): LogVal.lift(0.0), (1,1): LogVal.lift(0.5)}),
    ]

    check_metric_axioms(M, members)
    check_axioms_samples(M, members, star=False)


def test_star_doubling():
    # star_doubling must agree with star_fixpoint wherever the series converges.
    # Use each semiring's own metric via assert_equal rather than ad-hoc tolerances.

    def agree(x, *, expect=None, tol=1e-10):
        # Fixpoint's convergence check is itself metric-based, so its precision
        # is bounded by the semiring's own equality tolerance; callers can loosen
        # `tol` for floating-point semirings where that matters.
        S = type(x)
        fp = x.star_fixpoint()
        db = x.star_doubling()
        if expect is None:
            S.assert_equal(fp, db, tol=tol)
        else:
            S.assert_equal(expect, fp, db, tol=tol)
        return fp, db

    # LogVal: geometric series 1/(1-x), converges when x < 1.
    # Fixpoint convergence is ~1e-5 relative due to np.allclose in LogVal.__eq__.
    for v in [0.0, 0.1, 0.5, 0.9]:
        agree(LogVal.lift(v), expect=LogVal.lift(1/(1 - v)), tol=1e-3)

    # Boolean: idempotent; star(x) = one.
    for x in [Boolean.zero, Boolean.one]:
        agree(x, expect=Boolean.one)

    # MaxPlus: tropical; star converges when score <= 0.
    for v in [-5.0, -1.0, -0.1, 0.0]:
        agree(MaxPlus(v), expect=MaxPlus.one)

    # MinPlus: tropical; star converges when cost >= 0.
    for v in [0.0, 0.1, 1.0, 5.0]:
        agree(MinPlus(v), expect=MinPlus.one)

    # Łukasiewicz: + = max, so one absorbs — star(x) = one.
    for v in [0.0, 0.25, 0.5, 0.9, 1.0]:
        agree(Lukasiewicz(v), expect=Lukasiewicz.one)

    # Three-valued logic: + = max, * = min; star(x) = one.
    for x in [ThreeValuedLogic.zero, ThreeValuedLogic.unk, ThreeValuedLogic.one]:
        agree(x, expect=ThreeValuedLogic.one)

    # Why-provenance: idempotent set-of-sets.
    for x in [Why.zero, Why.one, Why.lift('a'), Why.lift('a') + Why.lift('b')]:
        agree(x)

    # Lineage: idempotent set semiring.
    for x in [Lineage.zero, Lineage.one, Lineage.lift('a'), Lineage.lift('a') + Lineage.lift('b')]:
        agree(x)

    # Bridge: intersection/union set semiring, idempotent.
    for x in [Bridge.zero, Bridge.one, Bridge.lift((0,1)), Bridge.lift((0,1)) * Bridge.lift((1,2))]:
        agree(x)

    # Division (gcd/lcm): gcd(1, _) = 1, so star(x) = Division(1) = one.
    for v in [0, 1, 2, 6, 15]:
        agree(Division(v), expect=Division.one)

    # Nilpotent matrix over LogVal: star(N) = I + N + N² terminates finitely.
    M = MatrixSemiring(LogVal, range(3))
    N = M({(0,1): LogVal.lift(0.5), (1,2): LogVal.lift(0.25), (0,2): LogVal.lift(0.1)})
    M.assert_equal(N.star_fixpoint(), N.star_doubling(), tol=1e-3)


def test_weighted_graph_transpose():
    g = WeightedGraph(Float)
    g[0, 1] = 0.3
    g[1, 2] = 0.5
    g[2, 0] = 0.7

    T = g.transpose()
    assert set(T.E) == {(1, 0), (2, 1), (0, 2)}
    for (i, j), w in g.E.items():
        Float.assert_equal(T[j, i], w)

    # Transpose-of-transpose recovers the original.
    TT = T.transpose()
    assert set(TT.E) == set(g.E)
    for k in g.E:
        Float.assert_equal(TT.E[k], g.E[k])

    # Lazy `incoming` matches a direct reverse-adjacency build.
    assert dict(g.incoming) == {1: {0}, 2: {1}, 0: {2}}


def test_weighted_graph_scc():
    # SCC-based closure should agree with the naive Floyd–Warshall reference.
    # Build a graph with two cycles connected by a bridge — nontrivial SCCs.
    g = WeightedGraph(Float)
    for (i, j) in [(0,1), (1,0),               # SCC {0,1}
                   (1,2),                       # bridge
                   (2,3), (3,4), (4,2),         # SCC {2,3,4}
                   (4,5)]:                      # sink {5}
        g[i, j] = 0.3

    scc_based = g.closure_scc_based()
    reference = g.closure_reference()
    for key in set(scc_based) | set(reference):
        Float.assert_equal(scc_based.get(key, Float.zero),
                           reference.get(key, Float.zero),
                           tol=1e-9)

    # solve_left: x = xA + b with b = e_0 should recover the 0-th row of star(A).
    b = Float.chart(); b[0] = Float.one
    sol = g.solve_left(b)
    for j in g.N:
        Float.assert_equal(sol.get(j, Float.zero),
                           scc_based.get((0, j), Float.zero),
                           tol=1e-9)

    # solve_right: x = Ax + b with b = e_5 should recover the 5-th column.
    b = Float.chart(); b[5] = Float.one
    sol = g.solve_right(b)
    for i in g.N:
        Float.assert_equal(sol.get(i, Float.zero),
                           scc_based.get((i, 5), Float.zero),
                           tol=1e-9)


if __name__ == '__main__':
    from arsenal import testing_framework
    testing_framework(globals())
