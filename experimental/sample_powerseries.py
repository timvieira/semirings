"""Truncated power series over an arbitrary scalar semiring.

This module provides the *generic* list-backed truncated-polynomial semiring
factory ``PowerSeriesOver(scalar, max_degree)``. Unlike the numpy-backed
``TruncatedPowerSeries`` (which is hardcoded to ``float`` scalars for FFT /
BLAS speed), this one accepts any ``Semiring`` scalar class and performs
multiplication in pure Python — O(N^2) — using the scalar's own ``+`` and
``*``. That buys us the full semiring-parsing generality: set the scalar
to whatever captures the information you want (samples, entropies,
derivations, provenance, …) and the DP automatically produces the
length-indexed refinement.

The canonical application in this notebook series is
``PowerSeriesOver(Sample, N)`` — truncated power series with
``Sample`` coefficients — which turns ``CFG.treesum()`` into an exact
length-conditional sampler: iterating the returned ``result.c[n]``
draws derivations of yield length ``n`` proportional to their PCFG
probability.

``SamplePowerSeries(N)`` is a convenience alias that supplies a
Sample-aware coefficient metric (``|a.w − b.w|`` per coefficient) so that
``CFG.agenda()`` converges even though ``Sample`` inherits the default
discrete ``==`` / ``metric``.

---

Why separate from ``TruncatedPowerSeries``:

- Different coefficient representation (Python list of scalar-semiring
  elements vs. numpy array of floats) and different multiplication
  machinery (Python nested loop in scalar arithmetic vs. numpy/BLAS/FFT).
  The numeric path optimizes for a specific common case; the generic
  path trades speed for generality.
- ``zero`` / ``one`` are scalar-valued here (``scalar.zero``,
  ``scalar.one``) rather than numeric primitives.

Both live under the same parent abstraction — commutative polynomial
semiring over a base ring — and could share a skeleton via a mixin if a
third variant arrives.
"""

from semirings.sample import Sample
from semirings.base import Semiring


_factory_cache = {}


def PowerSeriesOver(scalar, max_degree, *, coeff_metric=None):
    """Return a Semiring class for ``scalar[x]/(x^{max_degree+1})``.

    Parameters
    ----------
    scalar : a ``Semiring`` subclass
        The coefficient ring. Must provide ``zero``, ``one``, ``+``, ``*``.
        If the agenda-based treesum is to converge on your DP, the scalar
        should also provide a meaningful ``metric`` — or pass
        ``coeff_metric`` explicitly.
    max_degree : int
        Truncation degree ``N``. Elements hold coefficients for
        ``x^0 .. x^N``.
    coeff_metric : callable (a, b) -> float, optional
        How to measure "distance between two coefficients" for the
        purposes of agenda convergence. Default uses ``a.metric(b)`` if
        both are not the scalar's zero sentinel, else a discrete
        fallback. Override for scalars (like ``Sample``) whose default
        metric is discrete but whose "weight" field provides a
        meaningful numeric distance.
    """
    N = int(max_degree)
    if N < 0:
        raise ValueError(f'max_degree must be >= 0, got {max_degree}')

    key = (scalar, N, coeff_metric)
    cached = _factory_cache.get(key)
    if cached is not None:
        return cached

    _scalar_zero = scalar.zero
    _scalar_one = scalar.one

    def _default_coeff_metric(a, b):
        if a is b:
            return 0.0
        if a is _scalar_zero:
            return 1.0 if b is not _scalar_zero else 0.0
        if b is _scalar_zero:
            return 1.0
        if hasattr(a, 'metric'):
            return float(a.metric(b))
        return 0.0 if a == b else 1.0

    _coeff_metric = coeff_metric if coeff_metric is not None else _default_coeff_metric

    class _TPS(Semiring):
        max_degree = N
        scalar_type = scalar

        def __init__(self, coeffs):
            coeffs = list(coeffs)
            if len(coeffs) != N + 1:
                raise ValueError(
                    f'coeffs must be length {N + 1}, got {len(coeffs)}'
                )
            self.c = coeffs

        # ---------------------------- constructors
        @classmethod
        def monomial(cls, scalar_elem, k):
            "``scalar_elem · x^k`` (or ``zero`` if k > N)."
            c = [_scalar_zero] * (N + 1)
            if 0 <= k <= N:
                c[k] = scalar_elem
            return cls(c)

        @classmethod
        def constant(cls, scalar_elem):
            "``scalar_elem · x^0``."
            return cls.monomial(scalar_elem, 0)

        # ---------------------------- arithmetic
        def __add__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            out = []
            for a, b in zip(self.c, other.c):
                if a is _scalar_zero:
                    out.append(b)
                elif b is _scalar_zero:
                    out.append(a)
                else:
                    out.append(a + b)
            return _TPS(out)

        def __mul__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            out = [_scalar_zero] * (N + 1)
            for i, ai in enumerate(self.c):
                if ai is _scalar_zero:
                    continue
                for j, bj in enumerate(other.c):
                    if bj is _scalar_zero:
                        continue
                    k = i + j
                    if k > N:
                        break
                    out[k] = out[k] + ai * bj
            return _TPS(out)

        # ---------------------------- metric / equality
        def metric(self, other):
            if not isinstance(other, _TPS):
                return float('inf')
            return float(sum(
                _coeff_metric(a, b) for a, b in zip(self.c, other.c)
            ))

        def __eq__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            return self.metric(other) == 0.0

        def __hash__(self):
            # Any hash compatible with __eq__; metric==0 is a weaker relation
            # than "same instances" so we use the by-index tuple of ids/scalar
            # hashes when available, otherwise a constant. Agenda doesn't need
            # hashing; this is only present for dict keys in downstream code.
            try:
                return hash(tuple(hash(a) for a in self.c))
            except TypeError:
                return hash(id(self))

        def __repr__(self):
            return f'TPS[{scalar.__name__},N={N}]({self.c!r})'

    _TPS.zero = _TPS([_scalar_zero] * (N + 1))
    one_coeffs = [_scalar_zero] * (N + 1)
    one_coeffs[0] = _scalar_one
    _TPS.one = _TPS(one_coeffs)
    _TPS.__name__ = f'PowerSeriesOver_{scalar.__name__}_N{N}'
    _factory_cache[key] = _TPS
    return _TPS


# ---------------------------------------------------------------------------
# Sample-specific convenience.
# ---------------------------------------------------------------------------
def _sample_coeff_metric(a, b):
    """Weight-based metric for Sample coefficients.

    Sample inherits discrete ``metric``/``__eq__`` from Semiring/object; the
    weight field is the right numeric handle for agenda convergence.
    """
    wa = 0.0 if a is Sample.zero else a.w
    wb = 0.0 if b is Sample.zero else b.w
    return abs(wa - wb)


def SamplePowerSeries(max_degree):
    """``PowerSeriesOver(Sample, ...)`` with a weight-based coefficient metric.

    Iterating ``sps_element.c[n]`` draws exact samples from the conditional
    distribution over derivations whose polynomial index (e.g. yield length)
    equals ``n``.
    """
    return PowerSeriesOver(
        Sample, max_degree, coeff_metric=_sample_coeff_metric
    )


# ---------------------------------------------------------------------------
# CFG bridge (Sample-specific, since it constructs Sample rule-weights).
# ---------------------------------------------------------------------------
def build_sample_cfg(cfg, feature=None, max_value=None, *, use_cnf=True):
    """Wrap a Float-weighted CFG into a SamplePowerSeries-weighted CFG.

    Converts to CNF first (default) so that terminals only appear in
    preterminal rules — sampling otherwise drops positional terminals
    inside mixed-arity rule bodies (e.g. ``B -> A c``), because ``Prod``
    has no slot for a terminal between two non-terminal Samples.

    feature    : callable Rule -> int. Defaults to "number of terminals
                 in the body" (= yield length).
    max_value  : truncation degree; required.
    use_cnf    : bypass CNF conversion (only correct if the input grammar
                 is already CNF-like).
    """
    source = cfg.cnf if use_cnf else cfg
    if feature is None:
        feature = lambda r: sum(source.is_terminal(y) for y in r.body)
    if max_value is None:
        raise ValueError('max_value is required')
    SPS = SamplePowerSeries(max_value)
    new = source.__class__(R=SPS, S=source.S, V=source.V)
    for r in source.rules:
        k = feature(r)
        if len(r.body) == 1 and source.is_terminal(r.body[0]):
            data = (r.body[0],)
        else:
            data = ()
        new.add(SPS.monomial(Sample(r.w, data), k), r.head, *r.body)
    return new


def flatten_yield(d):
    """Walk a sampled ``d`` tree and return the yield as a tuple of leaves."""
    if isinstance(d, str):
        return (d,)
    if isinstance(d, (list, tuple)):
        out = []
        for x in d:
            out.extend(flatten_yield(x))
        return tuple(out)
    if d is None:
        return ()
    return (d,)


# ===========================================================================
# Tests
# ===========================================================================
import numpy as np
from semirings.base import check_axioms, check_metric_axioms
from semirings import Entropy


def _samples_for(TPS, scalar, rng, ctor=lambda w: w):
    """Small sample set for axiom checks of a generic TPS."""
    N = TPS.max_degree
    samples = [TPS.zero, TPS.one]
    for _ in range(3):
        k = int(rng.integers(0, N + 1))
        w = float(rng.random()) + 1e-3
        samples.append(TPS.monomial(ctor(w), k))
    return samples


def test_generic_tps_axioms_over_sample():
    "Truncated polynomial over Sample — axioms hold in the weight metric."
    rng = np.random.default_rng(0)
    for N in [0, 1, 4]:
        TPS = PowerSeriesOver(Sample, N, coeff_metric=_sample_coeff_metric)
        samples = _samples_for(TPS, Sample, rng, ctor=lambda w: Sample(w, (f'x',)))
        for A in samples:
            for B in samples:
                for C in samples:
                    check_axioms(TPS, A, B, C, star=False, hash_=False,
                                 check_multiplicity=False)


def test_generic_tps_axioms_over_entropy():
    "Truncated polynomial over Entropy — the generic factory is not Sample-specific."
    rng = np.random.default_rng(10)
    for N in [0, 1, 4]:
        TPS = PowerSeriesOver(Entropy, N)
        samples = _samples_for(TPS, Entropy, rng, ctor=Entropy.lift)
        for A in samples:
            for B in samples:
                for C in samples:
                    check_axioms(TPS, A, B, C, star=False, hash_=False,
                                 check_multiplicity=False)


def test_generic_tps_metric_axioms():
    rng = np.random.default_rng(1)
    TPS = PowerSeriesOver(Sample, 3, coeff_metric=_sample_coeff_metric)
    samples = _samples_for(TPS, Sample, rng, ctor=lambda w: Sample(w, (f'x',)))
    check_metric_axioms(TPS, samples)


def test_sample_power_series_is_alias():
    "SamplePowerSeries is just PowerSeriesOver(Sample, ...) with Sample's metric."
    A = SamplePowerSeries(4)
    B = PowerSeriesOver(Sample, 4, coeff_metric=_sample_coeff_metric)
    assert A is B
    # Different scalar => different class
    C = PowerSeriesOver(Entropy, 4)
    assert A is not C


def test_flatten_yield():
    assert flatten_yield(('a',)) == ('a',)
    assert flatten_yield([('a',), ('b',)]) == ('a', 'b')
    assert flatten_yield([[('a',), ('b',)], ('c',)]) == ('a', 'b', 'c')
    assert flatten_yield(()) == ()


def test_sample_power_series_cfg_end_to_end():
    """Exact sampling on a PCFG with multiple yields per length.

    Grammar:
      0.5: S -> A A
      0.5: S -> a
      0.5: A -> a
      0.5: A -> b

    At length 2, yields (aa, ab, ba, bb) are equiprobable.
    """
    import sys
    sys.path.insert(0, '/home/timv/projects/genparse')
    from genparse.cfg import CFG
    from genparse.semiring import Float as GenparseFloat

    cfg = CFG.from_string('''
    0.5: S -> A A
    0.5: S -> a
    0.5: A -> a
    0.5: A -> b
    ''', GenparseFloat)

    sps_cfg = build_sample_cfg(cfg, max_value=4)
    result = sps_cfg.treesum()

    # Total mass at length 2 = 0.5.
    assert abs(result.c[2].w - 0.5) < 1e-10

    np.random.seed(42)
    from collections import Counter
    counts = Counter()
    trials = 4000
    for _, (score, d) in zip(range(trials), iter(result.c[2])):
        counts[flatten_yield(d)] += 1
    expected = trials / 4
    for yld in [('a','a'), ('a','b'), ('b','a'), ('b','b')]:
        obs = counts[yld]
        z = (obs - expected) / np.sqrt(expected * 0.75)
        assert abs(z) < 4, (yld, obs, expected, z)
    assert len(counts) == 4, counts


if __name__ == '__main__':
    from arsenal import testing_framework
    testing_framework(globals())
