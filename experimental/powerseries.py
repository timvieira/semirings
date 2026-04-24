"""Power-series semirings — truncated and (stub) lazy.

This module provides a *general* polynomial / power-series semiring, suitable
for any dynamic program whose weights should carry a full distribution (or
generating function) rather than a single scalar. The canonical use case is
computing the PMF of an additive non-negative integer feature over a PCFG:
if the per-rule weight is ``w·x^k`` (``k`` = the feature's rule-level
contribution), then ``treesum()`` returns the generating function of the
feature's distribution over all derivations.

Mathematically, ``TruncatedPowerSeries(N)`` is the commutative semiring
``R[x]/(x^{N+1})`` — polynomials with coefficients in an underlying scalar
ring ``R`` (default ``float``), truncated at degree ``N``.

  addition:         (a + b)[n] = a[n] + b[n]             (pointwise)
  multiplication:   (a · b)[n] = sum_{i+j=n, i,j <= N} a[i]·b[j]   (convolution)
  zero:             the zero polynomial
  one:              1 + 0·x + 0·x^2 + ... = [1, 0, 0, ..., 0]
  metric:           sup-norm on coefficients  max_n |a[n] - b[n]|

All semiring laws (associativity, commutativity, distributivity, identity,
annihilation) hold for the truncated ring in the usual sense.

----------------------------------------------------------------------
Numerical precision
----------------------------------------------------------------------

PMF / GF coefficients can span wide dynamic ranges — e.g. the PCFG length
GF has coefficients that decay near the dominant singularity like
``rho^n · n^{-3/2}``, which crosses float64 underflow at modest ``N``. We
therefore expose three multiplication backends via the ``mul`` factory
argument, each with a documented precision story:

  mul='direct'      — np.convolve in float64. O(N^2). BLAS precision:
                      each output tap is a plain dot product, so error is
                      O(N) * eps * max_pointwise_product. Adequate for
                      non-negative, moderately-ranged PMFs.
  mul='longdouble'  — np.convolve in numpy's ``longdouble`` (80-bit on
                      x86, wider on some platforms), result cast back to
                      float64. O(N^2) with ~3-4 extra decimal digits in
                      the accumulator; catches most direct-path precision
                      issues with almost no code overhead.
  mul='fft'         — np.fft.rfft/irfft in float64, O(N log N). Roundoff
                      scales as O(log N * eps) in the *input norm* when
                      coefficients are of similar magnitude; when
                      coefficients span many orders of magnitude the FFT
                      mixes large and small entries and relative error
                      on the small entries can be much larger than for
                      'direct'. Prefer 'longdouble' or 'direct' when you
                      need accurate small coefficients.

Sign behaviour: for non-negative coefficients all three paths are free of
catastrophic cancellation. For sign-mixed inputs (signed measures,
derivatives of GFs, sub-critical subtraction), cancellation can amplify
error on any path; 'longdouble' is safest.

Range-scaling / tilting (the ``x -> alpha·x`` trick used in the PCFG
notebook to move coefficients away from the radius of convergence) is a
*caller-side* transformation, not the semiring's job. The semiring does
not assume anything about the provenance of its inputs.

----------------------------------------------------------------------
Example: length PMF of a PCFG
----------------------------------------------------------------------

    PS = TruncatedPowerSeries(40)

    # convert each rule 'w: A -> body' into PS element w * x^(#terminals in body)
    new_cfg = cfg.__class__(R=PS, S=cfg.S, V=cfg.V)
    for r in cfg.rules:
        k = sum(cfg.is_terminal(y) for y in r.body)
        new_cfg.add(PS.monomial(r.w, k), r.head, *r.body)

    pmf = new_cfg.treesum().c   # numpy array: pmf[n] = P(|yield(S)| = n)
"""
import numpy as np

from semirings import base


# ---------------------------------------------------------------------------
# factory cache: calling TruncatedPowerSeries(N) multiple times with the same
# arguments should return the *same* class, so elements compare / hash
# consistently across call sites.
# ---------------------------------------------------------------------------
_factory_cache = {}


def TruncatedPowerSeries(max_degree, *, scalar=float, mul='direct'):
    """Return a semiring class for ``R[x]/(x^{max_degree+1})``.

    Parameters
    ----------
    max_degree : int
        Truncation degree N (so elements hold coefficients for x^0 ... x^N).
    scalar : type, default ``float``
        Coefficient type. ``float`` is the only fully-tested option at
        present; the argument is carried through for forward-compatibility
        with custom scalar rings (``LogVal``, ``np.float128``, etc.).
    mul : {'direct', 'longdouble', 'fft'}
        Multiplication backend. See module docstring for precision notes.

    Returns
    -------
    A subclass of ``semirings.base.Semiring``. Each instance wraps a
    numpy array ``self.c`` of length ``max_degree + 1``.
    """
    if scalar is not float:
        # Until we ship a tested scalar-ring path, only accept float.
        raise NotImplementedError(
            f'TruncatedPowerSeries currently only supports scalar=float; '
            f'got {scalar!r}. (The factory arg is reserved for future use.)'
        )
    if mul not in ('direct', 'longdouble', 'fft'):
        raise ValueError(
            f"mul must be 'direct', 'longdouble', or 'fft'; got {mul!r}"
        )

    key = (max_degree, scalar, mul)
    cached = _factory_cache.get(key)
    if cached is not None:
        return cached

    N = int(max_degree)
    if N < 0:
        raise ValueError(f'max_degree must be >= 0, got {max_degree}')

    class _TPS(base.Semiring):
        # class-level metadata (useful for introspection / error messages)
        max_degree = N
        mul_mode = mul
        scalar_type = scalar

        def __init__(self, c):
            arr = np.asarray(c, dtype=float)
            if arr.ndim != 1 or arr.shape[0] != N + 1:
                raise ValueError(
                    f'coeffs must be a 1-D array of length {N+1}, '
                    f'got shape {arr.shape}'
                )
            # defensive copy: callers often mutate the source array
            self.c = arr.copy() if arr is c else arr

        # ------------------------ constructors ------------------------
        @classmethod
        def constant(cls, w):
            "Return w · x^0 = [w, 0, 0, ..., 0]."
            c = np.zeros(N + 1)
            c[0] = float(w)
            return cls(c)

        @classmethod
        def monomial(cls, w, k):
            """Return w · x^k. Returns ``zero`` when k > N (silent truncation).

            Silent truncation is the right default for PMF construction:
            callers typically pick ``max_degree`` as the maximum feature
            value they care about, and any rule whose feature exceeds the
            cap should contribute zero (its derivations are outside the
            range being modelled).
            """
            c = np.zeros(N + 1)
            if 0 <= k <= N:
                c[k] = float(w)
            return cls(c)

        @classmethod
        def from_pmf(cls, pmf):
            """Build a TPS element from a sequence or dict of coefficients.

            - sequence: ``pmf[i]`` becomes ``x^i``'s coefficient for i<=N.
            - dict:     ``pmf[k]`` becomes ``x^k``'s coefficient for k<=N.
            Entries beyond degree N are silently dropped.
            """
            c = np.zeros(N + 1)
            if hasattr(pmf, 'items'):
                for k, w in pmf.items():
                    if 0 <= k <= N:
                        c[k] += float(w)
            else:
                arr = np.asarray(pmf, dtype=float).ravel()
                m = min(arr.shape[0], N + 1)
                c[:m] = arr[:m]
            return cls(c)

        @classmethod
        def from_string(cls, s):
            """Parse a short string form.

            Accepted formats:
              'w'           -> constant(w)
              'w@k'         -> monomial(w, k)
              'w*x^k'       -> monomial(w, k)
              '[w0,w1,...]' -> from_pmf([w0, w1, ...])
            """
            s = s.strip()
            if s.startswith('['):
                import ast
                return cls.from_pmf(list(ast.literal_eval(s)))
            if '*x^' in s:
                wstr, kstr = s.split('*x^', 1)
                return cls.monomial(float(wstr), int(kstr))
            if '@' in s:
                wstr, kstr = s.split('@', 1)
                return cls.monomial(float(wstr), int(kstr))
            return cls.constant(float(s))

        # -------------------------- arithmetic ------------------------
        def __add__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            return _TPS(self.c + other.c)

        def __mul__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            if mul == 'direct':
                # np.convolve delegates to BLAS for the O(N^2) kernel; output
                # is length 2N+1, we truncate to the first N+1 coefficients.
                prod = np.convolve(self.c, other.c)
                return _TPS(prod[: N + 1].copy())
            elif mul == 'longdouble':
                a = self.c.astype(np.longdouble)
                b = other.c.astype(np.longdouble)
                prod = np.convolve(a, b)
                return _TPS(prod[: N + 1].astype(float))
            else:  # 'fft'
                size = 1
                while size < 2 * (N + 1):
                    size *= 2
                fa = np.fft.rfft(self.c, size)
                fb = np.fft.rfft(other.c, size)
                prod = np.fft.irfft(fa * fb, size)
                return _TPS(prod[: N + 1].copy())

        # ------------------------- misc -------------------------------
        def metric(self, other):
            """Sup-norm on coefficients (a true metric)."""
            if not isinstance(other, _TPS):
                return float('inf')
            return float(np.max(np.abs(self.c - other.c)))

        def __eq__(self, other):
            if not isinstance(other, _TPS):
                return NotImplemented
            return np.array_equal(self.c, other.c)

        def __hash__(self):
            return hash(self.c.tobytes())

        def __repr__(self):
            # compact coefficient display; full precision elsewhere
            coef_str = ', '.join(f'{x:g}' for x in self.c)
            return f'TPS[N={N},mul={mul}]([{coef_str}])'

        def __round__(self, precision):
            # Semiring.__round__ default returns self; override so chart
            # equality checks work sensibly after rounding.
            return _TPS(np.round(self.c, precision))

    _TPS.zero = _TPS(np.zeros(N + 1))
    one_coeffs = np.zeros(N + 1)
    one_coeffs[0] = 1.0
    _TPS.one = _TPS(one_coeffs)
    _TPS.__name__ = f'TruncatedPowerSeries_N{N}_{mul}'
    _TPS.__qualname__ = _TPS.__name__

    _factory_cache[key] = _TPS
    return _TPS


# ---------------------------------------------------------------------------
# LazyPowerSeries — TODO stub.
# ---------------------------------------------------------------------------
class LazyPowerSeries(base.Semiring):
    """Formal power series with lazy, on-demand coefficient generation.

    **TODO — not implemented yet.**

    Design sketch (for the eventual implementation):

    Representation
    --------------
    An element wraps a callable ``coeff: int -> scalar`` whose values are
    memoized on first access. No truncation: coefficients are generated
    only as they are asked for.

    Operations
    ----------
    addition:        c_new(n) = c_a(n) + c_b(n)
    multiplication:  c_new(n) = sum_{k=0..n} c_a(k) * c_b(n-k)   (Cauchy)
    zero / one:      constant-zero / indicator-at-0 coefficient functions.
    star:            1 / (1 - a)  when a(0) == 0 (implemented via the
                     Cauchy recurrence, which is well-founded under that
                     condition).
    truncate(N):     materialize the first N+1 coefficients as a
                     TruncatedPowerSeries(N) element.
    promote_to_N:    the reverse.
    __iter__:        yield c(0), c(1), c(2), ...

    Fixed points
    ------------
    ``rec(func)`` returns a series ``y`` satisfying ``y = func(y)``, for
    ``func`` that produces ``y[n]`` from ``y[0..n-1]`` (well-founded
    recursion). Implemented via memoized per-coefficient recursion — no
    explicit iteration loop needed.

    Acceleration
    ------------
    The naive Cauchy product is O(n) per coefficient, O(N^2) for the
    first N. van der Hoeven's relaxed / online FFT multiplication
    (ISSAC 2003; https://www.texmacs.org/joris/issac03/issac03.pdf) gives
    O(N log^2 N) — planned as a drop-in replacement for the multiplication
    backend.

    Why stub now?
    -------------
    The truncated variant covers every application needed by the current
    CFG tutorials. Lazy series are the right tool for open-ended GF work
    (e.g. asymptotic coefficient extraction, composition / substitution
    of series, Newton iteration on systems of series). Parking the
    scaffolding here keeps the shape of the eventual API visible while
    avoiding premature implementation.
    """

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            'LazyPowerSeries is a stub. See module docstring for the planned API. '
            'Use TruncatedPowerSeries for now.'
        )


# ---------------------------------------------------------------------------
# TruncatedPowerSeries — general polynomial / PMF semiring.
# ---------------------------------------------------------------------------
#from semirings.powerseries import TruncatedPowerSeries, LazyPowerSeries
from semirings.base import (
    check_axioms, check_metric_axioms
)
from arsenal import assert_throws


def _tps_random_samples(PS, rng, signed=True):
    """A modest sample set that exercises both structural and random elements."""
    N = PS.max_degree
    samples = [PS.zero, PS.one]
    # Random dense vectors (signed or non-negative)
    for _ in range(3):
        arr = rng.standard_normal(N + 1) if signed else rng.random(N + 1)
        samples.append(PS(arr))
    # Sparse structural elements
    samples.append(PS.constant(2.5))
    samples.append(PS.monomial(3.0, min(2, N)))
    samples.append(PS.from_pmf({0: 0.5, 1: 0.5}))
    return samples


def test_truncated_power_series_axioms():
    rng = np.random.default_rng(0)
    # Exercise several sizes and both backend variants.
    for N in [0, 1, 4, 16]:
        for mul in ('direct', 'longdouble', 'fft'):
            PS = TruncatedPowerSeries(N, mul=mul)
            samples = _tps_random_samples(PS, rng)
            # Commutative semiring: no star required. Disable the exact
            # multiplicity axioms for the 'fft' path, which introduces
            # O(log N) * eps roundoff on repeated self-addition; they are
            # checked exactly for 'direct' / 'longdouble'.
            for A in samples:
                for B in samples:
                    for C in samples:
                        check_axioms(
                            PS, A, B, C,
                            star=False,
                            check_multiplicity=(mul != 'fft'),
                        )


def test_truncated_power_series_convolution_reference():
    """Multiplication is polynomial convolution truncated at degree N."""
    rng = np.random.default_rng(1)
    N = 12
    for mul, atol in [('direct', 1e-14), ('longdouble', 1e-14), ('fft', 1e-10)]:
        PS = TruncatedPowerSeries(N, mul=mul)
        # monomial * monomial: w·x^k * u·x^j = (w u)·x^{k+j}   if k+j <= N
        for _ in range(10):
            w = float(rng.standard_normal())
            u = float(rng.standard_normal())
            k = int(rng.integers(0, N + 1))
            j = int(rng.integers(0, N + 1 - k))
            lhs = PS.monomial(w, k) * PS.monomial(u, j)
            rhs = PS.monomial(w * u, k + j)
            assert np.allclose(lhs.c, rhs.c, atol=atol), (mul, k, j, lhs, rhs)
        # Random vectors cross-checked against np.convolve.
        for _ in range(10):
            a = rng.standard_normal(N + 1)
            b = rng.standard_normal(N + 1)
            got = (PS(a) * PS(b)).c
            want = np.convolve(a, b)[: N + 1]
            assert np.allclose(got, want, atol=atol), (mul, np.max(np.abs(got - want)))


def test_truncated_power_series_metric_axioms():
    rng = np.random.default_rng(2)
    PS = TruncatedPowerSeries(6)
    samples = _tps_random_samples(PS, rng)
    check_metric_axioms(PS, samples)


def test_truncated_power_series_precision():
    """Compare backends on coefficients that span many orders of magnitude.

    Pins the precision claims in the module docstring:
    - 'direct' and 'longdouble' are accurate in sup-norm on the output.
    - 'fft' is accurate in sup-norm but can lose relative precision on
      coefficients that are orders of magnitude smaller than the largest.
    """
    rng = np.random.default_rng(3)
    N = 64
    # Reference: longdouble convolution, then down-cast.
    a = rng.random(N + 1)           # non-negative, O(1) coefficients
    b = np.logspace(0, -10, N + 1)  # coefficients spanning 10 orders of magnitude
    a_ld = a.astype(np.longdouble)
    b_ld = b.astype(np.longdouble)
    ref = np.convolve(a_ld, b_ld)[: N + 1].astype(float)

    for mul, abs_tol in [
        ('direct',     1e-12),
        ('longdouble', 1e-14),
        ('fft',        1e-10),
    ]:
        PS = TruncatedPowerSeries(N, mul=mul)
        got = (PS(a) * PS(b)).c
        err = float(np.max(np.abs(got - ref)))
        assert err < abs_tol, (mul, err)


def test_truncated_power_series_from_string_and_pmf():
    PS = TruncatedPowerSeries(5)
    assert PS.from_string('3.5') == PS.constant(3.5)
    assert PS.from_string('2.5*x^3') == PS.monomial(2.5, 3)
    assert PS.from_string('1.5@2') == PS.monomial(1.5, 2)
    assert PS.from_string('[1,2,3]') == PS.from_pmf([1, 2, 3])
    # Silent truncation beyond N.
    assert PS.monomial(1.0, 99) == PS.zero
    assert PS.from_pmf([1, 2, 3, 4, 5, 6, 7, 8, 9]) == PS.from_pmf([1, 2, 3, 4, 5, 6])


def test_truncated_power_series_factory_cached():
    "Repeated calls with the same args return the same class (for isinstance / equality)."
    A = TruncatedPowerSeries(7)
    B = TruncatedPowerSeries(7)
    assert A is B
    C = TruncatedPowerSeries(7, mul='fft')
    assert C is not A


def test_truncated_power_series_cfg_roundtrip():
    "Computing the length PMF of a simple PCFG via treesum matches closed form."
    try:
        # The genparse package provides the CFG class; skip if unavailable.
        import sys
        sys.path.insert(0, '/home/timv/projects/genparse')
        from genparse.cfg import CFG
        from genparse.semiring import Float as GenparseFloat
    except ImportError:  # pragma: no cover
        import pytest
        pytest.skip('genparse not available')

    grammar = '''
    0.5: S -> A B
    0.5: S -> a
    0.2: A -> A A
    0.8: A -> b
    0.4: B -> A c
    0.6: B -> c
    '''
    cfg = CFG.from_string(grammar, GenparseFloat)

    N = 40
    PS = TruncatedPowerSeries(N)
    new = cfg.__class__(R=PS, S=cfg.S, V=cfg.V)
    for r in cfg.rules:
        k = sum(cfg.is_terminal(y) for y in r.body)
        new.add(PS.monomial(r.w, k), r.head, *r.body)
    pmf = new.treesum().c

    # closed form: mu_S = 29/15
    mean_from_pmf = float(np.sum(np.arange(N + 1) * pmf))
    assert abs(mean_from_pmf - 29.0 / 15.0) < 1e-6, mean_from_pmf
    # PMF should normalize to 1 for a proper PCFG (up to truncation).
    assert abs(pmf.sum() - 1.0) < 1e-6


def test_lazy_power_series_is_stub():
    "LazyPowerSeries is intentionally a stub."
    with assert_throws(NotImplementedError):
        LazyPowerSeries()
