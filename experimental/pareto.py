"""
Pareto-front semiring sketch (work-in-progress; was a commented-out test_pareto
stub at the top of tests/test_semirings.py).

Elements are Pareto-optimal sets of points in R^d. With vector addition as *
and set-union-then-prune as +:

    X + Y = pareto(X ∪ Y)                   (keep only non-dominated points)
    X * Y = pareto({x + y : x in X, y in Y})
    zero  = []
    one   = [origin]

The original stub only sketched the 2D case with pointwise `<` domination:

    def pareto(Y):
        return np.array([y1 for y1 in Y
                         if not any(np.all(y2 < y1) for y2 in Y if y1 is not y2)])

Unresolved before this can be promoted to a real semiring:
  - Distributivity: probably holds (this is a well-known semiring in
    multi-objective optimization) but the stub never tested it.
  - `one = [np.zeros(2)]` hardcodes dimension 2; a real impl should be
    parametric over d.
  - Equality/hashing of Pareto fronts: point sets modulo ordering, with
    floating-point coordinates. Need a canonical form (sorted? rounded?).
  - The `show(points, c='r')` helper in the stub used pylab for visualization.
    Kept out of the main module since it's dev-time only.
"""

import numpy as np

from semirings import make_semiring


def pareto(Y):
    """2D Pareto front: keep points y1 such that no other y2 dominates y1."""
    return np.array([
        y1 for y1 in Y
        if not any(np.all(y2 < y1) for y2 in Y if y1 is not y2)
    ])


def build():
    return make_semiring(
        'Pareto',
        lambda X, Y: pareto(list(X) + list(Y)),
        lambda X, Y: pareto([x + y for x in X for y in Y]),
        [],
        [np.zeros(2)],
    )


def show(points, c='r'):  # pragma: no cover — dev-time pylab scatter plot
    import pylab as pl
    if hasattr(points, 'x'):
        points = points.x
    points = np.array(points)
    P = pareto(points)
    pl.scatter(points[:, 0], points[:, 1], alpha=0.5, c='b')
    pl.scatter(P[:, 0], P[:, 1], c=c, alpha=0.5)


if __name__ == '__main__':
    Pareto = build()
    points = np.random.uniform(size=(30, 2))
    print('pareto front:', pareto(points))
