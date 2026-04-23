"""
"Funky" semiring sketch (work-in-progress; was test_funky in tests/test_semirings.py).

The semiring requirements can be relaxed in many ways that still admit general
algorithms. For example, we can relax the requirement that * is a monoid. We
can allow it to be a multi-arity operator with no structure other than
closure:

    \\forall x_1, ..., x_K \\in \\mathcal{X}:  f(x_1, ..., x_K) \\in \\mathcal{X}

But we still want a kind of distributive property:

    \\forall k:  \\sum_{x_k} f(x_1, ..., x_k, ..., x_K)
               = f(x_1, ..., (\\sum_{x_k} x_k), ..., x_K)

Example from Gilea (2020), efficient outside computation:
https://aclanthology.org/2020.cl-4.2/

`F` below is a gestural placeholder for the K-arity operator — it's not wired
up, just a reminder of what the "multi-arity" part would look like.

The `S = make_semiring(...)` call currently builds an ordinary binary-`*`
semiring ("Funky" because the * is min(x, x+exp(y)), which is not really the
multi-arity operator this note is about). It was checking axioms with
`assoc=False, star=False, hash_=False` — i.e., disabling most of the
interesting axioms, so it wasn't really a test of anything specific. Keeping
it here as a reference point, not a test.
"""

import numpy as np
from semirings import make_semiring, check_axioms_samples


def F(x1, x2):
    """Placeholder for a K-arity operator; not implemented."""
    return


def build():
    S = make_semiring(
        'Funky',
        min,
        lambda x, y: x + np.exp(y),
        np.inf,
        1,   # ???
    )
    return S


if __name__ == '__main__':
    S = build()
    members = list(map(S, np.linspace(-5, 5, 11)))
    print(members)
    check_axioms_samples(S, members, hash_=False, assoc=False, star=False)
