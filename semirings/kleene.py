import numpy as np
from semirings.base import Semiring
from arsenal import Integerizer


def MatrixSemiring(base_semiring, domain):
    """
    Factory that returns a semiring class for |D|x|D| matrices
    over the given base semiring, where D is the domain.
    """

    S = base_semiring
    sigma = Integerizer(list(domain))
    sigma.freeze()
    n = len(sigma)

    class Matrix(Semiring):

        def __init__(self, data=None):
            if isinstance(data, np.ndarray):
                self.data = data
            elif isinstance(data, dict):
                self.data = np.empty((n, n), dtype=object)
                for i in range(n):
                    for j in range(n):
                        self.data[i, j] = S.zero
                for (di, dj), v in data.items():
                    self.data[sigma(di), sigma(dj)] = v
            elif data is None:
                self.data = np.empty((n, n), dtype=object)
                for i in range(n):
                    for j in range(n):
                        self.data[i, j] = S.zero
            else:
                raise TypeError(f'Expected np.ndarray or dict, got {type(data)}')

        def __add__(self, other):
            result = np.empty((n, n), dtype=object)
            for i in range(n):
                for j in range(n):
                    result[i, j] = self.data[i, j] + other.data[i, j]
            return Matrix(result)

        def __mul__(self, other):
            result = np.empty((n, n), dtype=object)
            for i in range(n):
                for j in range(n):
                    acc = S.zero
                    for k in range(n):
                        acc = acc + self.data[i, k] * other.data[k, j]
                    result[i, j] = acc
            return Matrix(result)

        def star(self):
            return Matrix(kleene(self.data, S, reflexive=True))

        def __eq__(self, other):
            if not isinstance(other, Matrix):
                return NotImplemented
            for i in range(n):
                for j in range(n):
                    if self.data[i, j] != other.data[i, j]:
                        return False
            return True

        def __hash__(self):
            return hash(tuple(
                tuple(self.data[i, j] for j in range(n))
                for i in range(n)
            ))

        def __getitem__(self, item):
            di, dj = item
            return self.data[sigma(di), sigma(dj)]

        def __setitem__(self, item, value):
            di, dj = item
            self.data[sigma(di), sigma(dj)] = value

        def metric(self, other):
            return max(
                self.data[i, j].metric(other.data[i, j])
                for i in range(n) for j in range(n)
            )

        def __repr__(self):
            rows = []
            for i in range(n):
                row = ', '.join(repr(self.data[i, j]) for j in range(n))
                rows.append(f'[{row}]')
            return f'Matrix([{", ".join(rows)}])'

    # Zero matrix
    zero_data = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            zero_data[i, j] = S.zero
    Matrix.zero = Matrix(zero_data)

    # Identity matrix
    one_data = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            one_data[i, j] = S.one if i == j else S.zero
    Matrix.one = Matrix(one_data)

    Matrix.__name__ = f'Matrix[{S.__name__}]'

    return Matrix


def kleene(A, semiring, reflexive=True):
    """
    Compute the matrix star
    """

    # initialization
    [N,_] = A.shape; zero = semiring.zero; one = semiring.one
    new = A.copy(); old = A.copy()
    for j in range(N):
        new[:,:] = zero
        sjj = semiring.star(old[j,j])
        for i in range(N):
            for k in range(N):
                # i ➙ j ⇝ j ➙ k
                new[i,k] = old[i,k] + old[i,j] * sjj * old[j,k]
        old, new = new, old   # swap to repurpose space
    if reflexive:  # post processing fix-up: add the identity matrix
        for i in range(N): old[i,i] += one
    return old
