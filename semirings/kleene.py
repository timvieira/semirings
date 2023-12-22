# TODO: Implement the semiring of square matrices


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
