def derivation(x):
    """Post-process a Viterbi/Point derivation (via `x.d` backpointer chain)
    into a nicely formatted NLTK Tree.

    Expects a semiring value carrying a `.d` field in the hypergraphs
    convention: either a leaf label, or a (body, label) pair where `body`
    is itself a derivation-carrying value or a tuple thereof.
    """
    from nltk import ImmutableTree
    if not isinstance(x.d, (list, tuple)):
        return x.d
    assert len(x.d) == 2
    label = x.d[1].d
    body = x.d[0].d
    if isinstance(body, (list, tuple)):
        return ImmutableTree(label, list(map(derivation, body)))
    else:
        return ImmutableTree(label, [body])


def post_process(x):
    """Convert a nested (body, label) tuple into an NLTK Tree."""
    from nltk import ImmutableTree
    if not isinstance(x, tuple):
        return x
    [body, label] = x
    if isinstance(body, tuple):
        return ImmutableTree(label, list(map(post_process, body)))
    else:
        return ImmutableTree(label, [body])
