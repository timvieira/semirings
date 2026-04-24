from sympy import Symbol, simplify  # noqa: F401

from semirings.base import Semiring


class Symbolic(Semiring):
    """Symbolic (SymPy) interpretation: lift(w, i) returns Symbol(i).

    Not a complete semiring — provides a lift bridging into SymPy so that
    sum-product computations can be carried out symbolically and simplified.
    """
    zero = 0
    one = 1

    @staticmethod
    def lift(w, i):
        return Symbol(i)

    @staticmethod
    def symbol(i):
        return Symbol(i)
