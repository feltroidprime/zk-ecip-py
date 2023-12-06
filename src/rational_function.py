from dataclasses import dataclass
from src.polynomial import Polynomial
from src.curve import Felt, Fp, INF


@dataclass
class RationalFunction:
    """
    A rational function is a quotient of two polynomials.
    Assumed to work in finite field Fp.
    """

    num: Polynomial
    den: Polynomial

    def evaluate(self, x: Felt) -> Felt:
        try:
            num_eval = self.num.evaluate(x)
            denom_eval = self.den.evaluate(x)
            res = num_eval / denom_eval
        except ZeroDivisionError:
            res = INF
        return res

    def __mul__(self, other) -> "RationalFunction":
        if isinstance(other, int):
            other = RationalFunction(
                Polynomial([Felt(other, Fp)]), Polynomial([Fp.one()])
            )
        if isinstance(other, Felt):
            other = RationalFunction(Polynomial([other]), Polynomial([Fp.one()]))
        if isinstance(other, Polynomial):
            other = RationalFunction(other, Polynomial([Fp.one()]))

        return RationalFunction(self.num * other.num, self.den * other.den)

    def __rmul__(self, other) -> "RationalFunction":
        return self.__mul__(other)

    def __neg__(self) -> "RationalFunction":
        return RationalFunction(-self.num, self.den)

    def __add__(self, other) -> "RationalFunction":
        if isinstance(other, int):
            other = RationalFunction(
                Polynomial([Felt(other, Fp)]), Polynomial([Fp.one()])
            )
        if isinstance(other, Felt):
            other = RationalFunction(Polynomial([other]), Polynomial([Fp.one()]))
        if isinstance(other, Polynomial):
            other = RationalFunction(other, Polynomial([Fp.one()]))

        return RationalFunction(
            self.num * other.den + self.den * other.num, self.den * other.den
        )

    def __radd__(self, other) -> "RationalFunction":
        return self.__add__(other)

    def __rsub__(self, other) -> "RationalFunction":
        return -self + other

    def __sub__(self, other) -> "RationalFunction":
        return self + (-other)
