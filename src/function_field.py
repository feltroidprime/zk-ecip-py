from dataclasses import dataclass
from random import randint as rint
from src.polynomial import Polynomial
from src.rational_function import RationalFunction
from src.field import BaseFieldElement, BaseField
from src.divisor import Divisor
from src.curve import P, Fp, A, B, G1Point, POINT_AT_INFINITY


@dataclass
class FunctionFelt:
    a: RationalFunction
    b: RationalFunction
    """
    A function field element as two rational functions a(x) and b(x)
    f(x,y) = a(x) - y*b(x) mod (y^2 - x^3 - A*x - B) where
    the curve's equation is y^2 = x^3 + A*x + B
    """

    def norm(self) -> RationalFunction:
        """
        Return the norm of this element as a rational function in x.
        N(f) = f(x,y) * f(x,-y) = N(x) = a(x)^2 - (x^3 + A*x + B) * b(x)^2
        See section 2.2.
        """
        ONE = Polynomial([Fp.one()])
        X_CUBE = RationalFunction(
            Polynomial([Fp.zero(), Fp.zero(), Fp.zero(), Fp.one()]), ONE
        )
        X = RationalFunction(Polynomial([Fp.zero(), Fp.one()]), ONE)

        return (
            self.a * self.a
            - (X_CUBE + A * X + Polynomial([BaseFieldElement(B, Fp)])) * self.b * self.b
        )

    def evaluate(self, pt: G1Point) -> BaseFieldElement:
        """
        Evaluate this function field element at the point pt(x,y).
        """
        return self.a.evaluate(pt.x) - pt.y * self.b.evaluate(pt.x)

    @staticmethod
    def gen_random(max_degree: int = 5) -> "FunctionFelt":
        """
        Generate a random function field element with numerator/denominator of maximum degree num_degree and den_degree respectively.
        """
        polys = []
        for _ in range(4):
            polys.append(
                Polynomial([Fp(rint(0, P - 1)) for _ in range(rint(1, max_degree + 1))])
            )
        return FunctionFelt(
            RationalFunction(polys[0], polys[1]), RationalFunction(polys[2], polys[3])
        )


def test_witness(f: FunctionFelt, d: Divisor) -> bool:
    """
    Test if the function field element f is correctly associated with the divisor d.
    """

    for p, np in d.points.items():
        if p == POINT_AT_INFINITY:
            continue
        if np < 0:
            raise ValueError(
                "Divisor must have points with non-negative multiplicities except for the point at infinity"
            )
        if f.evaluate(p) != Fp.zero():
            # Every point in the divisor must be a root of f
            return False

    pass


def incremental_witness(d: Divisor) -> FunctionFelt:
    """
    Compute the incremental witness for the divisor d.
    Uses incremental construction as per section 3.1.1
    """
    assert d.is_principal(), "Divisor must be principal"

    pass


def mumford_witness(d: Divisor) -> FunctionFelt:
    """
    Compute the function field element f assioated with the divisor d.
    Uses Mumford representation and Extended Euclidean Algorithm as per section 3.1.2
    """
    assert d.is_principal(), "Divisor must be principal"

    pass


if __name__ == "__main__":
    # Tests
    p = G1Point.gen_random_point()
    q = G1Point.gen_random_point()
    r = G1Point.gen_random_point()

    D = Divisor({p: 1, q: 2, r: 3, POINT_AT_INFINITY: -6})

    f1 = incremental_witness(D)
    f2 = mumford_witness(D)

    f = FunctionFelt.gen_random()

    f.evaluate(p)
    n = f.norm()
    n.evaluate(p.x)
