from dataclasses import dataclass
import random
from random import randint as rint
from src.polynomial import Polynomial
from src.rational_function import RationalFunction
from src.field import BaseFieldElement, BaseField
from src.divisor import Divisor
from src.curve import P, Fp, A, B, G1Point, POINT_AT_INFINITY


@dataclass
class FunctionFelt:
    a: Polynomial
    b: Polynomial
    """
    A function field element as two rational functions a(x) and b(x)
    f(x,y) = a(x) - y*b(x) mod (y^2 - x^3 - A*x - B) where
    the curve's equation is y^2 = x^3 + A*x + B
    """

    def __repr__(self) -> str:
        return f"FunctionFelt(a deg({self.a.degree()}) -y b deg({self.b.degree()}))"

    def norm(self) -> Polynomial:
        """
        Return the norm of this element as a rational function in x.
        N(f) = f(x,y) * f(x,-y) = N(x) = a(x)^2 - (x^3 + A*x + B) * b(x)^2
        See section 2.2.
        """
        ONE = Polynomial([Fp.one()])
        X_CUBE = Polynomial([Fp.zero(), Fp.zero(), Fp.zero(), Fp.one()])
        X = Polynomial([Fp.zero(), Fp.one()])

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
        Generate a random function field element with numerator/denominator of maximum degree max_degree.
        """
        polys = []
        for _ in range(2):
            polys.append(
                Polynomial([Fp(rint(0, P - 1)) for _ in range(rint(1, max_degree + 1))])
            )
        return FunctionFelt(a=polys[0], b=polys[1])

    def __mul__(self, other: "FunctionFelt") -> "FunctionFelt":
        """
        Multiply two function field elements. Susbtitutes y^2 for x^3 + A*x + B.

        (a(x) - yb(x)) * (a'(x) - yb'(x)) = a(x)a'(x) - y(a(x)b'(x) + a'(x)b(x)) + y^2b(x)b'(x)
        = (a(x)a'(x) + y^2*b(x)b'(x)) - y(a(x)b'(x) + a'(x)b(x))
        = (a(x)a'(x) + (x^3 + A*x + B)b(x)b'(x)) - y(a(x)b'(x) + a'(x)b(x))

        res_a = a(x)a'(x) + (x^3 + A*x + B)b(x)b'(x)
        res_b = a(x)b'(x) + a'(x)b(x)

        res = res_a - y*res_b
        """
        if not isinstance(other, FunctionFelt):
            raise TypeError("Can only multiply FunctionFelt by another FunctionFelt")
        # y2 = X^3 + A*X + B
        y2 = Polynomial(
            [BaseFieldElement(B, Fp), BaseFieldElement(A, Fp), Fp.zero(), Fp.one()]
        )
        res_b = self.a * other.b + self.b * other.a
        res_a = self.a * other.a + y2 * self.b * other.b
        return FunctionFelt(a=res_a, b=res_b)


def test_witness(f: FunctionFelt, d: Divisor) -> bool:
    """
    Test if the function field element f is correctly associated with the divisor d.
    """
    # TODO : check multpiplicity of roots

    for p, np in d.points.items():
        if p == POINT_AT_INFINITY:
            continue
        if np < 0:
            raise ValueError(
                "Divisor must have points with non-negative multiplicities except for the point at infinity"
            )
        if f.evaluate(p) == Fp.zero():
            print(f"f({p}) = 0")
        if f.evaluate(p) != Fp.zero():
            # Every point in the divisor must be a root of f
            print(f"f({p}) != 0")
            return False

    return True


def incremental_witness(d: Divisor) -> FunctionFelt:
    """
    Compute the incremental witness for the divisor d.
    Uses incremental construction as per section 3.1.1
    """
    assert d.is_principal(), "Divisor must be principal"

    pass


def mumford_witness(d: Divisor) -> FunctionFelt:
    """
    Compute the function field element f assiociated with the divisor d.
    Uses Mumford representation and Extended Euclidean Algorithm as per section 3.1.2
    """
    assert d.is_principal(), "Divisor must be principal"
    # Init res to neutral element for multiplication of function field elements
    res = FunctionFelt(a=Polynomial([Fp.one()]), b=Polynomial([Fp.zero()]))
    print(f"Divisor: {d}")
    X = Polynomial([Fp.zero(), Fp.one()])
    for point, m in d.points.items():
        print(point)
        print(m)
        if point == POINT_AT_INFINITY:
            continue
        if m < 0:
            raise ValueError(
                "Divisor must have points with non-negative multiplicities except for the point at infinity"
            )
        if m == 0:
            continue
        if m == 1:
            # u(x) = x - x_p => u(x_p) = 0
            u = X - Polynomial([point.x])
            # v(x_p) = y_p. Constant polynomial
            v = Polynomial([point.y])
            (_, b, a) = Polynomial.xgcd(u, v)
            # f(x,y) = a(x) - y*b(x)
            res = res * FunctionFelt(a=a, b=b)
            continue
        if m > 1:
            # Appendix 7.1.
            # Hensel Lifting, Appendix 7.1.
            v = []
            r = []
            q = []

            i = 0
            v.append(Polynomial([point.y]))
            print(f"v_{i}: {v[i].get_coeffs()}")

            r.append(
                (
                    Polynomial([point.y]) ** 2
                    - (X**3 + A * X + Polynomial([BaseFieldElement(B, Fp)]))
                )
                / (X - Polynomial([point.x]))
            )

            q.append(((-2 * r[i]) / v[i]) % ((X - Polynomial([point.x])) ** (2**i)))
            print(f"q_{i}: {q[i].get_coeffs()}")

            while True:
                check_poly = v[i] ** 2 - (
                    X**3 + A * X + Polynomial([BaseFieldElement(B, Fp)])
                )
                mod_check = (X - Polynomial([point.x])) ** (2 ** (i))
                if (check_poly % mod_check).is_zero() == False:
                    print(
                        f"Stage {i} not working, {check_poly.get_coeffs()} % {mod_check.get_coeffs()} != 0"
                    )
                    raise ValueError(f"Error in Hensel Lifting, stage {i}")
                else:
                    print(
                        f"Stage {i} working, {check_poly.get_coeffs()} % {mod_check.get_coeffs()} == 0"
                    )

                i += 1

                print(f"i: {i}")
                X_MIN_XP_POW_2I = (X - Polynomial([point.x])) ** (2 ** (i - 1))
                v.append(v[i - 1] + X_MIN_XP_POW_2I * q[i - 1])

                if 2 ** (i - 1) < m <= 2**i:
                    # Last iteration, v(x)  = v[i] % (X - x_p) ** m
                    print(f"Last iteration {i}")
                    break

                # r.append(
                #     ((r[i - 1] + 2 * v[i - 1] * q[i - 1]) / X_MIN_XP_POW_2I)
                #     + q[i - 1] ** 2
                # )

                # q_i = -2 * r_i / v_i
                # q_i = q_i.to_poly()
                # q_i = q_i % X_MIN_XP_POW_2I

            final_v = v[i] % (X - Polynomial([point.x])) ** m
            print(f"deg(final_v) = {final_v.degree()}")
            print(f"final_v = {final_v.get_coeffs()}")
            print(f"v(xp) = {final_v.evaluate(point.x)}")
            print(f"v'(xp) = {final_v.derivative().evaluate(point.x)}")
            raise NotImplementedError("Hensel Lifting not implemented")
    return res


if __name__ == "__main__":
    # Tests
    random.seed(0)
    p = G1Point.gen_random_point()
    q = G1Point.gen_random_point()
    r = G1Point.gen_random_point()
    s = G1Point.gen_random_point()

    D = Divisor({p: 1, q: 3, r: 2, POINT_AT_INFINITY: -6})

    # f1 = incremental_witness(D)
    # f2 = mumford_witness(D)

    f = FunctionFelt.gen_random()

    f.evaluate(p)
    n = f.norm()
    n.evaluate(p.x)

    D_single_multiplicities = Divisor({p: 1, q: 1, r: 1, s: 1, POINT_AT_INFINITY: -4})
    f = mumford_witness(D_single_multiplicities)
    print(f"Function field element: {f}")
    assert test_witness(f, D_single_multiplicities) == True, f"Wrong Mumford witness"

    print(f"Test multiplicity 2")
    f = mumford_witness(D)
