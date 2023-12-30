from src.field import BaseFieldElement as Felt, BaseField
from dataclasses import dataclass
import random

# ------------------------------------------------------------
# Generic Weierstrass curve:
# Equation : y^2 = x^3 + Ax + B
# BN254 Curve:
# Equation : y^2 = x^3 + 3
# A = 0, B = 3
A = 0
B = 3
# Curve prime field
P = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
# Curve order
N = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001
# ------------------------------------------------------------

# Encode special +inf points as non-rechable integer in the field
INF = -1

Fp = BaseField(P)

zero = Fp.zero()


@dataclass
class G1Point:
    x: Felt
    y: Felt

    def __str__(self) -> str:
        if self.is_identity():
            return "Identity"
        return f"X: {self.x.value}\nY: {self.y.value}"

    def __repr__(self) -> str:
        if self.is_identity():
            return "G1P(Identity)"
        x_hex = hex(self.x.value)
        y_hex = hex(self.y.value)
        x_formatted = x_hex[:7] + "..." + x_hex[-5:]
        y_formatted = y_hex[:7] + "..." + y_hex[-5:]

        return f"G1P(x={x_formatted}, y={y_formatted})"

    @classmethod
    def zero(cls) -> "G1Point":
        return cls(None, None)

    def is_identity(self) -> bool:
        return self.x is None and self.y is None

    def double(self) -> "G1Point":
        slope = (3 * self.x * self.x) / (2 * self.y)
        nx = slope * slope - 2 * self.x
        ny = slope * (self.x - nx) - self.y
        return G1Point(nx, ny)

    def add(self, other: "G1Point") -> "G1Point":
        slope = (other.y - self.y) / (other.x - self.x)
        nx = slope * slope - self.x - other.x
        ny = slope * (self.x - nx) - self.y
        return G1Point(nx, ny)

    def scalar_mul(self, scalar: int) -> "G1Point":
        if self.is_identity():
            return self
        if scalar == 0:
            return G1Point(None, None)
        if scalar < 0:
            return -self.scalar_mul(-scalar)
        result = G1Point(None, None)  # Identity
        addend = self
        while scalar:
            if scalar & 1:
                result += addend
            addend = addend.double()
            scalar >>= 1

        return result

    @staticmethod
    def gen_random_point() -> "G1Point":
        scalar = random.randint(1, N - 1)
        return G1.scalar_mul(scalar)

    def __eq__(self, other):
        if not isinstance(other, G1Point):
            raise TypeError("Cannot compare G1Point with non-G1Point")
        if self.is_identity() and other.is_identity():
            return True
        elif self.is_identity() or other.is_identity():
            return False
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __neg__(self):
        if self.is_identity():
            return self
        return G1Point(self.x, -self.y)

    def __add__(self, other):
        if not isinstance(other, G1Point):
            raise TypeError("Cannot add G1Point and non-G1Point")

        # Check for the identity element
        if self.is_identity():
            return other
        if other.is_identity():
            return self

        # Check for point doubling
        if self.x == other.x and self.y == other.y:
            return self.double()

        # Check for the additive inverse (result is the identity element)
        if self.x == other.x and self.y != other.y:
            return G1Point(None, None)

        return self.add(other)

    def __sub__(self, other):
        if not isinstance(other, G1Point):
            raise TypeError("Cannot subtract non-G1Point from G1Point")
        return self + (-other)

    def __mul__(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply G1Point by an integer")
        return self.scalar_mul(scalar)

    def __rmul__(self, scalar):
        return self.__mul__(scalar)


G1 = G1Point(Fp(1), Fp(2))
POINT_AT_INFINITY = G1Point(None, None)


def is_on_curve(pt: G1Point):
    left = pt.y**2
    right = pt.x**3 + A * pt.x + B
    return left == right, f"{left} != {right}"


if __name__ == "__main__":
    assert is_on_curve(G1)

    print(f"2*G1 :\n{G1.scalar_mul(2)}")
    assert G1.scalar_mul(2) == G1.double()

    assert G1.scalar_mul(N) == POINT_AT_INFINITY

    random_point = G1Point.gen_random_point()

    assert is_on_curve(random_point)
