from src.field import BaseFieldElement as Felt, BaseField
from dataclasses import dataclass

# ------------------------------------------------------------
# Equation : y^2 = x^3 + 3
# Curve prime field
P = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
# Curve order
N = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001
# ------------------------------------------------------------

Fp = BaseField(P)

zero = Fp.zero()


@dataclass
class G1Point:
    x: Felt
    y: Felt

    def __str__(self) -> str:
        return f"X: {self.x}\nY: {self.y}"

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
        result = G1Point(None, None)  # Identity
        addend = self
        while scalar:
            if scalar & 1:
                result += addend
            addend = addend.double()
            scalar >>= 1

        return result

    def __eq__(self, other):
        if not isinstance(other, G1Point):
            raise TypeError("Cannot compare G1Point with non-G1Point")
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
    right = pt.x**3 + 3
    return left == right, f"{left} != {right}"


assert is_on_curve(G1)

print(G1.scalar_mul(2))
print(G1.scalar_mul(N - 1))
assert G1.scalar_mul(N) == POINT_AT_INFINITY
