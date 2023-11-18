from src.field import BaseFieldElement as Felt, BaseField
from dataclasses import dataclass

# ------------------------------------------------------------
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
                result = result.add(addend)
            addend = addend.double()
            scalar >>= 1

        return result

    def __neg__(self):
        if self.is_identity():
            return self
        return G1Point(self.x, -self.y)

    def __add__(self, other):
        if not isinstance(other, G1Point):
            raise TypeError("Cannot add G1Point and non-G1Point")
        if self.is_identity():
            return other
        if other.is_identity():
            return self
        if self.x == other.x and self.y != other.y:
            return G1Point(None, None)  # Result is the identity element
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
    right = pt.x**3 + 3 * pt.x + 1
    return left == right, f"{left} != {right}"


assert is_on_curve(G1)
