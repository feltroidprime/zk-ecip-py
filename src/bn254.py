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

    def double(self) -> "G1Point":
        slope = (3 * self.x * self.x) / (2 * self.y)
        nx = slope * slope - 2 * self.x
        ny = slope * (self.x - nx) - self.y
        return G1Point(nx, ny)

    def add(self) -> "G1Point":
        slope = (self.y - self.y) / (self.x - self.x)
        nx = slope * slope - self.x - self.x
        ny = slope * (self.x - nx) - self.y
        return G1Point(nx, ny)


G1 = G1Point(Fp(1), Fp(2))


def is_on_curve(pt: G1Point):
    left = pt.y**2
    right = pt.x**3 + 3 * pt.x + 1
    return left == right, f"{left} != {right}"


assert is_on_curve(G1)
