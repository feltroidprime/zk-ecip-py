from src.curve import G1Point, G1, Fp


class Divisor:
    points: list[(int, G1Point)]
    value: G1Point
    degree: int

    def __init__(self, points: list[(int, G1Point)]) -> None:
        self.points = points
        self.value = G1Point.zero()
        for np, point in points:
            self.value += point.scalar_mul(np)
        self.degree = sum(np for np, _ in points)

    @staticmethod
    def empty():
        return Divisor([])

    def __neg__(self) -> "Divisor":
        neg = Divisor([])

        neg.value = -self.value
        neg.degree = -self.degree
        neg.points = [(-np, p) for np, p in self.points]

        return neg

    def __add__(self, other: "Divisor") -> "Divisor":
        sum = Divisor([])

        sum.value = self.value + other.value
        sum.degree = self.degree + other.degree
        sum.points = self.points + other.points

        return sum

    def __sub__(self, other: "Divisor") -> "Divisor":
        return self + (-other)
