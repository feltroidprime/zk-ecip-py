from src.curve import G1Point, G1, Fp


class Divisor:
    points: dict[G1Point, int]  # Dictionary mapping G1Points to their multiplicities
    degree: int  # Sum of multiplicities

    def __init__(self, points: dict[G1Point, int]) -> None:
        self.points = points
        self.degree = sum(points.values())

    @staticmethod
    def empty():
        return Divisor({})

    def is_principal(self) -> bool:
        return self.degree == 0

    def __eq__(self, other: "Divisor") -> bool:
        # remove points if they have multiplicity 0:
        x_points = {p: np for p, np in self.points.items() if np != 0}
        y_points = {p: np for p, np in other.points.items() if np != 0}

        return x_points == y_points

    def __neg__(self) -> "Divisor":
        neg = Divisor({})

        neg.points = {p: -np for p, np in self.points.items()}
        neg.degree = sum(neg.points.values())

        return neg

    def __add__(self, other: "Divisor") -> "Divisor":
        res = Divisor({})

        formal_sum = self.points.copy()

        for p, np in other.points.items():
            if p in formal_sum:
                formal_sum[p] += np
            else:
                formal_sum[p] = np

        res.points = formal_sum
        res.degree = sum(formal_sum.values())

        return res

    def __sub__(self, other: "Divisor") -> "Divisor":
        return self + (-other)


if __name__ == "__main__":
    # Test
    p = G1Point.gen_random_point()
    q = G1Point.gen_random_point()
    r = G1Point.gen_random_point()

    D = Divisor({p: 1, q: 2, r: 3})
    print(f"D: {D.points} \ndegree: {D.degree}\n")

    neg = -D
    print(f"-D: {neg.points} \ndegree: {neg.degree}\n")

    zero = D + neg
    print(f"D-D: {zero.points} \ndegree: {zero.degree}\n")

    assert zero == Divisor.empty()

    diff = D - neg
    print(f"D-(-D): {diff.points} \ndegree: {diff.degree}\n")

    double = D + D
    print(f"D+D: {double.points} \ndegree: {double.degree}\n")

    assert double == diff
