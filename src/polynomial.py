from src.field import *


class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = [c for c in coefficients]

    def degree(self):
        if self.coefficients == []:
            return -1
        zero = self.coefficients[0].field.zero()
        if self.coefficients == [zero] * len(self.coefficients):
            return -1
        maxindex = 0
        for i in range(len(self.coefficients)):
            if self.coefficients[i] != zero:
                maxindex = i
        return maxindex

    def get_coeffs(self):
        coeffs = [x.value % x.field.p for x in self.coefficients]
        while len(coeffs) > 0 and coeffs[-1] == 0:
            coeffs.pop()
        return coeffs

    def derivative(self):
        """Compute the derivative of the polynomial."""
        if self.degree() == 0:
            # The derivative of a constant is 0
            return Polynomial([self.coefficients[0].field.zero()])

        # Compute the derivative
        derivative_coeffs = [
            i * self.coefficients[i] for i in range(1, len(self.coefficients))
        ]
        return Polynomial(derivative_coeffs)

    def __neg__(self):
        return Polynomial([-c for c in self.coefficients])

    def __add__(self, other):
        if self.degree() == -1:
            return other
        elif other.degree() == -1:
            return self
        field = self.coefficients[0].field
        coeffs = [field.zero()] * max(len(self.coefficients), len(other.coefficients))
        for i in range(len(self.coefficients)):
            coeffs[i] = coeffs[i] + self.coefficients[i]
        for i in range(len(other.coefficients)):
            coeffs[i] = coeffs[i] + other.coefficients[i]
        return Polynomial(coeffs)

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other):
        if self.coefficients == [] or other.coefficients == []:
            return Polynomial([])
        zero = self.coefficients[0].field.zero()
        buf = [zero] * (len(self.coefficients) + len(other.coefficients) - 1)
        for i in range(len(self.coefficients)):
            if self.coefficients[i].is_zero():
                continue  # optimization for sparse polynomials
            for j in range(len(other.coefficients)):
                buf[i + j] = buf[i + j] + self.coefficients[i] * other.coefficients[j]
        return Polynomial(buf)

    def __truediv__(self, other):
        quo, rem = Polynomial.divide(self, other)
        assert (
            rem.is_zero()
        ), "cannot perform polynomial division because remainder is not zero"
        return quo

    def __floordiv__(self, other):
        quo, rem = Polynomial.divide(self, other)
        return quo

    def __mod__(self, other):
        quo, rem = Polynomial.divide(self, other)
        return rem

    def __eq__(self, other):
        assert type(self) == type(
            other
        ), f"type of self {type(self)} must be equal to type of other which is {type(other)}"
        if self.degree() != other.degree():
            return False
        if self.degree() == -1:
            return True
        return all(
            self.coefficients[i] == other.coefficients[i]
            for i in range(min(len(self.coefficients), len(other.coefficients)))
        )

    def __neq__(self, other):
        return not self.__eq__(other)

    def is_zero(self):
        if self.degree() == -1:
            return True
        return False

    def __str__(self):
        return "[" + ",".join(s.__str__() for s in self.coefficients) + "]"

    def leading_coefficient(self):
        return self.coefficients[self.degree()]

    def divide(numerator, denominator):
        if denominator.degree() == -1:
            return None
        if numerator.degree() < denominator.degree():
            return (Polynomial([]), numerator)
        field = denominator.coefficients[0].field
        remainder = Polynomial([n for n in numerator.coefficients])
        quotient_coefficients = [
            field.zero() for i in range(numerator.degree() - denominator.degree() + 1)
        ]
        for i in range(numerator.degree() - denominator.degree() + 1):
            if remainder.degree() < denominator.degree():
                break
            coefficient = (
                remainder.leading_coefficient() / denominator.leading_coefficient()
            )
            shift = remainder.degree() - denominator.degree()
            subtractee = (
                Polynomial([field.zero()] * shift + [coefficient]) * denominator
            )
            quotient_coefficients[shift] = coefficient
            remainder = remainder - subtractee
        quotient = Polynomial(quotient_coefficients)
        return quotient, remainder

    def is_zero(self):
        if self.coefficients == []:
            return True
        for c in self.coefficients:
            if not c.is_zero():
                return False
        return True

    @staticmethod
    def lagrange_interpolation(domain, values):
        assert len(domain) == len(
            values
        ), "number of elements in domain does not match number of values -- cannot interpolate"
        assert len(domain) > 0, "cannot interpolate between zero points"
        field = domain[0].field
        X = Polynomial([field.zero(), field.one()])
        acc = Polynomial([])
        for i in range(len(domain)):
            prod = Polynomial([values[i]])
            for j in range(len(domain)):
                if j == i:
                    continue
                prod = (
                    prod
                    * (X - Polynomial([domain[j]]))
                    * Polynomial([(domain[i] - domain[j]).inverse()])
                )
            acc = acc + prod
        return acc

    @staticmethod
    def hermite_interpolation(points, values, derivatives):
        n = len(points)
        assert (
            len(values) == n and len(derivatives) == n
        ), "Lengths of inputs must be equal"

        field = points[0].field
        X = Polynomial([field.zero(), field.one()])  # Polynomial x
        acc = Polynomial([field.zero()])  # Accumulator polynomial

        for i in range(n):
            # Construct the Lagrange basis polynomial for the ith point
            l_i = Polynomial([field.one()])
            for j in range(n):
                if j != i:
                    l_i *= (X - Polynomial([points[j]])) * Polynomial(
                        [(points[i] - points[j]).inverse()]
                    )

            q_i = l_i * l_i  # Square the Lagrange basis polynomial
            q_i_prime = q_i.derivative()
            p_i = Polynomial([values[i]]) + (X - Polynomial([points[i]])) * (
                Polynomial([derivatives[i] - q_i_prime.evaluate(points[i]) * values[i]])
            )

            acc += q_i * p_i

        return acc

    def zerofier_domain(domain):
        field = domain[0].field
        x = Polynomial([field.zero(), field.one()])
        acc = Polynomial([field.one()])
        for d in domain:
            acc = acc * (x - Polynomial([d]))
        return acc

    def evaluate(self, point):
        xi = point.field.one()
        value = point.field.zero()
        for c in self.coefficients:
            value = value + c * xi
            xi = xi * point
        return value

    def evaluate_domain(self, domain):
        return [self.evaluate(d) for d in domain]

    def __xor__(self, exponent):
        if self.is_zero():
            return Polynomial([])
        if exponent == 0:
            return Polynomial([self.coefficients[0].field.one()])
        acc = Polynomial([self.coefficients[0].field.one()])
        for i in reversed(range(len(bin(exponent)[2:]))):
            acc = acc * acc
            if (1 << i) & exponent != 0:
                acc = acc * self
        return acc

    def xgcd(x, y):
        one = Polynomial([x.coefficients[0].field.one()])
        zero = Polynomial([x.coefficients[0].field.zero()])
        old_r, r = (x, y)
        old_s, s = (one, zero)
        old_t, t = (zero, one)

        while not r.is_zero():
            quotient = old_r // r
            old_r, r = (r, old_r - quotient * r)
            old_s, s = (s, old_s - quotient * s)
            old_t, t = (t, old_t - quotient * t)

        lcinv = old_r.coefficients[old_r.degree()].inverse()

        # a, b, g
        return (
            Polynomial([c * lcinv for c in old_s.coefficients]),
            Polynomial([c * lcinv for c in old_t.coefficients]),
            Polynomial([c * lcinv for c in old_r.coefficients]),
        )


def test_colinearity(points):
    domain = [p[0] for p in points]
    values = [p[1] for p in points]
    polynomial = Polynomial.lagrange_interpolation(domain, values)
    return polynomial.degree() == 1


if __name__ == "__main__":
    from src.curve import P
    from random import randint as rint

    field = BaseField(P)
    MAX_DEGREE = 10
    N_TESTS = 100

    for _ in range(N_TESTS):
        # Random polynomial and its derivative
        F = Polynomial(
            [
                BaseFieldElement(rint(0, P - 1), field)
                for _ in range(rint(1, MAX_DEGREE + 1))
            ]
        )
        dF = F.derivative()

        # Points at which we want to interpolate
        points = [
            BaseFieldElement(rint(0, P - 1), field) for _ in range(2 * F.degree() + 2)
        ]

        # Corresponding values and derivatives
        values = [F.evaluate(x) for x in points]
        derivatives = [dF.evaluate(x) for x in points]

        # Perform Lagrange interpolation
        lagrange_poly = Polynomial.lagrange_interpolation(points, values)
        # Perform Hermite interpolation
        hermite_poly = Polynomial.hermite_interpolation(points, values, derivatives)
        dhermite_poly = hermite_poly.derivative()

        print(f"Testing the Lagrange interpolation")
        for x in points:
            assert lagrange_poly.evaluate(x) == F.evaluate(
                x
            ), f"Mismatch at x = {x}, f(x) = {F(x)}, p(x) = {lagrange_poly.evaluate(x)}"

        assert (
            lagrange_poly == F
        ), f"Lagrange and Hermite polynomials differ, {lagrange_poly} != {F}"

        print("Lagrange interpolation successful!")

        print("Testing Hermite Interpolation")
        for i, x in enumerate(points):
            hermite_val = hermite_poly.evaluate(x)
            expected_val = F.evaluate(x)
            hermite_deriv_val = dhermite_poly.evaluate(x)
            expected_deriv_val = dF.evaluate(x)

            assert (
                hermite_val == expected_val
            ), f"Mismatch at x = {x}, f(x) = {expected_val}, p(x) = {hermite_val}"

            assert (
                hermite_deriv_val == expected_deriv_val
            ), f"Mismatch at x = {x}, f'(x) = {expected_deriv_val}, p'(x) = {hermite_deriv_val}"

        print("Hermite interpolation successful!")

        print(
            f"lagrange_poly: {lagrange_poly.get_coeffs()}, degree : {lagrange_poly.degree()}"
        )
        print(
            f"hermite_poly: {hermite_poly.get_coeffs()}, degree : {hermite_poly.degree()}"
        )
        print(f"Original F: {F.get_coeffs()}, degree : {F.degree()}")
