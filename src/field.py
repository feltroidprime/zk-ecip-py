class BaseFieldElement:
    def __init__(self, value, field):
        self.value = value
        self.field = field

    def __str__(self) -> str:
        return str(self.value)

    def __repr__(self) -> str:
        return str(self.value)

    def __add__(self, right):
        if isinstance(right, int):
            right = BaseFieldElement(right, self.field)
        return self.field.add(self, right)

    def __mul__(self, right):
        if isinstance(right, int):
            right = BaseFieldElement(right, self.field)
        return self.field.multiply(self, right)

    def __sub__(self, right):
        if isinstance(right, int):
            right = BaseFieldElement(right, self.field)
        return self.field.subtract(self, right)

    def __truediv__(self, right):
        if isinstance(right, int):
            right = BaseFieldElement(right, self.field)
            assert not right.is_zero(), "divide by zero"
        return self.field.divide(self, right)

    def __radd__(self, left):
        if isinstance(left, int):
            left = BaseFieldElement(left, self.field)
        return self.field.add(left, self)

    def __rmul__(self, left):
        if isinstance(left, int):
            left = BaseFieldElement(left, self.field)
        return self.field.multiply(left, self)

    def __rsub__(self, left):
        if isinstance(left, int):
            left = BaseFieldElement(left, self.field)
        return self.field.subtract(left, self)

    def __rtruediv__(self, left):
        if isinstance(left, int):
            left = BaseFieldElement(left, self.field)
            assert not self.is_zero(), "divide by zero"
        return self.field.divide(left, self)

    def __neg__(self):
        return self.field.negate(self)

    def inverse(self):
        return self.field.inverse(self)

    # modular exponentiation -- be sure to encapsulate in parentheses!
    def __xor__(self, exponent):
        acc = BaseFieldElement(1, self.field)
        val = BaseFieldElement(self.value, self.field)
        for i in reversed(range(len(bin(exponent)[2:]))):
            acc = acc * acc
            if (1 << i) & exponent != 0:
                acc = acc * val
        return acc

    def __pow__(self, exponent):
        # If the exponent is another BaseFieldElement, use its value
        if isinstance(exponent, BaseFieldElement):
            exponent = exponent.value

        # Ensure the exponent is an integer
        if not isinstance(exponent, int):
            raise ValueError("Exponent must be an integer or BaseFieldElement")

        # Use modular exponentiation
        acc = BaseFieldElement(1, self.field)
        val = BaseFieldElement(self.value, self.field)
        for i in reversed(range(len(bin(exponent)[2:]))):
            acc = acc * acc
            if (1 << i) & exponent != 0:
                acc = acc * val
        return acc

    def __eq__(self, other):
        return self.value == other.value

    def __neq__(self, other):
        return self.value != other.value

    def __str__(self):
        return str(self.value)

    def __bytes__(self):
        return bytes(str(self).encode())

    def is_zero(self):
        if self.value == 0:
            return True
        else:
            return False

    def has_order_po2(self, order):
        assert order & (order - 1) == 0
        if self == self.field.one() and order == 1:
            return True
        return (self ^ order) == self.field.one() and not (
            self ^ (order // 2)
        ) == self.field.one()

    def __hash__(self):
        return self.value


class BaseField:
    def __init__(self, p):
        self.p = p

    def lift(self, bfe):
        return bfe

    def zero(self):
        return BaseFieldElement(0, self)

    def one(self):
        return BaseFieldElement(1, self)

    def multiply(self, left, right):
        return BaseFieldElement((left.value * right.value) % self.p, self)

    def add(self, left, right):
        return BaseFieldElement((left.value + right.value) % self.p, self)

    def subtract(self, left, right):
        return BaseFieldElement((self.p + left.value - right.value) % self.p, self)

    def negate(self, operand):
        return BaseFieldElement((self.p - operand.value) % self.p, self)

    def inverse(self, operand):
        return BaseFieldElement(pow(operand.value, -1, self.p), self)

    def divide(self, left, right):
        if right.is_zero():
            raise ZeroDivisionError("Cannot divide by zero")
        return left * right.inverse()

    def __call__(self, integer):
        return BaseFieldElement(integer % self.p, self)
