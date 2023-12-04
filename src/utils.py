import random


def neg_3_base_le(scalar):
    """
    Decomposes a scalar into base -3 representation.
    :param scalar: The integer to be decomposed.
    :return: A list of coefficients in base -3 representation. (Least significant bit first),
    with digits [-1, 0, 1]
    """
    if scalar == 0:
        return [0]

    digits = []
    while scalar != 0:
        remainder = scalar % 3
        if (
            remainder == 2
        ):  # if the remainder is 2, we set it to -1 and add 1 to the next digit
            remainder = -1
            scalar += 1
        # For remainder 1 and 0, no change is required
        digits.append(remainder)
        scalar = -(scalar // 3)  # divide by -3 for the next digit

    return digits


if __name__ == "__main__":
    random.seed(0)
    rscalars = [random.randint(0, 100000) for _ in range(1000)]

    def test_neg3():
        for scalar in rscalars:
            decomposed_scalar = neg_3_base_le(scalar)
            eval = [digit * (-3) ** i for i, digit in enumerate(decomposed_scalar)]

            print(scalar, decomposed_scalar, eval)
            assert scalar == sum(eval)

    test_neg3()
