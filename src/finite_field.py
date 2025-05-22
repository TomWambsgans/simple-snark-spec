P = 2130706433  # Prime field of our finite field
P_BITS = 31  # Log2(P)
TWO_ADICITY = 24  # Largest n such that 2^n divides P-1
TWO_ADIC_GENERATOR = 1791270792  # A generator of the multiplicative subgroup of size 2^TWO_ADICITY.


class Fp:
    def __init__(self, value: int):
        self.value = value % P

    def __add__(self, other):
        return Fp((self.value + other.value))

    def __sub__(self, other):
        return Fp((self.value - other.value))

    def __mul__(self, other):
        return Fp((self.value * other.value))

    def __pow__(self, exponent: int):
        return Fp(pow(self.value, exponent, P))

    def __eq__(self, other):
        return self.value == other.value

    def two_addic_generator(bits: int):
        assert bits <= TWO_ADICITY
        return Fp(TWO_ADIC_GENERATOR) ** (2**(TWO_ADICITY - bits))

    def __repr__(self):
        return f"Fp({self.value})"
