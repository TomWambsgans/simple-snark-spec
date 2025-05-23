from typing import List, Union

P = 2130706433  # Prime field of our finite field
P_BITS = 31  # Log2(P)
TWO_ADICITY = 24  # Largest n such that 2^n divides P-1
TWO_ADIC_GENERATOR = 1791270792  # A generator of the multiplicative subgroup of size 2^TWO_ADICITY (in Fp ).

DEG = 4  # Extension degree
W = 3  # such that (x^DEG - W) is irreducible over Fp


# Base field
class F:
    def __init__(self, value: int):
        self.value = value % P

    def __add__(self, other: "F") -> "F": return F(self.value + other.value)
    def __sub__(self, other: "F") -> "F": return F(self.value - other.value)
    def __mul__(self, other: "F") -> "F": return F(self.value * other.value)
    def __pow__(self, exponent: int) -> "F": return F(pow(self.value, exponent, P))
    def __eq__(self, other: object) -> bool: return self.value == other.value
    def __repr__(self) -> str: return f"F({self.value})"

    @staticmethod
    def zero() -> "F":
        return F(0)

    @staticmethod
    def two_addic_generator(bits: int) -> "F":
        assert bits <= TWO_ADICITY
        return F(TWO_ADIC_GENERATOR) ** (2**(TWO_ADICITY - bits))


# Extension field Fq (q = P^DEG)
class EF:
    def __init__(self, value: List[F]):
        assert len(value) == DEG
        self.value = value

    def from_base(value: F) -> "EF":
        return EF([value] + [F(0) for _ in range(DEG - 1)])

    def __add__(self, other: "EF") -> "EF":
        return EF([self.value[i] + other.value[i] for i in range(DEG)])

    def __sub__(self, other: "EF") -> "EF":
        return EF([self.value[i] - other.value[i] for i in range(DEG)])

    def __mul__(self, other: "EF") -> "EF":
        result = [F(0) for _ in range(DEG)]
        for i in range(DEG):
            for j in range(DEG):
                if i + j < DEG:
                    result[i + j] += self.value[i] * other.value[j]
                else:
                    result[i + j - DEG] += self.value[i] * other.value[j] * F(W)
        return EF(result)

    def __eq__(self, other: object) -> bool:
        return all(self.value[i] == other.value[i] for i in range(DEG))

    def __repr__(self) -> str:
        return f"EF({self.value})"

    @staticmethod
    def zero() -> "F":
        return EF([F(0) for _ in range(DEG)])
    
    @staticmethod
    def one() -> "F":
        return EF([F(1)] + [F(0) for _ in range(DEG - 1)])


def list_to_base_field(list: Union[List[F], List[EF]]) -> List[F]:
    if len(list) == 0:
        return []
    if isinstance(list[0], EF):
        return sum([e.value for e in list], [])
    else:
        return list


def list_to_ext_field(list: List[F]) -> EF:
    assert len(list) % DEG == 0
    return [EF(list[i:i + DEG]) for i in range(0, len(list), DEG)]
