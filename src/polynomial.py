from dataclasses import dataclass
from typing import Literal, Optional, List, Sequence
from finite_field import *


@dataclass
class Evaluation:
    point: List[EF]
    value: EF


class UnivariatePolynomial:
    def __init__(self, coefficients: List[EF]):
        self.coefficients = coefficients

    def evaluate(self, x: EF):
        result = 0
        for i, coeff in enumerate(self.coefficients):
            result += coeff * (x ** i)
        return result


class MultilinearCoeffs:
    # a multilinear polynomial defined by its coefficients (canonical form)

    def __init__(self, coefficients: List[EF]):
        self.coefficients = coefficients

    def evaluate(self, x: List[EF]) -> EF:
        result = EF.zero()
        for i, coeff in enumerate(self.coefficients):
            term = coeff
            for j, x_j in enumerate(x):
                term *= x_j ** ((i >> j) & 1)
            result += term
        return result


class MultilinearEvals:
    # a multilinear polynomial defined by its evaluations (on the hypercube)

    def __init__(self, evals: List[EF]):
        self.evals = evals

    def evaluate(self, x: List[EF]) -> EF:
        result = EF.zero()
        for i in range(len(self.evals)):
            basis_term = EF.one()
            for j in range(len(x)):
                bit = (i >> j) & 1
                if bit == 1:
                    basis_term *= x[j]
                else:
                    basis_term *= (EF.one() - x[j])
            result += self.evals[i] * basis_term
        return result


def multilinear_point_from_univariate(point: EF, num_variables: int) -> List[EF]:
    res = []
    cur = point
    for _ in range(num_variables):
        res.append(cur)
        cur = cur * cur
    return res


def eq_extension(s1: List[EF], s2: List[EF]) -> EF:
    assert len(s1) == len(s2)
    if not s1:
        return EF.one()
    result = EF.one()
    for i in range(len(s1)):
        result *= s2[i] * s1[i] + (EF.one() - s2[i]) * (EF.one() - s1[i])
    return result


Op = Literal["const", "input", "add", "mul"]


@dataclass
class ArithmeticCircuit:
    op: Op
    children: Sequence["ArithmeticCircuit"] = ()  # for "add" / "mul"
    const_value: Optional[F] = None
    input_index: Optional[int] = None

    @staticmethod
    def const(value: F) -> "ArithmeticCircuit": return ArithmeticCircuit("const", const_value=value)

    @staticmethod
    def var(i: int) -> "ArithmeticCircuit": return ArithmeticCircuit("input", input_index=i)

    def __add__(self, other: "ArithmeticCircuit") -> "ArithmeticCircuit": return ArithmeticCircuit("add", children=(self, other))
    def __mul__(self, other: "ArithmeticCircuit") -> "ArithmeticCircuit": return ArithmeticCircuit("mul", children=(self, other))

    def evaluate(self, inputs: List[EF]) -> EF:
        if self.op == "const":
            return EF.from_base(self.const_value)

        if self.op == "input":
            idx = self.input_index
            return inputs[idx]

        if self.op == "add":
            acc = EF.zero()
            for child in self.children:
                acc += child.evaluate(inputs)
            return acc

        if self.op == "mul":
            acc = EF.one()
            for child in self.children:
                acc *= child.evaluate(inputs)
            return acc

    @staticmethod
    def eq_extension_2n_vars(n: int) -> "ArithmeticCircuit":
        # eq(Xs, Ys) = ((Y0 X0 + (Y0 - 1) (X0 - 1)) * ((Y1 X1 + (Y1 - 1) (X1 - 1)) ...
        left = [ArithmeticCircuit.var(i) for i in range(n)]
        right = [ArithmeticCircuit.var(i + n) for i in range(n)]
        return ArithmeticCircuit._eq_extension(left, right)

    def eq_extension_n_scalars(self, scalars: List[EF]) -> "ArithmeticCircuit":
        # eq(scalars, Xs) = ((scalars[0] X0 + (scalars[0] - 1) (X0 - 1)) * ((scalars[1] X1 + (scalars[1] - 1) (X1 - 1)) ...
        left = [ArithmeticCircuit.var(i) for i in range(len(scalars))]
        right = [ArithmeticCircuit.const(scalar) for scalar in scalars]
        return ArithmeticCircuit._eq_extension(left, right)

    @staticmethod
    def _eq_extension(left: List["ArithmeticCircuit"], right: List["ArithmeticCircuit"]) -> "ArithmeticCircuit":
        assert len(left) == len(right)
        result = ArithmeticCircuit.const(EF.one())
        for l, r in zip(left, right):
            result *= (l * r) + ((ArithmeticCircuit.const(EF.one()) - l) * (ArithmeticCircuit.const(EF.one()) - r))
        return result

    @staticmethod
    def next(n: int) -> "ArithmeticCircuit":
        # returns a polynomial P in 2n vars, where P(x, y) = 1 iif y = x + 1 in big endian (both numbers are n bits)

        def factor(l, r):
            return ArithmeticCircuit.var(l) * (ArithmeticCircuit.const(EF.one()) + ArithmeticCircuit.var(r) * ArithmeticCircuit.const(EF(-1)))

        def g(k):
            factors = []
            # Add factors for high bits that should match
            for i in range(n - k, n):
                factors.append(factor(i, i + n))
            factors.append(factor(2 * n - 1 - k, n - 1 - k))
            if k < n - 1:
                left = [ArithmeticCircuit.var(i) for i in range(n - k - 1)]
                right = [ArithmeticCircuit.var(i + n) for i in range(n - k - 1)]
                factors.append(ArithmeticCircuit._eq_extension(left, right))
            result = ArithmeticCircuit.const(EF.one())
            for f in factors:
                result *= f
            return result

        result = ArithmeticCircuit.const(EF.zero())
        for k in range(n):
            result += g(k)
        return result

    def matrix_up_lde(n: int) -> "ArithmeticCircuit":
        return ArithmeticCircuit.eq_extension_2n_vars(n) + ArithmeticCircuit.eq_extension_n_scalars([EF.one() for _ in (2 * n - 1)]) * (ArithmeticCircuit.const(EF.one()) -
                                                                                                                                     ArithmeticCircuit.var(2 * n - 1) * ArithmeticCircuit.const(EF.from_base(F(2))))

    def matrix_down_lde(n: int) -> "ArithmeticCircuit":
        return ArithmeticCircuit.next(n) + ArithmeticCircuit.eq_extension_n_scalars([EF.one() for _ in (2 * n)])
