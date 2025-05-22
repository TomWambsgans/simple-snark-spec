from typing import List
from finite_field import Fp

class UnivariatePolynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
    
    def evaluate(self, x):
        result = 0
        for i, coeff in enumerate(self.coefficients):
            result += coeff * (x ** i)
        return result

class MultilinearPolynomial:
    def __init__(self, coefficients: List[Fp]):
        self.coefficients = coefficients
    
    def evaluate(self, x: List[Fp]) -> Fp:
        result = Fp(0)
        for i, coeff in enumerate(self.coefficients):
            term = coeff
            for j, x_j in enumerate(x):
                term *= x_j ** ((i >> j) & 1)
            result += term
        return result
    
  
def multilinear_point_from_univariate(point: Fp, num_variables: int) -> List[Fp]:
    res = []
    cur = point
    for _ in range(num_variables):
        res.append(cur)
        cur = cur * cur
    return res

def eq_extension(s1: List[Fp], s2: List[Fp]) -> Fp:
    assert len(s1) == len(s2)
    if not s1:
        return Fp(1)
    result = Fp(1)
    for i in range(len(s1)):
        result *= s2[i] * s1[i] + (Fp(1) - s2[i]) * (Fp(1) - s1[i])
    return result