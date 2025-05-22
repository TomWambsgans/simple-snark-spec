from typing import Tuple
from finite_field import *

PermutationState = Tuple[Fp, Fp, Fp, Fp]


def poseidon2_permutation(state: PermutationState) -> PermutationState:
    # TODO
    # Dummy permutation for now:
    s1, s2, s3, s4 = state
    o1 = s1 + s2 * Fp(785) + s3 * Fp(123) + s4 * Fp(456)
    o2 = s1 * Fp(789) + s2 + s3 * Fp(321) + s4 * Fp(654)
    o3 = s1 * Fp(159) + s2 * Fp(753) + s3 + s4 * Fp(951)
    o4 = s1 * Fp(357) + s2 * Fp(159) + s3 * Fp(753) + s4
    return o1*o2, o3*o4, o1*o1, o4*o4
