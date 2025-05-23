from finite_field import *
from typing import Annotated, List, Literal

POSEIDON_WIDTH = 16
DIGEST_LEN = POSEIDON_WIDTH // 2

PermutationState = Annotated[List[F], Literal[POSEIDON_WIDTH]]
Digest = Annotated[List[F], Literal[DIGEST_LEN]]


def poseidon2_permutation(state: PermutationState) -> PermutationState:
    # TODO
    pass