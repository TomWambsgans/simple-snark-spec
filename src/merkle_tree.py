from typing import Annotated, List, Literal, Union
from finite_field import *
from poseidon2 import *

Digest = Annotated[List[F], Literal[DIGEST_LEN]]


def verify_merkle_path(root: Digest, index: int, leaf: Union[List[F], List[EF]], auth_path: List[Digest], height: int):
    # 1. height consistency check
    assert len(auth_path) == height
    # 2. hash the leaf
    leaf_base = list_to_base_field(leaf)
    state: PermutationState = [F(0) for _ in range(POSEIDON_WIDTH)]
    for i in range(0, len(leaf), DIGEST_LEN):
        for j in range(DIGEST_LEN):
            state[j] = leaf_base[i+j] if i + j < len(leaf_base) else F(0)
        state = poseidon2_permutation(state)
    # 3. walk up the path
    for i in range(height):
        is_left = (index >> i) & 1  # doable via bit-decomposition for the recursion
        if is_left:
            for j in range(DIGEST_LEN):
                state[j + DIGEST_LEN] = auth_path[i][j]
        else:
            for j in range(DIGEST_LEN):
                state[j] = auth_path[i][j]
        state = poseidon2_permutation(state)
    # 4. check the root
    assert state[:DIGEST_LEN] == root
