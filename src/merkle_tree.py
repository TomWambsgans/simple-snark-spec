from typing import Tuple, List
from finite_field import *
from poseidon2 import poseidon2_permutation

Digest = Tuple[Fp, Fp]

def verify_merkle_path(root: Digest, index: int, leaf: List[Fp], auth_path: List[Digest], height: int):
    # 1. height consistency check
    assert len(auth_path) == height
    # 2. hash the leaf
    state = [Fp(0) for _ in range(4)]
    for i in range(0, len(leaf), 2):
        state[2] = leaf[i]
        state[3] = leaf[i+1] if i+1 < len(leaf) else Fp(0)
        state = poseidon2_permutation(state)
    # 3. walk up the path
    for i in range(height):
        is_left = (index >> i) & 1 # doable via bit-decomposition for the recursion
        if is_left:
            state[0] = auth_path[i][0]
            state[1] = auth_path[i][1]
        else:
            state[2] = auth_path[i][0]
            state[3] = auth_path[i][1]
        state = poseidon2_permutation(state)
    # 4. check the root
    assert (state[0], state[1]) == root

