from typing import Tuple, List
from finite_field import *
from poseidon2 import poseidon2_permutation

Digest = Tuple[Fp, Fp]

def verify_merkle_path(root: Digest, leaf: List[Fp], auth_path: List[(Digest, bool)], height: int):
    # 0) height consistency check
    assert len(auth_path) == height
    # 1) hash the leaf
    state = [Fp(0) for _ in range(4)]
    for i in range(0, len(leaf), 2):
        state[2] = leaf[i]
        state[3] = leaf[i+1] if i+1 < len(leaf) else Fp(0)
        state = poseidon2_permutation(state)
    # 2) walk up the path
    for i in range(height):
        digest, is_left = auth_path[i]
        if is_left:
            state[0] = digest[0]
            state[1] = digest[1]
        else:
            state[2] = digest[0]
            state[3] = digest[1]
        state = poseidon2_permutation(state)
    # 3) check the root
    assert (state[0], state[1]) == root

