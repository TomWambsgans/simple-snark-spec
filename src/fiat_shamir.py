from typing import List
from finite_field import *
from poseidon2 import *


class FiatShamirVerifier:
    def __init__(self, transcript: List[Fp]):
        self.transcript = transcript
        self.cursor = 0
        self.state = (Fp(0), Fp(0))

    def _update_state(self, scalars: List[Fp]) -> None:
        for i in range(0, len(scalars), 2):
            inputs = (self.state[0], self.state[1],  scalars[i], scalars[i+1] if i+1 < len(scalars) else Fp(0))
            self.state = poseidon2_permutation(inputs)[:2]

    def receive_scalars(self, n: int) -> List[Fp]:
        scalars = self.transcript[self.cursor:self.cursor + n]
        self.cursor += n
        self._update_state(scalars)
        return scalars

    def random_scalar(self) -> Fp:
        challenge = self.state[0]
        self.state = poseidon2_permutation((self.state[0], self.state[1], Fp(0), Fp(0)))[:2]
        return challenge

    def random_bits(self, n: int) -> List[bool]:
        assert (n < P_BITS)
        scalar = self.random_scalar()
        # Not very recursion friendly, requires to decompose a field element into P_BITS individual bits.
        [(scalar.value >> i) & 1 for i in range(n)]

    def pow_grinding(self, bits: int):
        _ = self.random_scalar()  # nonce
        assert self.random_bits(bits) == [0] * bits
