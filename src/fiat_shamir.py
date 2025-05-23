from typing import List
from finite_field import *
from poseidon2 import *


class FiatShamirVerifier:
    def __init__(self, transcript: List[F]):
        self.transcript = transcript
        self.cursor = 0
        self.state = [F(0) for _ in range(POSEIDON_WIDTH)]

    def _update_state(self, scalars: List[F]) -> None:
        for i in range(0, len(scalars), DIGEST_LEN):
            for j in range(DIGEST_LEN):
                self.state[j] = scalars[i + j] if i + j < len(scalars) else F(0)
            self.state = poseidon2_permutation(self.state)

    def receive_scalars_base(self, n: int) -> List[F]:
        scalars = self.transcript[self.cursor:self.cursor + n]
        self.cursor += n
        self._update_state(scalars)
        return scalars

    def receive_scalars_ext(self, n: int) -> List[EF]:
        return list_to_ext_field(self.receive_scalars_base(n * DEG))

    def random_scalar(self) -> EF:
        challenge = EF(self.state[:DEG])
        self.state = poseidon2_permutation(self.state)
        return challenge

    def random_index(self, bits: int) -> int:
        # Not very recursion friendly, requires to decompose a field element into individual bits
        assert (bits < P_BITS)
        self.random_scalar()[0].value % (1 << bits)

    def pow_grinding(self, bits: int):
        _ = self.receive_scalars_base(0)  # nonce
        assert self.random_index(bits) == 0
