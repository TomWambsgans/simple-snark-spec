from typing import Literal, Optional, List, Sequence
from finite_field import *
from polynomial import *
from whir import *
import copy


class AirTable:
    n_columns: int
    log_n_rows: int
    constraints: List[ArithmeticCircuit]  # each circuit has 2 * n_columns inputs
    max_constraint_degree: int
    preprocessed_columns: List[List[F]]
    whir_params: WhirParams

    def n_witness_columns(self) -> int:
        return self.n_columns - len(self.preprocessed_columns)

    def log_n_witness_columns(self) -> int:
        # rounded up
        return self.n_witness_columns().bit_length() - 1


def column_up(col: List[F]) -> List[F]:
    up = copy.deepcopy(col)
    up[-1] = up[-2]
    return up


def column_down(col: List[F]) -> List[F]:
    down = copy.deepcopy(col[1:])
    down.append(down[-1])
    return down


def sumcheck_verify(
    fs: FiatShamirVerifier,
    degree: int,
    n_vars: int
) -> Tuple[EF, Evaluation]:  # (sum, delayed evaluation)
    challenges = []
    sum = None
    target = None
    for i in range(n_vars):
        poly = UnivariatePolynomial(fs.receive_scalars_ext(degree + 1))
        if i == 0:
            sum = poly.evaluate(EF.zero()) + poly.evaluate(EF.one())
        else:
            assert target == poly.evaluate(EF.zero()) + poly.evaluate(EF.one())
        challenge = fs.random_scalar()
        challenges.append(challenge)
        target = poly.evaluate(challenge)
    return sum, Evaluation(challenges, target)


def piop_verify(table: AirTable, proof_transcript: List[F]):
    fs = FiatShamirVerifier(proof_transcript)
    whir_commitment = whir_parse_commitment(table.whir_params, fs)
    constraints_batching_scalar = fs.random_scalar()
    zerocheck_challenges = [fs.random_scalar() for _ in range(table.log_n_rows)]
    (zero_sum, zerocheck_eval) = sumcheck_verify(fs, table.max_constraint_degree+1, table.log_n_rows)
    assert zero_sum == 0
    witness_shifted_evals = fs.receive_scalars_ext(table.n_witness_columns() * 2)
    witness_up = witness_shifted_evals[:table.n_witness_columns()]
    witness_down = witness_shifted_evals[table.n_witness_columns():]
    preprocessed_up = [MultilinearEvals(column_up(col)).evaluate(zerocheck_eval.point) for col in table.preprocessed_columns]
    preprocessed_down = [MultilinearEvals(column_down(col)).evaluate(zerocheck_eval.point) for col in table.preprocessed_columns]
    global_point = preprocessed_up + witness_up + preprocessed_down + witness_down
    global_constraint_eval = EF.zero()
    for i, constraint in enumerate(table.constraints):
        global_constraint_eval += constraints_batching_scalar ** i * constraint.evaluate(global_point)
    assert global_constraint_eval * eq_extension(zerocheck_challenges, zerocheck_eval.point) == zerocheck_eval.value

    secondary_sumcheck_batching_scalar = fs.random_scalar()
    batched_inner_sum, inner_sumcheck_challenge = sumcheck_verify(fs, 3, table.log_n_rows + 1)
    assert batched_inner_sum == sum([e * secondary_sumcheck_batching_scalar ** i for i, e in enumerate(witness_shifted_evals)])

    matrix_lde_point = inner_sumcheck_challenge[:1] + zerocheck_eval.point[1:] + inner_sumcheck_challenge[1:]
    matrix_up_eval = ArithmeticCircuit.matrix_up_lde(table.log_n_rows).evaluate(matrix_lde_point)
    matrix_down_eval = ArithmeticCircuit.matrix_down_lde(table.log_n_rows).evaluate(matrix_lde_point)

    final_inner_claims = fs.receive_scalars_ext(table.n_witness_columns())
    batched_inner_value = EF.zero()
    for u in range(table.n_witness_columns()):
        batched_inner_value += final_inner_claims[u] * (secondary_sumcheck_batching_scalar ** u * matrix_up_eval +
                                                        secondary_sumcheck_batching_scalar ** (u + table.n_witness_columns()) * matrix_down_eval)
    batched_inner_value *= MultilinearEvals([F(1) - zerocheck_eval[0], zerocheck_eval[0]]).evaluate(inner_sumcheck_challenge.point[:1])
    assert batched_inner_value == inner_sumcheck_challenge.value

    final_random_scalars = [fs.random_scalar() for _ in range(table.log_n_witness_columns())]
    final_point = final_random_scalars + inner_sumcheck_challenge.point[1:]
    packed_value = MultilinearEvals(final_inner_claims + [EF.zero()
                                    for _ in 2**table.log_n_witness_columns - table.n_witness_columns]).evaluate(final_random_scalars)

    whir_verify(table.whir_params, fs, whir_commitment, Evaluation(final_point, packed_value))
