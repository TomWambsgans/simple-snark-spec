from typing import List, Tuple
from dataclasses import dataclass

from fiat_shamir import FiatShamirVerifier
from finite_field import Fp
from polynomial import *
from merkle_tree import verify_merkle_path


@dataclass
class RoundParams:
    n_variables: int  # before folding
    domain_size: int  # before folding, in log2
    folding_factor: int
    ood_samples: int  # on the folded polynomial
    num_queries: int  # on the original polynomial
    combination_pow_bits: int
    folding_pow_bits: int


@dataclass
class WhirParams:
    initial_ood_samples: int
    rounds: List[RoundParams]
    final_queries: int
    final_sumcheck_rounds: int
    final_combination_pow_bits: int
    final_folding_pow_bits: int


@dataclass
class ParsedCommitment:
    merkle_root: Tuple[Fp, Fp]
    ood_points: List[List[Fp]]
    ood_answers: List[Fp]


def whir_parse_commitment(params: WhirParams, fs: FiatShamirVerifier) -> ParsedCommitment:
    merkle_root = fs.receive_scalars(2)
    ood_points = [multilinear_point_from_univariate(fs.random_scalar()) for _ in params.initial_ood_samples]
    ood_answers = [fs.receive_scalars(1)[0] for _ in params.initial_ood_samples]
    return ParsedCommitment(merkle_root, ood_points, ood_answers)


def whir_verify(params: WhirParams, fs: FiatShamirVerifier, commitment: ParsedCommitment, eval: Evaluation):
    assert len(eval.point) == params.rounds[0].n_variables
    evaluation_points = [commitment.ood_points + [eval.point]]
    combination_randomness = []
    expected_evals = commitment.ood_answers + [eval.value]
    merkle_root = commitment.merkle_root
    expected_sumcheck_output = Fp(0)
    all_folding_randomness = []

    for round in params.rounds:

        # 0. Combination randomness
        fs.pow_grinding(round.combination_pow_bits)
        combination_randomness_gen = fs.random_scalar()
        expected_sumcheck_output += [r * combination_randomness_gen ** i for i, r in enumerate(expected_evals)]
        folding_randomness = []

        # 1. Sumcheck rounds
        for _ in range(round.folding_factor):
            sumcheck_poly = UnivariatePolynomial(fs.receive_scalars(3))
            assert sumcheck_poly.evaluate(0) + sumcheck_poly.evaluate(1) == expected_sumcheck_output
            randomness = fs.random_scalar()
            expected_sumcheck_output = sumcheck_poly.evaluate(randomness)
            folding_randomness.append(randomness)
            fs.pow_grinding(round.folding_pow_bits)

        # 2. Receive folded function
        folded_merkle_root = fs.receive_scalars(2)

        # 3. Out-of-domain sample
        ood_points = [multilinear_point_from_univariate(fs.random_scalar()) for _ in round.ood_samples]

        # 4. Out-of-domain answers
        ood_answers = [fs.receive_scalars(1)[0] for _ in round.ood_samples]

        # 5. Shift queries
        group_gen = Fp.two_addic_generator(round.domain_size - round.folding_factor)
        z_is = []
        folded_evals = []
        for _ in range(round.num_queries):
            index = fs.random_index(round.domain_size - round.folding_factor)
            z_i = group_gen ** index
            merkle_branch = fs.receive_scalars(2 ** round.folding_factor + round.domain_size)
            (leaf, auth_path) = merkle_branch[:2 ** round.folding_factor], merkle_branch[2 ** round.folding_factor:]
            verify_merkle_path(merkle_root, index, leaf, auth_path, round.domain_size)
            folded_eval = MultilinearCoeffs(leaf).evaluate(folding_randomness)
            z_is.append(multilinear_point_from_univariate(z_i))
            folded_evals.append(folded_eval)

        merkle_root = folded_merkle_root
        expected_evals = ood_answers + folded_evals
        all_folding_randomness += folding_randomness
        evaluation_points.append(ood_points + z_is)
        combination_randomness.append(combination_randomness_gen)

    # For simplicity, we do not use the trick of sending the polynomial earlier
    # We assume that it is folded until it becomes constant

    claimed_constant_poly = fs.receive_scalars(1)[0]
    verify_merkle_path(merkle_root, 0, [claimed_constant_poly], [], 0)

    expected_constant_poly = Fp(0)
    for eval_points, combination_randomness_gen in zip(evaluation_points, combination_randomness):
        for i, eval_point in enumerate(eval_points):
            expected_constant_poly += eq_extension(eval_point, all_folding_randomness[-len(eval_point):]) * combination_randomness_gen ** i

    assert expected_constant_poly == claimed_constant_poly
