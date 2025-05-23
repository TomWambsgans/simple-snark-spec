from typing import List, Tuple
from dataclasses import dataclass

from fiat_shamir import FiatShamirVerifier
from finite_field import F
from polynomial import *
from merkle_tree import verify_merkle_path
from poseidon2 import *


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


@dataclass
class ParsedCommitment:
    merkle_root: Digest
    ood_points: List[List[EF]]
    ood_answers: List[EF]


def whir_parse_commitment(params: WhirParams, fs: FiatShamirVerifier) -> ParsedCommitment:
    merkle_root = fs.receive_scalars_base(DIGEST_LEN)
    ood_points = [multilinear_point_from_univariate(fs.random_scalar()) for _ in params.initial_ood_samples]
    ood_answers = [fs.receive_scalars_ext(1)[0] for _ in params.initial_ood_samples]
    return ParsedCommitment(merkle_root, ood_points, ood_answers)


def whir_verify(params: WhirParams, fs: FiatShamirVerifier, commitment: ParsedCommitment, eval: Evaluation):
    assert len(eval.point) == params.rounds[0].n_variables
    evaluation_points = [commitment.ood_points + [eval.point]]
    combination_randomness = []
    expected_evals = commitment.ood_answers + [eval.value]
    merkle_root = commitment.merkle_root
    expected_sumcheck_output = EF.zero()
    all_folding_randomness = []

    for round in params.rounds:

        # 0. Combination randomness
        fs.pow_grinding(round.combination_pow_bits)
        combination_randomness_gen = fs.random_scalar()
        expected_sumcheck_output += [r * combination_randomness_gen ** i for i, r in enumerate(expected_evals)]
        folding_randomness = []

        # 1. Sumcheck rounds
        for _ in range(round.folding_factor):
            sumcheck_poly = UnivariatePolynomial(fs.receive_scalars_ext(3))
            assert sumcheck_poly.evaluate(EF.zero()) + sumcheck_poly.evaluate(EF.one()) == expected_sumcheck_output
            randomness = fs.random_scalar()
            expected_sumcheck_output = sumcheck_poly.evaluate(randomness)
            folding_randomness.append(randomness)
            fs.pow_grinding(round.folding_pow_bits)

        # 2. Receive folded function
        folded_merkle_root = fs.receive_scalars_base(DIGEST_LEN)

        # 3. Out-of-domain sample
        ood_points = [multilinear_point_from_univariate(fs.random_scalar()) for _ in round.ood_samples]

        # 4. Out-of-domain answers
        ood_answers = [fs.receive_scalars_ext(1)[0] for _ in round.ood_samples]

        # 5. Shift queries
        query_domain = round.domain_size - round.folding_factor
        group_gen = EF.from_base(F.two_addic_generator(query_domain))
        z_is = []
        folded_evals = []
        for _ in range(round.num_queries):
            index = fs.random_index(query_domain)
            z_i = group_gen ** index
            leaf: Union[List[F], List[EF]] = fs.receive_scalars_base(
                2 ** round.folding_factor) if round == 0 else fs.receive_scalars_ext(2 ** round.folding_factor)
            auth_path = fs.receive_scalars_base(query_domain * DIGEST_LEN)
            verify_merkle_path(merkle_root, index, leaf, auth_path, query_domain)
            folded_eval = MultilinearCoeffs(list_to_ext_field(leaf)).evaluate(folding_randomness)
            z_is.append(multilinear_point_from_univariate(z_i))
            folded_evals.append(folded_eval)

        merkle_root = folded_merkle_root
        expected_evals = ood_answers + folded_evals
        all_folding_randomness += folding_randomness
        evaluation_points.append(ood_points + z_is)
        combination_randomness.append(combination_randomness_gen)

    # For simplicity, we do not use the trick of sending the polynomial earlier
    # We assume that it is folded until it becomes constant

    claimed_constant_poly = fs.receive_scalars_ext(1)[0]
    verify_merkle_path(merkle_root, 0, [claimed_constant_poly], [], 0)

    expected_constant_poly = EF.zero()
    for eval_points, combination_randomness_gen in zip(evaluation_points, combination_randomness):
        for i, eval_point in enumerate(eval_points):
            expected_constant_poly += eq_extension(eval_point, all_folding_randomness[-len(eval_point):]) * combination_randomness_gen ** i

    assert expected_constant_poly == claimed_constant_poly
