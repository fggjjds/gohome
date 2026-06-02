from fractions import Fraction

import pytest

from homcomp.ckks import CKKSParameters, NativeCKKSContext
from homcomp.evalcomp import EvalCompBootstrapper, EvalCompParameters
from homcomp.sign_methods import (
    OptimizedCipherSign,
    generate_dynamic_polynomial_stages,
    optimized_sign_polynomial,
    solve_optimized_sign_coefficients,
)


def context():
    return NativeCKKSContext(
        CKKSParameters(
            poly_degree=8,
            coefficient_modulus_bits=4096,
            scale_bits=16,
            error_bound=0,
        ),
        seed=2026,
    )


def test_optimized_cipher_sign_coefficients_are_solved_not_hardcoded():
    cubic = solve_optimized_sign_coefficients(3)
    quintic = solve_optimized_sign_coefficients(5)

    assert cubic == (Fraction(3, 2), Fraction(-1, 2))
    assert quintic == (Fraction(15, 8), Fraction(-10, 8), Fraction(3, 8))


def test_optimized_cipher_sign_jointly_uses_dynamic_cubic_and_quintic():
    method = OptimizedCipherSign()

    assert method.degrees == (3, 5)
    assert [stage.degree for stage in method.stages] == [3, 5]
    assert method.stages[-1].target_error < method.stages[0].source_error


def test_dynamic_stage_generator_reduces_error_to_alpha_threshold():
    stages = generate_dynamic_polynomial_stages(0.5, precision_alpha=8)

    assert [stage.degree for stage in stages] == [3, 5]
    assert stages[-1].target_error <= 2 ** (1 - 8)


def test_optimized_cipher_sign_polynomial_is_endpoint_preserving():
    for degree in (3, 5):
        polynomial = optimized_sign_polynomial(degree)

        assert polynomial.degree == degree
        assert polynomial(1) == 1
        assert polynomial(-1) == -1
        assert polynomial(0) == 0


def test_evalcomp_hcf_baseline_comparison_function():
    ctx = context()
    bootstrapper = EvalCompBootstrapper(
        ctx,
        EvalCompParameters(
            integer_bound=1,
            sign_iterations=1,
            optimized_sign=False,
            hcf_domain_radius=1.0,
        ),
    )

    assert bootstrapper.hcf(ctx.encrypt(0.75), 0.5).decrypt() > 0.65
    assert bootstrapper.hcf(ctx.encrypt(0.25), 0.5).decrypt() < 0.35


def test_optimized_sign_replaces_evalcomp_hcf_component():
    ctx = context()
    baseline = EvalCompBootstrapper(
        ctx,
        EvalCompParameters(
            integer_bound=1,
            sign_iterations=1,
            optimized_sign=False,
            hcf_domain_radius=1.0,
        ),
    )
    optimized = EvalCompBootstrapper(
        ctx,
        EvalCompParameters(
            integer_bound=1,
            sign_iterations=1,
            optimized_sign=True,
            optimized_sign_degrees=(3, 5),
            hcf_domain_radius=1.0,
        ),
    )

    baseline_gap = baseline.hcf(ctx.encrypt(0.75), 0.5).decrypt() - baseline.hcf(
        ctx.encrypt(0.25), 0.5
    ).decrypt()
    optimized_gap = optimized.hcf(ctx.encrypt(0.75), 0.5).decrypt() - optimized.hcf(
        ctx.encrypt(0.25), 0.5
    ).decrypt()

    assert optimized_gap > baseline_gap


def test_evalcomp_extracts_integer_multiple_and_bootstraps_ciphertext():
    ctx = context()
    bootstrapper = EvalCompBootstrapper(
        ctx,
        EvalCompParameters(
            modulus_step=1.0,
            integer_bound=1,
            sign_iterations=1,
            optimized_sign=True,
            optimized_sign_degrees=(3, 5),
            hcf_domain_radius=1.5,
        ),
    )
    ciphertext = ctx.encrypt(0.25)

    assert bootstrapper.eval_round(ciphertext).decrypt() == pytest.approx(
        0.25, abs=2e-2
    )
    assert bootstrapper.bootstrap(ciphertext).decrypt() == pytest.approx(
        0.0, abs=2e-2
    )
