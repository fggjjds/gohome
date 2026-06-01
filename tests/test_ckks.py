import pytest

from homcomp import compare
from homcomp.ckks import CKKSParameters, NativeCKKSContext


def context(scale_bits=18, modulus_bits=4096):
    return NativeCKKSContext(
        CKKSParameters(
            poly_degree=8,
            coefficient_modulus_bits=modulus_bits,
            scale_bits=scale_bits,
            error_bound=1,
        ),
        seed=1234,
    )


def test_native_ckks_encrypts_adds_and_multiplies_ciphertexts():
    ctx = context()
    a = ctx.encrypt(0.375)
    b = ctx.encrypt(0.25)

    assert (a + b).decrypt() == pytest.approx(0.625, abs=1e-4)
    assert (a - b).decrypt() == pytest.approx(0.125, abs=1e-4)
    assert (a * b).decrypt() == pytest.approx(0.09375, abs=1e-4)


def test_native_ckks_evaluates_encrypted_newcomp():
    ctx = context(scale_bits=20)
    encrypted = ctx.compare(ctx.encrypt(0.75), ctx.encrypt(0.25), n=1, d=2)
    plaintext = compare(0.75, 0.25, n=1, d=2)

    assert encrypted.decrypt() == pytest.approx(plaintext, abs=1e-4)


def test_native_ckks_vector_comparison_runs_slotwise():
    ctx = context(scale_bits=20)
    left = ctx.encrypt_vector([0.75, 0.25])
    right = ctx.encrypt_vector([0.25, 0.75])

    decrypted = left.compare(right, n=1, d=2).decrypt()

    assert decrypted[0] > 0.9
    assert decrypted[1] < 0.1


def test_native_ckks_can_run_degree_nine_paper_polynomial_once():
    ctx = NativeCKKSContext(
        CKKSParameters(
            poly_degree=16,
            coefficient_modulus_bits=8192,
            scale_bits=16,
            error_bound=1,
        ),
        seed=99,
    )

    encrypted = ctx.compare(ctx.encrypt(0.75), ctx.encrypt(0.25), n=4, d=1)
    plaintext = compare(0.75, 0.25, n=4, d=1)

    assert encrypted.decrypt() == pytest.approx(plaintext, abs=1e-3)
