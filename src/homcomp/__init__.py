"""Composite-polynomial homomorphic comparison routines.

The package includes the paper polynomial formulas plus a native CKKS/HEAAN-style
ciphertext backend for encrypted end-to-end reproduction.
"""

from .ckks import (
    CKKSParameters,
    CKKSCiphertext,
    CKKSEncryptedVector,
    NativeCKKSContext,
)
from .evalcomp import EvalCompBootstrapper, EvalCompParameters
from .paper_results import format_paper_report
from .tenseal_backend import TenSEALCKKSContext, TenSEALCKKSParameters, TenSEALCiphertext
from .sign_methods import (
    CompositePolynomialSign,
    OptimizedCipherSign,
    solve_optimized_sign_coefficients,
)
from .core import (
    OddPolynomial,
    comparison_error,
    compare,
    compare_integers,
    compose_sign,
    estimate_newcomp_iterations,
    f_polynomial,
    sign_approx,
)

__all__ = [
    "CKKSParameters",
    "CKKSCiphertext",
    "CKKSEncryptedVector",
    "NativeCKKSContext",
    "TenSEALCKKSContext",
    "TenSEALCKKSParameters",
    "TenSEALCiphertext",
    "EvalCompBootstrapper",
    "EvalCompParameters",
    "CompositePolynomialSign",
    "OptimizedCipherSign",
    "solve_optimized_sign_coefficients",
    "OddPolynomial",
    "comparison_error",
    "compare",
    "compare_integers",
    "compose_sign",
    "estimate_newcomp_iterations",
    "f_polynomial",
    "format_paper_report",
    "sign_approx",
]
