"""EvalComp-style bootstrapping built from homomorphic comparison functions.

EvalComp replaces the EvalMod/EvalRound modular-function approximation in CKKS
bootstrapping with a rounding function assembled from Homomorphic Comparison
Functions (HCF).  In this compact reproduction, a ciphertext encodes a scalar
``x = m + q * I`` and EvalComp removes the encrypted integer multiple ``q * I``
by homomorphically estimating ``round(x / q)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, TYPE_CHECKING

from .sign_methods import (
    CipherSignMethod,
    CompositePolynomialSign,
    OptimizedCipherSign,
    sign_method_from_name,
)

if TYPE_CHECKING:  # pragma: no cover - imported only for type checkers
    from .ckks import CKKSCiphertext, NativeCKKSContext


@dataclass(frozen=True)
class EvalCompParameters:
    """Parameters for HCF-based EvalComp rounding.

    ``integer_bound`` is the maximum absolute integer multiple to remove from
    ``x / modulus_step``.  For example, bound 2 supports encrypted values close
    to ``m + qI`` for ``I ∈ {-2,-1,0,1,2}``.
    """

    modulus_step: float = 1.0
    integer_bound: int = 2
    sign_iterations: int = 2
    sign_n: int = 1
    optimized_sign: bool = False
    optimized_sign_degrees: tuple[int, ...] = (3, 5)
    optimized_precision_alpha: int | None = 8
    hcf_domain_radius: float | None = None

    def __post_init__(self) -> None:
        if self.modulus_step <= 0:
            raise ValueError("modulus_step must be positive")
        if self.integer_bound < 1:
            raise ValueError("integer_bound must be at least one")
        if self.sign_iterations < 1:
            raise ValueError("sign_iterations must be at least one")
        if self.sign_n < 1:
            raise ValueError("sign_n must be at least one")
        if not self.optimized_sign_degrees:
            raise ValueError("optimized_sign_degrees cannot be empty")
        invalid_degrees = [
            degree for degree in self.optimized_sign_degrees if degree not in {3, 5}
        ]
        if invalid_degrees:
            raise ValueError(
                f"optimized_sign_degrees must contain only 3 and/or 5: "
                f"{invalid_degrees}"
            )
        if self.optimized_precision_alpha is not None and self.optimized_precision_alpha < 2:
            raise ValueError("optimized_precision_alpha must be at least 2")
        if self.hcf_domain_radius is not None and self.hcf_domain_radius <= 0:
            raise ValueError("hcf_domain_radius must be positive")

    @property
    def domain_radius(self) -> float:
        if self.hcf_domain_radius is not None:
            return self.hcf_domain_radius
        return self.integer_bound + 0.5

    def sign_method(self) -> CipherSignMethod:
        if self.optimized_sign:
            return OptimizedCipherSign(
                iterations=self.sign_iterations,
                degrees=self.optimized_sign_degrees,
                precision_alpha=self.optimized_precision_alpha,
            )
        return CompositePolynomialSign(n=self.sign_n, iterations=self.sign_iterations)


class EvalCompBootstrapper:
    """HCF/EvalComp evaluator for the native CKKS backend."""

    def __init__(
        self,
        context: "NativeCKKSContext",
        params: EvalCompParameters | None = None,
        sign_method: CipherSignMethod | None = None,
    ) -> None:
        self.context = context
        self.params = params or EvalCompParameters()
        self.sign_method = sign_method or self.params.sign_method()

    @classmethod
    def with_sign_method(
        cls,
        context: "NativeCKKSContext",
        method: str,
        params: EvalCompParameters | None = None,
    ) -> "EvalCompBootstrapper":
        params = params or EvalCompParameters()
        sign_method = sign_method_from_name(
            method,
            iterations=params.sign_iterations,
            n=params.sign_n,
            optimized_degrees=params.optimized_sign_degrees,
            optimized_alpha=params.optimized_precision_alpha,
        )
        return cls(context, params=params, sign_method=sign_method)

    def hcf(
        self, ciphertext: "CKKSCiphertext", threshold: float = 0.0
    ) -> "CKKSCiphertext":
        """Evaluate HCF(x, threshold) ≈ 1[x >= threshold] on ciphertexts."""

        shifted = ciphertext.add_plain(-threshold).mul_plain(
            1.0 / self.params.domain_radius
        )
        return self.sign_method.sign(self.context, shifted).mul_plain(
            Fraction(1, 2)
        ).add_plain(Fraction(1, 2))

    def endpoint_nodes(self) -> list[float]:
        """Return EvalComp endpoint vector Node from Algorithm 1."""

        step = self.params.modulus_step
        bound = self.params.integer_bound
        return [(index - 0.5) * step for index in range(-bound, bound + 2)]

    def midpoint_constants(self) -> list[float]:
        """Return EvalComp midpoint vector Mid from Algorithm 2."""

        step = self.params.modulus_step
        bound = self.params.integer_bound
        return [index * step for index in range(-bound, bound + 1)]

    def endpoint_comparison_vector(
        self, ciphertext: "CKKSCiphertext"
    ) -> list["CKKSCiphertext"]:
        """Compute comp(ct, endpoint_i) for all EvalComp endpoints."""

        return [self.hcf(ciphertext, endpoint) for endpoint in self.endpoint_nodes()]

    def interval_judgement_vector(
        self, ciphertext: "CKKSCiphertext"
    ) -> list["CKKSCiphertext"]:
        """Compute redge_i = comp_i - comp_{i+1} from Algorithm 1."""

        comparisons = self.endpoint_comparison_vector(ciphertext)
        return [
            comparisons[index] - comparisons[index + 1]
            for index in range(len(comparisons) - 1)
        ]

    def encrypted_round(self, ciphertext: "CKKSCiphertext") -> "CKKSCiphertext":
        """Approximate qI via EvalComp interval judgement and Mid inner product."""

        redge = self.interval_judgement_vector(ciphertext)
        midpoints = self.midpoint_constants()
        weighted = [
            selector.mul_plain(midpoint)
            for selector, midpoint in zip(redge, midpoints)
        ]
        return _sum_ciphertexts(weighted)

    def eval_round(self, ciphertext: "CKKSCiphertext") -> "CKKSCiphertext":
        """Return an encrypted approximation of q * round(x / q)."""

        return self.encrypted_round(ciphertext)

    def bootstrap(self, ciphertext: "CKKSCiphertext") -> "CKKSCiphertext":
        """Remove the estimated q-multiple from a ciphertext: x - q round(x/q)."""

        return ciphertext - self.eval_round(ciphertext)

    def bootstrap_vector(
        self, ciphertexts: Iterable["CKKSCiphertext"]
    ) -> list["CKKSCiphertext"]:
        """Apply EvalComp independently to scalar-slot ciphertexts."""

        return [self.bootstrap(ciphertext) for ciphertext in ciphertexts]


def _sum_ciphertexts(ciphertexts: list["CKKSCiphertext"]) -> "CKKSCiphertext":
    if not ciphertexts:
        raise ValueError("cannot sum an empty ciphertext list")
    out = ciphertexts[0]
    for ciphertext in ciphertexts[1:]:
        out = out + ciphertext
    return out


__all__ = ["EvalCompBootstrapper", "EvalCompParameters"]
