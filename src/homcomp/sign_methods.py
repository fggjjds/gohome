"""Composable ciphertext sign-function evaluators.

The EvalComp bootstrapping paper builds its rounding step from a Homomorphic
Comparison Function (HCF).  This module separates the HCF sign approximation
from the surrounding rounding logic so that the comparison-function component can
be replaced by alternative ciphertext sign-function methods.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import prod
from typing import Protocol, TYPE_CHECKING

from .core import OddPolynomial

if TYPE_CHECKING:  # pragma: no cover - imported only for type checkers
    from .ckks import CKKSCiphertext, NativeCKKSContext


class CipherSignMethod(Protocol):
    """Protocol for encrypted sign approximators used by EvalComp HCF."""

    name: str

    def sign(
        self, context: "NativeCKKSContext", ciphertext: "CKKSCiphertext"
    ) -> "CKKSCiphertext":
        """Return an encrypted approximation to sign(ciphertext)."""


@dataclass(frozen=True)
class CompositePolynomialSign:
    """EvalComp baseline HCF sign using the Cheon-Kim-Kim polynomial family."""

    n: int = 1
    iterations: int = 2
    name: str = "evalcomp-composite"

    def sign(
        self, context: "NativeCKKSContext", ciphertext: "CKKSCiphertext"
    ) -> "CKKSCiphertext":
        return context.sign(ciphertext, n=self.n, d=self.iterations)


@dataclass(frozen=True)
class OptimizedCipherSign:
    """Optimized cipher-symbol sign method with solved polynomial coefficients.

    The optimized paper derives odd cubic/quintic sign polynomials by solving a
    small polynomial-equation system rather than by tuning coefficients manually.
    The optimized evaluator jointly uses the cubic and quintic stages by default:
    it solves both systems and composes the resulting polynomials on ciphertexts.
    """

    iterations: int = 2
    degrees: tuple[int, ...] = field(default_factory=lambda: (3, 5))
    precision_alpha: int | None = 8
    initial_error: float = 0.5
    cubic_count: int = 1
    name: str = "optimized-cipher-sign"

    def __post_init__(self) -> None:
        if not self.degrees:
            raise ValueError("at least one optimized sign degree is required")
        invalid = [degree for degree in self.degrees if degree not in {3, 5}]
        if invalid:
            raise ValueError(f"unsupported optimized sign degree(s): {invalid}")
        if self.precision_alpha is not None and self.precision_alpha < 2:
            raise ValueError("precision_alpha must be at least 2")
        if not 0 < self.initial_error < 1:
            raise ValueError("initial_error must lie in (0, 1)")
        if self.cubic_count < 0:
            raise ValueError("cubic_count cannot be negative")

    @property
    def stages(self) -> tuple[DynamicPolynomialStage, ...]:
        if self.precision_alpha is None:
            return tuple(
                DynamicPolynomialStage(
                    degree=degree,
                    source_error=0.0,
                    target_error=0.0,
                    polynomial=optimized_sign_polynomial(degree),
                    extrema=(),
                )
                for degree in self.degrees
            )
        return generate_dynamic_polynomial_stages(
            self.initial_error, self.precision_alpha, self.cubic_count
        )

    @property
    def polynomials(self) -> tuple[OddPolynomial, ...]:
        return tuple(stage.polynomial for stage in self.stages)

    def sign(
        self, context: "NativeCKKSContext", ciphertext: "CKKSCiphertext"
    ) -> "CKKSCiphertext":
        out = ciphertext
        polynomials = self.polynomials
        for _ in range(self.iterations):
            for polynomial in polynomials:
                out = context.eval_polynomial(out, polynomial)
        return out


def solve_optimized_sign_coefficients(degree: int = 5) -> tuple[Fraction, ...]:
    """Solve the cubic/quintic optimized sign polynomial coefficient system.

    For an odd polynomial ``p(x)=Σ a_i x^(2i+1)``, the coefficient equations are

    * ``p(1)=1`` so the positive endpoint is fixed;
    * ``p^(r)(1)=0`` for ``r=1..m`` where ``degree=2m+1``, flattening the curve
      at ±1 and making repeated composition saturate toward the sign function.

    The paper discusses the cubic and quintic cases.  This function supports
    exactly those two systems and solves them by exact Gaussian elimination.
    """

    if degree not in {3, 5}:
        raise ValueError("optimized cipher sign supports only degree 3 or 5")

    max_power_index = (degree - 1) // 2
    powers = [2 * i + 1 for i in range(max_power_index + 1)]
    matrix: list[list[Fraction]] = []
    rhs: list[Fraction] = []

    for derivative_order in range(max_power_index + 1):
        matrix.append(
            [
                Fraction(_falling_factorial(power, derivative_order), 1)
                for power in powers
            ]
        )
        rhs.append(Fraction(1 if derivative_order == 0 else 0, 1))

    return tuple(_solve_linear_system(matrix, rhs))


def optimized_sign_polynomial(degree: int = 5) -> OddPolynomial:
    """Return the solved optimized ciphertext sign polynomial."""

    return OddPolynomial(solve_optimized_sign_coefficients(degree))


def sign_method_from_name(
    name: str,
    iterations: int = 2,
    n: int = 1,
    optimized_degrees: tuple[int, ...] = (3, 5),
    optimized_alpha: int | None = 8,
) -> CipherSignMethod:
    """Build a sign-method object from a CLI-friendly name."""

    normalized = name.strip().lower().replace("_", "-")
    if normalized in {"evalcomp", "baseline", "composite", "newcomp"}:
        return CompositePolynomialSign(n=n, iterations=iterations)
    if normalized in {"optimized", "optimized-cipher-sign", "cipher-symbol"}:
        return OptimizedCipherSign(
            iterations=iterations,
            degrees=optimized_degrees,
            precision_alpha=optimized_alpha,
        )
    raise ValueError(f"unknown sign method: {name}")


def _falling_factorial(value: int, length: int) -> int:
    if length == 0:
        return 1
    return prod(value - offset for offset in range(length))


def _solve_linear_system(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> list[Fraction]:
    """Solve a square linear system by exact Gauss-Jordan elimination."""

    size = len(rhs)
    augmented = [row[:] + [rhs_value] for row, rhs_value in zip(matrix, rhs)]

    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column] != 0),
            None,
        )
        if pivot is None:
            raise ValueError("singular coefficient system")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]

        pivot_value = augmented[column][column]
        augmented[column] = [value / pivot_value for value in augmented[column]]

        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]

    return [row[-1] for row in augmented]

@dataclass(frozen=True)
class DynamicPolynomialStage:
    """One dynamically solved cubic/quintic stage from the optimized method."""

    degree: int
    source_error: float
    target_error: float
    polynomial: OddPolynomial
    extrema: tuple[float, ...]


def solve_dynamic_cubic_stage(error: float) -> DynamicPolynomialStage:
    """Solve the paper's cubic dynamic-polynomial system for one error band."""

    if not 0 < error < 1:
        raise ValueError("error must lie in (0, 1)")
    left = 1 - error
    right = 1 + error
    sum_pair = left * left + left * right + right * right
    extremum = (sum_pair / 3) ** 0.5
    denominator = (extremum**3 - sum_pair * extremum) + (
        left**3 - sum_pair * left
    )
    a = 2 / denominator
    b = -a * sum_pair
    target_error = a * extremum**3 + b * extremum - 1
    return DynamicPolynomialStage(
        degree=3,
        source_error=error,
        target_error=abs(target_error),
        polynomial=_float_odd_polynomial((b, a)),
        extrema=(extremum,),
    )


def solve_dynamic_quintic_stage(error: float) -> DynamicPolynomialStage:
    """Solve the paper's quintic dynamic-polynomial system for one error band."""

    if not 0 < error < 1:
        raise ValueError("error must lie in (0, 1)")
    solution = _solve_quintic_system(error)
    a, b, c, x0, x1, target_error = solution
    return DynamicPolynomialStage(
        degree=5,
        source_error=error,
        target_error=abs(target_error),
        polynomial=_float_odd_polynomial((c, b, a)),
        extrema=(x0, x1),
    )


def generate_dynamic_polynomial_stages(
    initial_error: float,
    precision_alpha: int,
    cubic_count: int = 1,
) -> tuple[DynamicPolynomialStage, ...]:
    """Generate finite cubic+quintic stages following Algorithm 1 in the image.

    The loop first applies a bounded number of cubic stages and then uses
    quintic stages until the dynamic error reaches ``2 ** (1 - alpha)``.
    """

    if precision_alpha < 2:
        raise ValueError("precision_alpha must be at least 2")
    if cubic_count < 0:
        raise ValueError("cubic_count cannot be negative")
    threshold = 2 ** (1 - precision_alpha)
    error = initial_error
    stages: list[DynamicPolynomialStage] = []

    for _ in range(cubic_count):
        if error <= threshold:
            break
        stage = solve_dynamic_cubic_stage(error)
        stages.append(stage)
        error = stage.target_error

    while error > threshold:
        stage = solve_dynamic_quintic_stage(error)
        stages.append(stage)
        error = stage.target_error

    return tuple(stages)


def _float_odd_polynomial(coefficients: tuple[float, ...]) -> OddPolynomial:
    return OddPolynomial(
        tuple(Fraction(value).limit_denominator(10**12) for value in coefficients)
    )


def _solve_quintic_system(error: float) -> tuple[float, float, float, float, float, float]:
    left = 1 - error
    right = 1 + error
    variables = [
        100.0,
        -300.0,
        201.0,
        left + (right - left) * 0.25,
        left + (right - left) * 0.75,
        error / 2,
    ]

    for _ in range(80):
        residual = _quintic_residual(variables, error)
        norm = max(abs(value) for value in residual)
        if norm < 1e-12:
            return tuple(variables)  # type: ignore[return-value]
        jacobian = _finite_difference_jacobian(variables, residual, error)
        step = _solve_float_linear_system(jacobian, [-value for value in residual])
        damping = 1.0
        while damping > 1e-6:
            candidate = [
                value + damping * delta for value, delta in zip(variables, step)
            ]
            if left < candidate[3] < candidate[4] < right:
                candidate_norm = max(
                    abs(value) for value in _quintic_residual(candidate, error)
                )
                if candidate_norm < norm:
                    variables = candidate
                    break
            damping *= 0.5
        else:
            variables = [value + delta for value, delta in zip(variables, step)]
    raise ValueError("quintic dynamic polynomial solve did not converge")


def _quintic_residual(variables: list[float], error: float) -> list[float]:
    a, b, c, x0, x1, target_error = variables
    left = 1 - error
    right = 1 + error

    def h(value: float) -> float:
        return a * value**5 + b * value**3 + c * value

    def derivative(value: float) -> float:
        return 5 * a * value**4 + 3 * b * value**2 + c

    return [
        h(left) - (1 - target_error),
        h(right) - (1 + target_error),
        h(x0) - (1 + target_error),
        h(x1) - (1 - target_error),
        derivative(x0),
        derivative(x1),
    ]


def _finite_difference_jacobian(
    variables: list[float], residual: list[float], error: float
) -> list[list[float]]:
    columns: list[list[float]] = []
    for index, value in enumerate(variables):
        step = max(1e-7, abs(value) * 1e-7)
        shifted = variables[:]
        shifted[index] += step
        shifted_residual = _quintic_residual(shifted, error)
        columns.append(
            [
                (shifted_value - base_value) / step
                for shifted_value, base_value in zip(shifted_residual, residual)
            ]
        )
    return [list(row) for row in zip(*columns)]


def _solve_float_linear_system(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    size = len(rhs)
    augmented = [row[:] + [rhs_value] for row, rhs_value in zip(matrix, rhs)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-18:
            raise ValueError("singular dynamic polynomial system")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [value / pivot_value for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[-1] for row in augmented]


__all__ = [
    "CipherSignMethod",
    "CompositePolynomialSign",
    "OptimizedCipherSign",
    "optimized_sign_polynomial",
    "sign_method_from_name",
    "solve_optimized_sign_coefficients",
    "DynamicPolynomialStage",
    "generate_dynamic_polynomial_stages",
    "solve_dynamic_cubic_stage",
    "solve_dynamic_quintic_stage",
]
