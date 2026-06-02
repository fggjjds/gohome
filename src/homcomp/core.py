"""Core composite-polynomial comparison algorithms.

The paper's main comparison method is based on the odd polynomial family

    f_n(x) = c_n ∫_0^x (1 - s^2)^n ds,

where c_n is chosen so that f_n(1)=1.  Repeated self-composition of f_n
approximates the sign function on [-1,-eps] ∪ [eps,1], and comparison follows
from comp(a,b) ≈ (f_n^(d)(a-b)+1)/2.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, comb, log
from typing import Iterable


@dataclass(frozen=True)
class OddPolynomial:
    """Odd polynomial represented as x * P(x^2).

    ``coefficients[k]`` is the coefficient of x^(2k+1).  Storing only odd
    powers mirrors the sign-approximation polynomials used by the paper and
    makes evaluation compact for homomorphic backends.
    """

    coefficients: tuple[Fraction, ...]

    @property
    def degree(self) -> int:
        """Return the ordinary polynomial degree."""

        return 2 * (len(self.coefficients) - 1) + 1

    def __call__(self, x: float | Fraction) -> float | Fraction:
        """Evaluate the polynomial by Horner's rule in x²."""

        y = x * x
        acc: float | Fraction = self.coefficients[-1]
        for coeff in reversed(self.coefficients[:-1]):
            acc = coeff + y * acc
        return x * acc

    def as_float_coefficients(self) -> list[float]:
        """Return coefficients as floats for display or plotting."""

        return [float(c) for c in self.coefficients]

    def format(self, precision: int = 12) -> str:
        """Return a human-readable polynomial expression."""

        terms: list[str] = []
        for k, coeff in enumerate(self.coefficients):
            power = 2 * k + 1
            value = float(coeff)
            if abs(value) < 10 ** (-(precision - 2)):
                continue
            variable = "x" if power == 1 else f"x^{power}"
            terms.append(f"{value:.{precision}g}*{variable}")
        return " + ".join(terms).replace("+ -", "- ") or "0"


def f_polynomial(n: int) -> OddPolynomial:
    """Construct the paper's normalized polynomial f_n.

    Expanding (1-s²)^n and integrating gives
    f_n(x)=c_n Σ_{i=0}^n (-1)^i binom(n,i) x^(2i+1)/(2i+1), with
    c_n chosen as the reciprocal of the same sum evaluated at x=1.
    Fractions keep the small-n reference coefficients exact.
    """

    if n < 1:
        raise ValueError("n must be at least 1")

    integral_at_one = sum(
        Fraction((-1) ** i * comb(n, i), 2 * i + 1) for i in range(n + 1)
    )
    normalizer = Fraction(1, 1) / integral_at_one
    coefficients = tuple(
        normalizer * Fraction((-1) ** i * comb(n, i), 2 * i + 1)
        for i in range(n + 1)
    )
    return OddPolynomial(coefficients)


def compose_sign(x: float, n: int = 4, d: int = 8) -> float:
    """Approximate sign(x) by d self-compositions of f_n."""

    if not -1.0 <= x <= 1.0:
        raise ValueError("x must lie in [-1, 1] for the reproduced guarantee")
    if d < 0:
        raise ValueError("d must be non-negative")

    poly = f_polynomial(n)
    y = x
    for _ in range(d):
        y = float(poly(y))
    return y


def sign_approx(x: float, n: int = 4, d: int = 8) -> float:
    """Alias for :func:`compose_sign` with a descriptive name."""

    return compose_sign(x, n=n, d=d)


def compare(a: float, b: float, n: int = 4, d: int = 8) -> float:
    """Return the NewComp approximation to comp(a,b).

    The inputs should be scaled into [0,1].  The exact comparison convention is
    1 when a>b, 0 when a<b, and 1/2 when a=b; this routine returns a smooth
    approximation of that value.
    """

    if not (0.0 <= a <= 1.0 and 0.0 <= b <= 1.0):
        raise ValueError("a and b must be scaled into [0, 1]")
    return 0.5 * (compose_sign(a - b, n=n, d=d) + 1.0)


def compare_integers(
    a: int, b: int, bits: int, n: int = 4, d: int | None = None
) -> float:
    """Compare two unsigned integers after scaling them to [0,1].

    If ``d`` is omitted, a conservative default is selected for the grid gap
    epsilon = 1/(2^bits-1) and an alpha equal to ``bits``.
    """

    if bits < 1:
        raise ValueError("bits must be positive")
    upper = (1 << bits) - 1
    if not (0 <= a <= upper and 0 <= b <= upper):
        raise ValueError(f"a and b must be unsigned {bits}-bit integers")
    if d is None:
        d = estimate_newcomp_iterations(epsilon=1 / upper, alpha=bits, n=n)
    return compare(a / upper, b / upper, n=n, d=d)


def exact_compare(a: float, b: float) -> float:
    """Exact comparison with the paper's 1/2 tie convention."""

    if a > b:
        return 1.0
    if a < b:
        return 0.0
    return 0.5


def comparison_error(
    pairs: Iterable[tuple[float, float]], n: int = 4, d: int = 8
) -> float:
    """Return the maximum absolute NewComp error over sample pairs."""

    max_error = 0.0
    for a, b in pairs:
        error = abs(compare(a, b, n=n, d=d) - exact_compare(a, b))
        max_error = max(max_error, error)
    return max_error


def estimate_newcomp_iterations(epsilon: float, alpha: int, n: int = 4) -> int:
    """Heuristic iteration count matching the paper's asymptotic theorem.

    The theorem states that d is O(log(1/epsilon)/log(n) + log(alpha)/log(n)).
    Hidden constants depend on the detailed convergence proof; this practical
    estimator uses the leading terms plus two safety compositions.  It is meant
    for experiments, not for formal parameter selection.
    """

    if not (0.0 < epsilon < 1.0):
        raise ValueError("epsilon must lie in (0, 1)")
    if alpha < 1:
        raise ValueError("alpha must be positive")
    if n < 2:
        # log(1) degenerates in the theorem; n=1 still works but converges more slowly.
        base = 2.0
    else:
        base = float(n)
    estimate = 2 * log(1 / epsilon, base) + log(max(alpha, 2), base) + 2
    return max(1, ceil(estimate))


def sample_grid(size: int) -> list[float]:
    """Return a uniform grid on [0,1] including endpoints."""

    if size < 2:
        raise ValueError("size must be at least 2")
    return [i / (size - 1) for i in range(size)]


def separated_pairs(epsilon: float, grid_size: int = 129) -> list[tuple[float, float]]:
    """Generate [0,1]² grid pairs separated by at least epsilon."""

    grid = sample_grid(grid_size)
    return [(a, b) for a in grid for b in grid if abs(a - b) >= epsilon]
