"""Optional TenSEAL backend for running the circuits on a real CKKS library.

This module is intentionally separate from :mod:`homcomp.ckks`: the native
backend is kept for deterministic unit tests, while this adapter delegates
ciphertext storage and arithmetic to TenSEAL (Microsoft SEAL underneath) when the
``tenseal`` package is installed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence

from .core import OddPolynomial, f_polynomial


try:  # pragma: no cover - availability depends on the user's environment
    import tenseal as ts
except ImportError:  # pragma: no cover
    ts = None  # type: ignore[assignment]


@dataclass(frozen=True)
class TenSEALCKKSParameters:
    """Configuration for the optional TenSEAL CKKS backend.

    The defaults are common TenSEAL demo parameters and should be adjusted for
    serious experiments according to the target security level and multiplicative
    depth.
    """

    poly_modulus_degree: int = 8192
    coeff_mod_bit_sizes: tuple[int, ...] = (60, 40, 40, 40, 60)
    global_scale_bits: int = 40
    generate_galois_keys: bool = True
    generate_relin_keys: bool = True


class TenSEALCKKSContext:
    """CKKS context adapter backed by the external TenSEAL library."""

    def __init__(self, params: TenSEALCKKSParameters | None = None) -> None:
        if ts is None:
            raise ImportError(
                "TenSEAL is not installed. Install it with `pip install tenseal` "
                "in an environment that can access Python packages."
            )
        self.params = params or TenSEALCKKSParameters()
        self.context = ts.context(
            ts.SCHEME_TYPE.CKKS,
            poly_modulus_degree=self.params.poly_modulus_degree,
            coeff_mod_bit_sizes=list(self.params.coeff_mod_bit_sizes),
        )
        self.context.global_scale = 2**self.params.global_scale_bits
        if self.params.generate_galois_keys:
            self.context.generate_galois_keys()
        if self.params.generate_relin_keys:
            self.context.generate_relin_keys()

    def encrypt(self, value: float | Sequence[float]) -> "TenSEALCiphertext":
        """Encrypt a scalar or packed vector with TenSEAL CKKSVector."""

        values = [float(value)] if isinstance(value, (int, float)) else list(value)
        return TenSEALCiphertext(self, ts.ckks_vector(self.context, values))

    def encrypt_vector(self, values: Iterable[float]) -> "TenSEALCiphertext":
        """Encrypt a packed CKKS vector."""

        return self.encrypt(list(values))

    def decrypt(self, ciphertext: "TenSEALCiphertext") -> list[float]:
        """Decrypt to a list of approximate slot values."""

        self._check_context(ciphertext)
        return [float(value) for value in ciphertext.vector.decrypt()]

    def decrypt_scalar(self, ciphertext: "TenSEALCiphertext") -> float:
        """Decrypt the first slot as a scalar."""

        return self.decrypt(ciphertext)[0]

    def plaintext(
        self, value: float | Fraction, scale: int | None = None
    ) -> "TenSEALCiphertext":
        """Return an encrypted constant vector with one slot.

        TenSEAL's high-level Python API does not expose a first-class plaintext
        object for all operations used here, so the adapter encrypts constants in
        helper circuits.  This still runs on a real CKKS backend, but users who
        need plaintext-only constants can port this adapter to lower-level SEAL.
        """

        return self.encrypt(float(value))

    def eval_polynomial(
        self, ciphertext: "TenSEALCiphertext", polynomial: OddPolynomial
    ) -> "TenSEALCiphertext":
        """Evaluate an odd polynomial on a TenSEAL ciphertext."""

        y = ciphertext * ciphertext
        acc = y.mul_plain(polynomial.coefficients[-1])
        for coeff in reversed(polynomial.coefficients[:-1]):
            acc = (y * acc).add_plain(coeff)
        return ciphertext * acc

    def sign(
        self, ciphertext: "TenSEALCiphertext", n: int = 4, d: int = 8
    ) -> "TenSEALCiphertext":
        """Evaluate the NewComp sign approximation on TenSEAL ciphertexts."""

        out = ciphertext
        polynomial = f_polynomial(n)
        for _ in range(d):
            out = self.eval_polynomial(out, polynomial)
        return out

    def compare(
        self,
        a: "TenSEALCiphertext",
        b: "TenSEALCiphertext",
        n: int = 4,
        d: int = 8,
    ) -> "TenSEALCiphertext":
        """Return encrypted approximation to 1[a>b]."""

        return self.sign(a - b, n=n, d=d).mul_plain(Fraction(1, 2)).add_plain(
            Fraction(1, 2)
        )

    def compare_plain(self, a: float, b: float, n: int = 4, d: int = 8) -> float:
        """Encrypt, compare, and decrypt two scalar values."""

        return self.decrypt_scalar(self.compare(self.encrypt(a), self.encrypt(b), n=n, d=d))

    def _check_context(self, ciphertext: "TenSEALCiphertext") -> None:
        if ciphertext.context is not self:
            raise ValueError("ciphertext belongs to a different TenSEAL context")


@dataclass(frozen=True)
class TenSEALCiphertext:
    """Small wrapper exposing the arithmetic interface used by EvalComp."""

    context: TenSEALCKKSContext
    vector: object

    def __add__(self, other: "TenSEALCiphertext") -> "TenSEALCiphertext":
        self.context._check_context(other)
        return TenSEALCiphertext(self.context, self.vector + other.vector)

    def __sub__(self, other: "TenSEALCiphertext") -> "TenSEALCiphertext":
        self.context._check_context(other)
        return TenSEALCiphertext(self.context, self.vector - other.vector)

    def __mul__(self, other: "TenSEALCiphertext") -> "TenSEALCiphertext":
        self.context._check_context(other)
        return TenSEALCiphertext(self.context, self.vector * other.vector)

    def add_plain(self, value: float | Fraction) -> "TenSEALCiphertext":
        return TenSEALCiphertext(self.context, self.vector + float(value))

    def mul_plain(self, value: float | Fraction) -> "TenSEALCiphertext":
        return TenSEALCiphertext(self.context, self.vector * float(value))

    @property
    def scale(self) -> int:
        return 1

    def decrypt(self) -> list[float]:
        return self.context.decrypt(self)

    def decrypt_scalar(self) -> float:
        return self.context.decrypt_scalar(self)


__all__ = ["TenSEALCKKSContext", "TenSEALCKKSParameters", "TenSEALCiphertext"]
