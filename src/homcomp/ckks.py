"""Native CKKS/HEAAN-style ciphertext backend.

The implementation in this module is intentionally small and auditable: it
implements a scalar-slot, leveled RLWE approximate-arithmetic scheme over
``Z_q[X]/(X^N+1)`` with CKKS-style fixed-point scaling.  It is suitable for
reproducing the encrypted comparison circuit end-to-end in tests and examples.
For production deployments you should port the same circuit to a hardened CKKS
library (HEAAN, SEAL, OpenFHE, Lattigo, TenSEAL, ...), because this educational
backend does not implement side-channel hardening, SIMD canonical embedding, or
security-parameter validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from random import Random
from typing import Iterable, Sequence

from .core import OddPolynomial, f_polynomial


Polynomial = tuple[int, ...]


def _center(value: int, modulus: int) -> int:
    """Return the centered representative of value modulo modulus."""

    value %= modulus
    if value > modulus // 2:
        value -= modulus
    return value


def _trim(poly: Sequence[int], degree: int) -> Polynomial:
    """Normalize a polynomial to exactly ``degree`` coefficients."""

    return tuple(int(poly[i]) if i < len(poly) else 0 for i in range(degree))


def _poly_add(a: Polynomial, b: Polynomial, modulus: int) -> Polynomial:
    return tuple((x + y) % modulus for x, y in zip(a, b))


def _poly_sub(a: Polynomial, b: Polynomial, modulus: int) -> Polynomial:
    return tuple((x - y) % modulus for x, y in zip(a, b))


def _poly_scalar_mul(a: Polynomial, scalar: int, modulus: int) -> Polynomial:
    return tuple((x * scalar) % modulus for x in a)


def _poly_mul(a: Polynomial, b: Polynomial, modulus: int) -> Polynomial:
    """Negacyclic multiplication in Z_q[X]/(X^N+1)."""

    degree = len(a)
    out = [0] * degree
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            k = i + j
            if k >= degree:
                out[k - degree] -= ai * bj
            else:
                out[k] += ai * bj
    return tuple(x % modulus for x in out)


def _round_fraction(numerator: int, denominator: int) -> int:
    """Round numerator/denominator to the nearest integer without floats."""

    if denominator <= 0:
        raise ValueError("denominator must be positive")
    sign = -1 if numerator < 0 else 1
    numerator = abs(numerator)
    quotient, remainder = divmod(numerator, denominator)
    if 2 * remainder >= denominator:
        quotient += 1
    return sign * quotient

def _poly_round_div(a: Polynomial, divisor: int, modulus: int) -> Polynomial:
    """Divide centered coefficients by ``divisor`` and round."""

    return tuple(_round_fraction(_center(x, modulus), divisor) % modulus for x in a)


def _zero_poly(degree: int) -> Polynomial:
    return (0,) * degree


@dataclass(frozen=True)
class CKKSParameters:
    """Parameters for the native scalar-slot CKKS backend.

    ``poly_degree`` must be a power of two.  ``coefficient_modulus_bits`` should
    be large enough for the chosen comparison depth; the defaults prioritize
    reproducible correctness for small paper-style experiments over speed.
    """

    poly_degree: int = 16
    coefficient_modulus_bits: int = 8192
    scale_bits: int = 30
    error_bound: int = 1

    def __post_init__(self) -> None:
        if self.poly_degree < 2 or self.poly_degree & (self.poly_degree - 1):
            raise ValueError("poly_degree must be a power of two and at least 2")
        if self.coefficient_modulus_bits <= self.scale_bits * 4:
            raise ValueError("coefficient modulus is too small for CKKS scaling")
        if self.scale_bits < 10:
            raise ValueError("scale_bits must be at least 10")
        if self.error_bound < 0:
            raise ValueError("error_bound cannot be negative")

    @property
    def modulus(self) -> int:
        return 1 << self.coefficient_modulus_bits

    @property
    def scale(self) -> int:
        return 1 << self.scale_bits


class NativeCKKSContext:
    """Secret-key CKKS context with encrypt/decrypt/evaluate operations."""

    def __init__(self, params: CKKSParameters | None = None, seed: int | None = 0):
        self.params = params or CKKSParameters()
        self._rng = Random(seed)
        self.secret_key = self._sample_secret()

    @property
    def modulus(self) -> int:
        return self.params.modulus

    @property
    def degree(self) -> int:
        return self.params.poly_degree

    @property
    def scale(self) -> int:
        return self.params.scale

    def encrypt(self, value: float, scale: int | None = None) -> "CKKSCiphertext":
        """Encrypt one approximate scalar value."""

        scale = scale or self.scale
        plaintext = self._encode(value, scale)
        a = self._sample_uniform_poly()
        error = self._sample_error_poly()
        c0 = _poly_add(_poly_mul(a, self.secret_key, self.modulus), error, self.modulus)
        c0 = _poly_add(c0, plaintext, self.modulus)
        c1 = _poly_scalar_mul(a, -1, self.modulus)
        return CKKSCiphertext(self, (c0, c1), scale)

    def encrypt_vector(
        self, values: Iterable[float], scale: int | None = None
    ) -> "CKKSEncryptedVector":
        """Encrypt a vector as independent scalar-slot ciphertexts."""

        return CKKSEncryptedVector([self.encrypt(value, scale=scale) for value in values])

    def decrypt(self, ciphertext: "CKKSCiphertext") -> float:
        """Decrypt one scalar ciphertext to a CKKS approximate value."""

        if ciphertext.context is not self:
            raise ValueError("ciphertext belongs to a different context")
        accum = _zero_poly(self.degree)
        secret_power = _trim((1,), self.degree)
        for component in ciphertext.components:
            accum = _poly_add(
                accum, _poly_mul(component, secret_power, self.modulus), self.modulus
            )
            secret_power = _poly_mul(secret_power, self.secret_key, self.modulus)
        return _center(accum[0], self.modulus) / ciphertext.scale

    def decrypt_vector(self, ciphertext: "CKKSEncryptedVector") -> list[float]:
        """Decrypt a vector produced by :meth:`encrypt_vector`."""

        return [self.decrypt(slot) for slot in ciphertext.slots]

    def plaintext(self, value: float, scale: int | None = None) -> "CKKSCiphertext":
        """Create a plaintext constant in ciphertext form for circuit building."""

        scale = scale or self.scale
        return CKKSCiphertext(self, (self._encode(value, scale),), scale)

    def eval_polynomial(
        self, ciphertext: "CKKSCiphertext", polynomial: OddPolynomial
    ) -> "CKKSCiphertext":
        """Homomorphically evaluate an odd polynomial on a ciphertext."""

        y = ciphertext * ciphertext
        acc = self.plaintext(polynomial.coefficients[-1], scale=ciphertext.scale)
        for coeff in reversed(polynomial.coefficients[:-1]):
            acc = (y * acc).add_plain(coeff)
        return ciphertext * acc

    def sign(self, ciphertext: "CKKSCiphertext", n: int = 4, d: int = 8) -> "CKKSCiphertext":
        """Evaluate the encrypted NewComp sign approximation."""

        if d < 0:
            raise ValueError("d must be non-negative")
        polynomial = f_polynomial(n)
        out = ciphertext
        for _ in range(d):
            out = self.eval_polynomial(out, polynomial)
        return out

    def compare(
        self, a: "CKKSCiphertext", b: "CKKSCiphertext", n: int = 4, d: int = 8
    ) -> "CKKSCiphertext":
        """Return an encrypted approximation of 1[a>b]."""

        return self.sign(a - b, n=n, d=d).mul_plain(Fraction(1, 2)).add_plain(
            Fraction(1, 2)
        )

    def compare_plain(self, a: float, b: float, n: int = 4, d: int = 8) -> float:
        """Encrypt, compare, and decrypt two scaled values for demos/tests."""

        return self.decrypt(self.compare(self.encrypt(a), self.encrypt(b), n=n, d=d))

    def _encode(self, value: float | Fraction, scale: int) -> Polynomial:
        if isinstance(value, Fraction):
            coeff = _round_fraction(value.numerator * scale, value.denominator)
        else:
            coeff = round(value * scale)
        return _trim((coeff % self.modulus,), self.degree)

    def _sample_secret(self) -> Polynomial:
        return tuple(
            self._rng.choice((-1, 0, 1)) % self.modulus
            for _ in range(self.degree)
        )

    def _sample_uniform_poly(self) -> Polynomial:
        return tuple(self._rng.randrange(self.modulus) for _ in range(self.degree))

    def _sample_error_poly(self) -> Polynomial:
        bound = self.params.error_bound
        if bound == 0:
            return _zero_poly(self.degree)
        return tuple(
            self._rng.randint(-bound, bound) % self.modulus
            for _ in range(self.degree)
        )


@dataclass(frozen=True)
class CKKSCiphertext:
    """Ciphertext with one or more RLWE components and a CKKS scale."""

    context: NativeCKKSContext
    components: tuple[Polynomial, ...]
    scale: int

    def __post_init__(self) -> None:
        if self.scale <= 0:
            raise ValueError("scale must be positive")
        normalized = tuple(
            _trim(component, self.context.degree) for component in self.components
        )
        object.__setattr__(self, "components", normalized)

    def __add__(self, other: "CKKSCiphertext") -> "CKKSCiphertext":
        left_ct, right_ct, scale = self._align_pair(other)
        left, right = self._pad_components(left_ct.components, right_ct.components)
        components = tuple(
            _poly_add(a, b, self.context.modulus) for a, b in zip(left, right)
        )
        return CKKSCiphertext(self.context, components, scale)

    def __sub__(self, other: "CKKSCiphertext") -> "CKKSCiphertext":
        left_ct, right_ct, scale = self._align_pair(other)
        left, right = self._pad_components(left_ct.components, right_ct.components)
        components = tuple(
            _poly_sub(a, b, self.context.modulus) for a, b in zip(left, right)
        )
        return CKKSCiphertext(self.context, components, scale)

    def __mul__(self, other: "CKKSCiphertext") -> "CKKSCiphertext":
        if self.context is not other.context:
            raise ValueError("ciphertexts belong to different contexts")
        out = [_zero_poly(self.context.degree)] * (
            len(self.components) + len(other.components) - 1
        )
        for i, left in enumerate(self.components):
            for j, right in enumerate(other.components):
                product = _poly_mul(left, right, self.context.modulus)
                out[i + j] = _poly_add(out[i + j], product, self.context.modulus)
        return CKKSCiphertext(self.context, tuple(out), self.scale * other.scale)

    def add_plain(self, value: float | Fraction) -> "CKKSCiphertext":
        """Add a plaintext constant encoded at the ciphertext scale."""

        plain = self.context._encode(value, self.scale)
        components = (
            _poly_add(self.components[0], plain, self.context.modulus),
        ) + self.components[1:]
        return CKKSCiphertext(self.context, components, self.scale)

    def mul_plain(
        self, value: float | Fraction, plain_scale: int | None = None
    ) -> "CKKSCiphertext":
        """Multiply by a CKKS-encoded plaintext scalar.

        The plaintext is encoded as ``round(value * plain_scale)`` and the
        ciphertext scale is multiplied by ``plain_scale``, matching CKKS
        plaintext multiplication semantics and preserving RLWE mask cancellation.
        """

        plain_scale = plain_scale or self.context.scale
        if isinstance(value, Fraction):
            encoded = _round_fraction(value.numerator * plain_scale, value.denominator)
        else:
            encoded = round(value * plain_scale)
        components = tuple(
            _poly_scalar_mul(component, encoded, self.context.modulus)
            for component in self.components
        )
        return CKKSCiphertext(self.context, components, self.scale * plain_scale)

    def rescale(self, divisor: int | None = None) -> "CKKSCiphertext":
        """CKKS-style rescale by coefficient rounding and scale reduction."""

        divisor = divisor or self.context.scale
        if divisor <= 1:
            raise ValueError("rescale divisor must be greater than one")
        components = tuple(
            _poly_round_div(component, divisor, self.context.modulus)
            for component in self.components
        )
        return CKKSCiphertext(
            self.context, components, max(1, round(self.scale / divisor))
        )

    def decrypt(self) -> float:
        """Convenience wrapper around the owning context's decrypt method."""

        return self.context.decrypt(self)

    def _align_pair(
        self, other: "CKKSCiphertext"
    ) -> tuple["CKKSCiphertext", "CKKSCiphertext", int]:
        if self.context is not other.context:
            raise ValueError("ciphertexts belong to different contexts")
        target_scale = max(self.scale, other.scale)
        return (
            self._raise_scale(target_scale),
            other._raise_scale(target_scale),
            target_scale,
        )

    def _raise_scale(self, target_scale: int) -> "CKKSCiphertext":
        if self.scale == target_scale:
            return self
        if target_scale % self.scale != 0:
            raise ValueError("target scale must be an integer multiple of current scale")
        factor = target_scale // self.scale
        return CKKSCiphertext(
            self.context,
            tuple(
                _poly_scalar_mul(component, factor, self.context.modulus)
                for component in self.components
            ),
            target_scale,
        )

    @staticmethod
    def _pad_components(
        left: tuple[Polynomial, ...], right: tuple[Polynomial, ...]
    ) -> tuple[tuple[Polynomial, ...], tuple[Polynomial, ...]]:
        degree = len(left[0])
        size = max(len(left), len(right))
        zero = _zero_poly(degree)
        return (
            left + (zero,) * (size - len(left)),
            right + (zero,) * (size - len(right)),
        )


@dataclass(frozen=True)
class CKKSEncryptedVector:
    """Vector convenience wrapper using one native CKKS ciphertext per slot."""

    slots: list[CKKSCiphertext]

    def __post_init__(self) -> None:
        if not self.slots:
            raise ValueError("encrypted vector must contain at least one slot")
        first = self.slots[0].context
        if any(slot.context is not first for slot in self.slots):
            raise ValueError("all vector slots must belong to the same context")

    def decrypt(self) -> list[float]:
        return [slot.decrypt() for slot in self.slots]

    def compare(
        self, other: "CKKSEncryptedVector", n: int = 4, d: int = 8
    ) -> "CKKSEncryptedVector":
        if len(self.slots) != len(other.slots):
            raise ValueError("vectors must have the same length")
        context = self.slots[0].context
        return CKKSEncryptedVector(
            [
                context.compare(a, b, n=n, d=d)
                for a, b in zip(self.slots, other.slots)
            ]
        )
