# Efficient Homomorphic Comparison Methods — CKKS Reproduction

This repository reproduces the comparison circuit from:

> Jung Hee Cheon, Dongwoo Kim, and Duhyeong Kim, **Efficient Homomorphic
> Comparison Methods with Optimal Complexity**, ASIACRYPT 2020.

The project now contains two layers:

1. **Reference arithmetic layer**: exact construction of the paper's composite
   polynomials for fast plaintext experimentation.
2. **Native CKKS/HEAAN-style ciphertext layer**: a dependency-free RLWE
   approximate-arithmetic backend over `Z_q[X]/(X^N+1)` with encryption,
   decryption, ciphertext addition/subtraction, ciphertext multiplication,
   plaintext multiplication, rescaling, scalar-slot vectors, and encrypted
   evaluation of the `NewComp` comparison circuit.

The native backend is intentionally compact and auditable so the complete
ciphertext flow can be inspected in one repository.  It is appropriate for
reproduction, teaching, and circuit debugging.  For production security, port the
same evaluator to a hardened CKKS library such as HEAAN, SEAL, OpenFHE, Lattigo,
or TenSEAL, because this repository does not attempt side-channel hardening,
SIMD canonical embedding, or formal parameter validation.

This branch also adds an EvalRound/EvalComp reproduction path: EvalComp replaces
the modular-function approximation in CKKS bootstrapping with HCF-based
homomorphic rounding, then lets the comparison-function component be swapped for
the optimized cipher-sign method.

## Implemented comparison method

For `n >= 1`, the paper defines the normalized odd polynomial

```text
f_n(x) = c_n ∫_0^x (1 - s^2)^n ds,
```

where `c_n` makes `f_n(1)=1`.  Expanding the integral gives

```text
f_n(x)=c_n Σ_{i=0}^n (-1)^i binom(n,i) x^(2i+1)/(2i+1).
```

The comparison approximation is then

```text
NewComp(a,b;n,d) = (f_n composed d times at (a-b) + 1) / 2,
```

for scaled inputs `a,b ∈ [0,1]`.  The exact comparison convention is `1` for
`a>b`, `0` for `a<b`, and `1/2` for ties.

## Quick start

```bash
python -m venv .venv
. .venv/bin/activate
pip install -e . pytest
pytest
```

Print the degree-9 polynomial `f_4`:

```bash
homcomp coeffs -n 4
```

Run a plaintext comparison:

```bash
homcomp compare 0.75 0.25 -n 4 -d 6
```

Run an encrypted CKKS comparison end to end.  The command encrypts both inputs,
evaluates the polynomial comparison on ciphertexts, and decrypts only the final
result:

```bash
homcomp ckks-demo 0.75 0.25 -n 1 -d 2 --scale-bits 20 --modulus-bits 4096
```

Run the encrypted vector demo:

```bash
PYTHONPATH=src python examples/encrypted_ckks_demo.py
```

Run the EvalComp demo with the optimized cipher-sign replacement:

```bash
homcomp evalcomp-demo 0.25 --sign-method optimized --optimized-degrees 3,5 --integer-bound 1
PYTHONPATH=src python examples/evalcomp_replacement_demo.py
```

Run a plaintext grid sweep:

```bash
homcomp sweep --epsilon 0.0625 --alpha 8 -n 4 --grid-size 65
```

Print the transcribed paper-table comparison report:

```bash
homcomp paper-report --combined-alpha 8 --hcf-share 1.0
PYTHONPATH=src python examples/paper_comparison_report.py
```

## Native CKKS API example

```python
from homcomp.ckks import CKKSParameters, NativeCKKSContext

ctx = NativeCKKSContext(
    CKKSParameters(poly_degree=8, coefficient_modulus_bits=4096, scale_bits=20),
    seed=7,
)

enc_a = ctx.encrypt(0.75)
enc_b = ctx.encrypt(0.25)
enc_cmp = ctx.compare(enc_a, enc_b, n=1, d=2)
print(enc_cmp.decrypt())
```

The `NativeCKKSContext.encrypt_vector` helper stores each scalar as an
independent CKKS ciphertext slot and exposes `CKKSEncryptedVector.compare` for
slotwise encrypted comparisons.

## EvalComp + optimized sign replacement

`homcomp.evalcomp` implements the two requested reproduction steps:

1. **EvalComp on an EvalRound-style path.**  A ciphertext encrypting
   `x = m + qI` is processed by `eval_round`, which estimates `q round(x/q)`
   from the endpoint comparison vector, interval judgement vector, midpoint
   vector, and inner product.  `bootstrap` subtracts that estimated integer
   multiple from the original ciphertext.
2. **Comparison-function replacement.**  The HCF sign component is configurable.
   `CompositePolynomialSign` is the baseline HCF sign evaluator, while
   `OptimizedCipherSign` replaces it with the optimized odd cipher-symbol
   pipeline.  The optimized path dynamically generates finite cubic/quintic
   stages, keeps the comparison polynomial wave range in the intended bounded
   interval, and jointly composes the solved cubic and quintic stages by default
   (`--optimized-degrees 3,5`).  The limiting coefficients are not hand-entered:
   `solve_optimized_sign_coefficients(3)` and
   `solve_optimized_sign_coefficients(5)` solve the cubic/quintic equation
   systems exactly, while `generate_dynamic_polynomial_stages` solves the new
   error-band equations used by the dynamic composite construction.

## Reported paper-table trend

The supplied paper tables report that EvalComp is faster and more precise than
EvalMod/EvalRound for the matched parameter families, and that NewCompH3&5
reduces both comparison complexity and HEAAN runtime versus NewCompG/H.  The
`paper-report` command prints those reported improvements and also estimates a
combined **NewEvalCompH3&5** scheme: EvalComp's HCF-based bootstrapping flow plus
the NewCompH3&5 comparison-function replacement.  The estimate keeps EvalComp's
reported precision/modulus metrics and applies the H3&5 comparison-component
speedup to a configurable HCF time share (`--hcf-share`).  The native Python CKKS
backend remains an educational scalar-slot backend, so its wall-clock time is not
directly comparable to the C++/HEAAN laboratory timings in the papers.

## Scope

- Implemented: exact `f_n` generation, plaintext `NewComp`, native CKKS-style
  ciphertexts, encrypted `NewComp`, scalar-slot encrypted vectors, EvalComp/HCF
  rounding, optimized cipher-sign replacement, CLI demos, and tests.
- Not implemented: `NewCompG` minimax/Remez acceleration, full RNS modulus
  chains, packed SIMD bootstrapping transforms, and production security audits.
