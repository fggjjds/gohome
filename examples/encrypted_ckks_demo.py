"""End-to-end encrypted comparison demo with the native CKKS backend."""

from homcomp.ckks import CKKSParameters, NativeCKKSContext


ctx = NativeCKKSContext(
    CKKSParameters(poly_degree=16, coefficient_modulus_bits=8192, scale_bits=20),
    seed=2026,
)

left = ctx.encrypt_vector([0.75, 0.25, 0.50])
right = ctx.encrypt_vector([0.25, 0.75, 0.50])
result = left.compare(right, n=1, d=2)

print([round(value, 6) for value in result.decrypt()])
