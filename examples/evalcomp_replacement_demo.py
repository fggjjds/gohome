"""Compare EvalComp HCF with the optimized cipher-sign replacement."""

from homcomp.ckks import CKKSParameters, NativeCKKSContext
from homcomp.evalcomp import EvalCompBootstrapper, EvalCompParameters
from homcomp.sign_methods import (
    generate_dynamic_polynomial_stages,
    solve_optimized_sign_coefficients,
)


ctx = NativeCKKSContext(
    CKKSParameters(poly_degree=8, coefficient_modulus_bits=4096, scale_bits=16),
    seed=2026,
)

print("cubic coefficients:", solve_optimized_sign_coefficients(3))
print("quintic coefficients:", solve_optimized_sign_coefficients(5))
print(
    "dynamic stages:",
    [
        (stage.degree, round(stage.target_error, 8))
        for stage in generate_dynamic_polynomial_stages(0.5, 8)
    ],
)

for optimized in (False, True):
    params = EvalCompParameters(
        integer_bound=1,
        sign_iterations=1,
        optimized_sign=optimized,
        hcf_domain_radius=1.0,
    )
    bootstrapper = EvalCompBootstrapper(ctx, params)
    high = bootstrapper.hcf(ctx.encrypt(0.75), 0.5).decrypt()
    low = bootstrapper.hcf(ctx.encrypt(0.25), 0.5).decrypt()
    rounded = bootstrapper.eval_round(ctx.encrypt(0.25)).decrypt()
    print(
        f"optimized={optimized} hcf_high={high:.6f} "
        f"hcf_low={low:.6f} eval_round={rounded:.6f}"
    )
