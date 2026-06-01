"""Command line interface for the reproduction package."""

from __future__ import annotations

import argparse

from .ckks import CKKSParameters, NativeCKKSContext
from .evalcomp import EvalCompBootstrapper, EvalCompParameters
from .paper_results import format_paper_report
from .core import (
    compare,
    compare_integers,
    comparison_error,
    estimate_newcomp_iterations,
    f_polynomial,
    separated_pairs,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Reproduce composite-polynomial homomorphic comparison experiments."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    coeffs = sub.add_parser("coeffs", help="print f_n coefficients")
    coeffs.add_argument("-n", type=int, default=4)

    cmp_parser = sub.add_parser("compare", help="compare two scaled values in [0,1]")
    cmp_parser.add_argument("a", type=float)
    cmp_parser.add_argument("b", type=float)
    cmp_parser.add_argument("-n", type=int, default=4)
    cmp_parser.add_argument("-d", type=int, default=8)

    int_parser = sub.add_parser("compare-int", help="compare unsigned integers")
    int_parser.add_argument("a", type=int)
    int_parser.add_argument("b", type=int)
    int_parser.add_argument("--bits", type=int, required=True)
    int_parser.add_argument("-n", type=int, default=4)
    int_parser.add_argument("-d", type=int)

    bench = sub.add_parser("sweep", help="print max error on a separated grid")
    bench.add_argument("--epsilon", type=float, default=1 / 16)
    bench.add_argument("--alpha", type=int, default=8)
    bench.add_argument("--grid-size", type=int, default=65)
    bench.add_argument("-n", type=int, default=4)
    bench.add_argument("-d", type=int)

    ckks = sub.add_parser(
        "ckks-demo", help="encrypt, compare, and decrypt two scaled values"
    )
    ckks.add_argument("a", type=float)
    ckks.add_argument("b", type=float)
    ckks.add_argument("-n", type=int, default=1)
    ckks.add_argument("-d", type=int, default=2)
    ckks.add_argument("--poly-degree", type=int, default=8)
    ckks.add_argument("--modulus-bits", type=int, default=4096)
    ckks.add_argument("--scale-bits", type=int, default=20)
    ckks.add_argument("--seed", type=int, default=0)

    evalcomp = sub.add_parser(
        "evalcomp-demo",
        help="run EvalComp-style encrypted rounding/bootstrap with HCF",
    )
    evalcomp.add_argument("value", type=float)
    evalcomp.add_argument("--modulus-step", type=float, default=1.0)
    evalcomp.add_argument("--integer-bound", type=int, default=1)
    evalcomp.add_argument("--sign-method", choices=["baseline", "optimized"], default="optimized")
    evalcomp.add_argument("--sign-iterations", type=int, default=1)
    evalcomp.add_argument("--sign-n", type=int, default=1)
    evalcomp.add_argument(
        "--optimized-degrees",
        default="3,5",
        help="comma-separated optimized sign stages, e.g. 3,5",
    )
    evalcomp.add_argument("--optimized-alpha", type=int, default=8)
    evalcomp.add_argument("--hcf-domain-radius", type=float)
    evalcomp.add_argument("--poly-degree", type=int, default=8)
    evalcomp.add_argument("--modulus-bits", type=int, default=4096)
    evalcomp.add_argument("--scale-bits", type=int, default=16)
    evalcomp.add_argument("--seed", type=int, default=0)

    report = sub.add_parser(
        "paper-report",
        help="print reported paper-table improvements for EvalComp and H3&5",
    )
    report.add_argument("--combined-alpha", type=int, default=8)
    report.add_argument("--hcf-share", type=float, default=1.0)
    report.add_argument("--baseline-comparison", choices=["G", "H"], default="H")

    args = parser.parse_args(argv)

    if args.cmd == "coeffs":
        poly = f_polynomial(args.n)
        print(f"degree: {poly.degree}")
        print(poly.format())
        return 0

    if args.cmd == "compare":
        print(compare(args.a, args.b, n=args.n, d=args.d))
        return 0

    if args.cmd == "compare-int":
        print(compare_integers(args.a, args.b, bits=args.bits, n=args.n, d=args.d))
        return 0

    if args.cmd == "sweep":
        d = (
            args.d
            if args.d is not None
            else estimate_newcomp_iterations(args.epsilon, args.alpha, args.n)
        )
        err = comparison_error(separated_pairs(args.epsilon, args.grid_size), n=args.n, d=d)
        print(
            f"n={args.n} d={d} epsilon={args.epsilon:g} "
            f"grid={args.grid_size} max_error={err:.6g}"
        )
        return 0

    if args.cmd == "ckks-demo":
        ctx = NativeCKKSContext(
            CKKSParameters(
                poly_degree=args.poly_degree,
                coefficient_modulus_bits=args.modulus_bits,
                scale_bits=args.scale_bits,
            ),
            seed=args.seed,
        )
        result = ctx.compare_plain(args.a, args.b, n=args.n, d=args.d)
        print(result)
        return 0

    if args.cmd == "paper-report":
        print(
            format_paper_report(
                combined_alpha=args.combined_alpha,
                hcf_share=args.hcf_share,
                baseline_comparison=args.baseline_comparison,
            )
        )
        return 0

    if args.cmd == "evalcomp-demo":
        ctx = NativeCKKSContext(
            CKKSParameters(
                poly_degree=args.poly_degree,
                coefficient_modulus_bits=args.modulus_bits,
                scale_bits=args.scale_bits,
            ),
            seed=args.seed,
        )
        params = EvalCompParameters(
            modulus_step=args.modulus_step,
            integer_bound=args.integer_bound,
            sign_iterations=args.sign_iterations,
            sign_n=args.sign_n,
            optimized_sign=args.sign_method == "optimized",
            optimized_sign_degrees=_parse_optimized_degrees(args.optimized_degrees),
            optimized_precision_alpha=args.optimized_alpha,
            hcf_domain_radius=args.hcf_domain_radius,
        )
        bootstrapper = EvalCompBootstrapper(ctx, params)
        ciphertext = ctx.encrypt(args.value)
        rounded = bootstrapper.eval_round(ciphertext).decrypt()
        bootstrapped = bootstrapper.bootstrap(ciphertext).decrypt()
        print(f"eval_round={rounded}")
        print(f"evalcomp={bootstrapped}")
        return 0

    raise AssertionError("unreachable")


def _parse_optimized_degrees(value: str) -> tuple[int, ...]:
    degrees = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    if not degrees:
        raise ValueError("--optimized-degrees must contain at least one degree")
    invalid = [degree for degree in degrees if degree not in {3, 5}]
    if invalid:
        raise ValueError(f"--optimized-degrees supports only 3 and 5: {invalid}")
    return degrees


if __name__ == "__main__":
    raise SystemExit(main())
