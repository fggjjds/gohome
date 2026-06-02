"""Small plaintext experiment for the NewComp approximation."""

from homcomp import comparison_error, estimate_newcomp_iterations
from homcomp.core import separated_pairs


for bits in (4, 8, 12):
    epsilon = 1 / ((1 << bits) - 1)
    for n in (2, 4, 8):
        d = estimate_newcomp_iterations(epsilon=epsilon, alpha=bits, n=n)
        err = comparison_error(separated_pairs(epsilon, grid_size=65), n=n, d=d)
        print(f"bits={bits:2d} n={n:2d} d={d:2d} max_error={err:.3e}")
