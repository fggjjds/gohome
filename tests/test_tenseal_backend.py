import pytest

from homcomp.tenseal_backend import TenSEALCKKSContext, TenSEALCKKSParameters, ts


def test_tenseal_backend_reports_missing_dependency_or_runs_compare():
    if ts is None:
        with pytest.raises(ImportError):
            TenSEALCKKSContext(TenSEALCKKSParameters())
        return

    ctx = TenSEALCKKSContext(
        TenSEALCKKSParameters(
            poly_modulus_degree=8192,
            coeff_mod_bit_sizes=(60, 40, 40, 60),
            global_scale_bits=40,
        )
    )
    assert ctx.compare_plain(0.75, 0.25, n=1, d=1) > 0.5
