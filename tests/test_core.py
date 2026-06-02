from fractions import Fraction

import pytest

from homcomp import compare, compare_integers, compose_sign, f_polynomial
from homcomp.core import comparison_error, separated_pairs


def test_f1_matches_newton_sign_polynomial():
    poly = f_polynomial(1)
    assert poly.coefficients == (Fraction(3, 2), Fraction(-1, 2))
    assert poly(Fraction(1, 1)) == 1
    assert poly(Fraction(-1, 1)) == -1


def test_f4_is_odd_and_endpoint_normalized():
    poly = f_polynomial(4)
    for x in [Fraction(0), Fraction(1, 4), Fraction(3, 4), Fraction(1)]:
        assert poly(-x) == -poly(x)
    assert poly(Fraction(1, 1)) == 1


def test_composition_moves_away_from_zero_toward_sign():
    assert compose_sign(0.2, n=4, d=2) > 0.2
    assert compose_sign(-0.2, n=4, d=2) < -0.2
    assert compose_sign(0.75, n=4, d=4) > 0.999
    assert compose_sign(-0.75, n=4, d=4) < -0.999


def test_compare_convention_and_accuracy_on_separated_grid():
    assert compare(0.5, 0.5, n=4, d=6) == pytest.approx(0.5)
    assert compare(0.75, 0.25, n=4, d=6) > 0.999
    assert compare(0.25, 0.75, n=4, d=6) < 0.001
    err = comparison_error(separated_pairs(0.25, grid_size=9), n=4, d=6)
    assert err < 1e-3


def test_unsigned_integer_helper():
    assert compare_integers(7, 3, bits=3, n=4, d=7) > 0.99
    assert compare_integers(3, 7, bits=3, n=4, d=7) < 0.01
    with pytest.raises(ValueError):
        compare_integers(8, 0, bits=3)
