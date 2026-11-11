import math

import numpy as np
import pytest

from thinfilm import (
    bk7,
    bragg_mirror,
    coh_tmm,
    fused_silica,
    quarter_wave_stack_reflectance,
    single_layer_ar,
)

INF = math.inf


def test_bare_glass_normal_incidence_is_four_percent():
    res = coh_tmm("s", [1.0, 1.5], [INF, INF], 0.0, 0.55)
    assert res.R == pytest.approx(0.04, abs=1e-12)
    assert res.T == pytest.approx(0.96, abs=1e-12)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_zero_thickness_layer_has_no_effect(pol):
    bare = coh_tmm(pol, [1.0, 1.5], [INF, INF], 0.3, 0.6)
    coated = coh_tmm(pol, [1.0, 2.1, 1.5], [INF, 0.0, INF], 0.3, 0.6)
    assert coated.R == pytest.approx(bare.R, abs=1e-12)


def test_ideal_quarter_wave_ar_coating_is_reflectionless():
    n_list, d_list = single_layer_ar(1.0, 1.52, lam0=0.55)
    res = coh_tmm("s", n_list, d_list, 0.0, 0.55)
    assert res.R == pytest.approx(0.0, abs=1e-12)


def test_half_wave_layer_is_absentee():
    n_layer, lam = 1.38, 0.55
    d = lam / (2 * n_layer)
    bare = coh_tmm("s", [1.0, 1.52], [INF, INF], 0.0, lam)
    half = coh_tmm("s", [1.0, n_layer, 1.52], [INF, d, INF], 0.0, lam)
    assert half.R == pytest.approx(bare.R, abs=1e-12)


def test_brewster_angle_p_reflectance_vanishes():
    theta_b = math.atan(1.5)
    res = coh_tmm("p", [1.0, 1.5], [INF, INF], theta_b, 0.6)
    assert res.R == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("theta", [0.0, 0.4, 0.9, 1.3])
def test_energy_conservation_lossless(pol, theta):
    n_list = [1.0, 2.35, 1.46, 2.35, 1.52]
    d_list = [INF, 0.06, 0.1, 0.07, INF]
    res = coh_tmm(pol, n_list, d_list, theta, 0.633)
    assert res.R + res.T == pytest.approx(1.0, abs=1e-10)


def test_absorbing_layer_absorbs():
    res = coh_tmm("s", [1.0, 0.2 + 3.0j, 1.5], [INF, 0.02, INF], 0.0, 0.55)
    assert 0 < res.A < 1
    assert res.R + res.T < 1


def test_total_internal_reflection():
    # glass -> air beyond the critical angle
    theta = math.radians(60)
    for pol in ("s", "p"):
        res = coh_tmm(pol, [1.5, 1.0], [INF, INF], theta, 0.55)
        assert res.R == pytest.approx(1.0, abs=1e-12)
        assert res.T == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("pairs", [1, 3, 8])
def test_bragg_mirror_matches_closed_form(pairs):
    n_list, d_list = bragg_mirror(1.0, 2.35, 1.46, 1.52, pairs, lam0=1.064)
    res = coh_tmm("s", n_list, d_list, 0.0, 1.064)
    expected = quarter_wave_stack_reflectance(1.0, 2.35, 1.46, 1.52, pairs)
    assert res.R == pytest.approx(expected, rel=1e-10)


def test_s_and_p_agree_at_normal_incidence():
    n_list = [1.0, 2.0, 1.4, 1.6]
    d_list = [INF, 0.1, 0.2, INF]
    s = coh_tmm("s", n_list, d_list, 0.0, 0.5)
    p = coh_tmm("p", n_list, d_list, 0.0, 0.5)
    assert s.R == pytest.approx(p.R, abs=1e-12)


def test_input_validation():
    with pytest.raises(ValueError):
        coh_tmm("x", [1.0, 1.5], [INF, INF], 0.0, 0.5)
    with pytest.raises(ValueError):
        coh_tmm("s", [1.0, 1.5], [0.0, INF], 0.0, 0.5)
    with pytest.raises(ValueError):
        coh_tmm("s", [1.0, 1.5, 1.4], [INF, INF], 0.0, 0.5)


def test_sellmeier_reference_values():
    # Sodium d-line 587.6 nm reference values
    assert fused_silica(0.5876) == pytest.approx(1.4585, abs=2e-4)
    assert bk7(0.5876) == pytest.approx(1.5168, abs=2e-4)
    # Normal dispersion: index decreases with wavelength in the visible
    wl = np.linspace(0.4, 0.8, 20)
    assert np.all(np.diff(fused_silica(wl)) < 0)
