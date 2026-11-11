"""Helpers for building classic thin-film designs."""

from __future__ import annotations

import math

__all__ = [
    "quarter_wave_thickness",
    "single_layer_ar",
    "bragg_mirror",
    "quarter_wave_stack_reflectance",
]


def quarter_wave_thickness(n_layer: float, lam0: float) -> float:
    """Physical thickness of a quarter-wave layer at design wavelength ``lam0``."""
    if n_layer <= 0 or lam0 <= 0:
        raise ValueError("index and wavelength must be positive")
    return lam0 / (4.0 * n_layer)


def single_layer_ar(n0: float, n_sub: float, lam0: float, n_coat: float | None = None):
    """Single quarter-wave anti-reflection coating.

    If ``n_coat`` is omitted the ideal index ``sqrt(n0 * n_sub)`` is used.

    Returns ``(n_list, d_list)`` ready for :func:`thinfilm.coh_tmm`.
    """
    if n_coat is None:
        n_coat = math.sqrt(n0 * n_sub)
    d = quarter_wave_thickness(n_coat, lam0)
    return [n0, n_coat, n_sub], [math.inf, d, math.inf]


def bragg_mirror(n0: float, n_high: float, n_low: float, n_sub: float, pairs: int, lam0: float):
    """Quarter-wave (HL)^N stack, high-index layer facing the incident medium.

    Returns ``(n_list, d_list)``.
    """
    if pairs < 1:
        raise ValueError("pairs must be >= 1")
    dh = quarter_wave_thickness(n_high, lam0)
    dl = quarter_wave_thickness(n_low, lam0)
    n_list = [n0] + [n_high, n_low] * pairs + [n_sub]
    d_list = [math.inf] + [dh, dl] * pairs + [math.inf]
    return n_list, d_list


def quarter_wave_stack_reflectance(n0: float, n_high: float, n_low: float, n_sub: float, pairs: int) -> float:
    """Closed-form peak reflectance of an (HL)^N quarter-wave stack at normal incidence.

    Each quarter-wave layer maps the optical admittance ``Y -> n_layer**2 / Y``,
    so after N pairs on substrate ``n_sub`` the admittance is
    ``Y = (n_high / n_low)**(2N) * n_sub``.
    """
    y = (n_high / n_low) ** (2 * pairs) * n_sub
    return ((n0 - y) / (n0 + y)) ** 2
