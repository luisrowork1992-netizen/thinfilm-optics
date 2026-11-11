"""Simple dispersion models for common optical materials.

Wavelengths are in micrometres unless stated otherwise.
"""

from __future__ import annotations

import numpy as np

__all__ = ["sellmeier", "cauchy", "fused_silica", "bk7"]


def sellmeier(lam_um, B, C):
    """Sellmeier equation ``n^2 = 1 + sum B_i lam^2 / (lam^2 - C_i)``.

    ``C`` coefficients are in um^2.
    """
    lam2 = np.asarray(lam_um, dtype=float) ** 2
    n2 = 1.0 + sum(b * lam2 / (lam2 - c) for b, c in zip(B, C))
    return np.sqrt(n2)


def cauchy(lam_um, A, B=0.0, C=0.0):
    """Cauchy model ``n = A + B / lam^2 + C / lam^4``."""
    lam = np.asarray(lam_um, dtype=float)
    return A + B / lam**2 + C / lam**4


# Malitson (1965) coefficients for fused silica, valid ~0.21-3.71 um.
_SILICA_B = (0.6961663, 0.4079426, 0.8974794)
_SILICA_C = (0.0684043**2, 0.1162414**2, 9.896161**2)

# Schott N-BK7 catalogue coefficients, valid ~0.3-2.5 um.
_BK7_B = (1.03961212, 0.231792344, 1.01046945)
_BK7_C = (0.00600069867, 0.0200179144, 103.560653)


def fused_silica(lam_um):
    """Refractive index of fused silica (SiO2)."""
    return sellmeier(lam_um, _SILICA_B, _SILICA_C)


def bk7(lam_um):
    """Refractive index of Schott N-BK7 glass."""
    return sellmeier(lam_um, _BK7_B, _BK7_C)
