"""Coherent transfer-matrix method (TMM) for planar multilayer stacks.

Convention
----------
* Complex refractive index ``n = n' + i k`` with ``k >= 0`` for absorbing media
  (time dependence ``exp(-i w t)``).
* ``n_list`` runs from the incident medium (index 0) to the exit medium
  (index -1). Both outer media are semi-infinite.
* ``d_list`` holds the physical thicknesses in the same length unit as the
  wavelength; the first and last entries must be ``inf``.
* ``theta0`` is the angle of incidence in radians, measured in the incident
  medium, which must be lossless.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["TMMResult", "coh_tmm", "spectrum", "fresnel"]


@dataclass(frozen=True)
class TMMResult:
    """Amplitude and power coefficients for one wavelength / angle."""

    r: complex
    t: complex
    R: float
    T: float

    @property
    def A(self) -> float:
        """Absorbed fraction ``1 - R - T``."""
        return 1.0 - self.R - self.T


def _forward_cos(n: complex, n0: complex, theta0: float) -> complex:
    """cos(theta) in a layer of index ``n`` via Snell's law.

    Picks the branch of the square root for a forward-travelling wave: the
    z-component of the wavevector must decay (Im > 0) or, if purely real,
    point forward (Re > 0).
    """
    sin_t = n0 * np.sin(theta0) / n
    cos_t = np.sqrt(1.0 - sin_t * sin_t + 0j)
    nc = n * cos_t
    if nc.imag < -1e-14 or (abs(nc.imag) <= 1e-14 and nc.real < 0):
        cos_t = -cos_t
    return complex(cos_t)


def fresnel(pol: str, n_i: complex, n_f: complex, cos_i: complex, cos_f: complex):
    """Fresnel amplitude coefficients (r, t) at a single interface."""
    if pol == "s":
        denom = n_i * cos_i + n_f * cos_f
        return (n_i * cos_i - n_f * cos_f) / denom, 2 * n_i * cos_i / denom
    if pol == "p":
        denom = n_f * cos_i + n_i * cos_f
        return (n_f * cos_i - n_i * cos_f) / denom, 2 * n_i * cos_i / denom
    raise ValueError("pol must be 's' or 'p'")


def coh_tmm(pol: str, n_list, d_list, theta0: float, lam: float) -> TMMResult:
    """Reflection and transmission of a coherent multilayer stack.

    Parameters
    ----------
    pol : {'s', 'p'}
        Polarisation (TE = 's', TM = 'p').
    n_list : sequence of complex
        Refractive indices, incident medium first, exit medium last.
    d_list : sequence of float
        Thicknesses; ``inf`` for the two outer media.
    theta0 : float
        Angle of incidence in radians.
    lam : float
        Vacuum wavelength, same unit as ``d_list``.
    """
    n = np.asarray(n_list, dtype=complex)
    d = np.asarray(d_list, dtype=float)
    if n.shape != d.shape or n.size < 2:
        raise ValueError("n_list and d_list must have the same length (>= 2)")
    if not (np.isinf(d[0]) and np.isinf(d[-1])):
        raise ValueError("first and last thickness must be inf")
    if abs(n[0].imag) > 1e-12:
        raise ValueError("incident medium must be lossless")
    if pol not in ("s", "p"):
        raise ValueError("pol must be 's' or 'p'")

    cos = np.array([_forward_cos(nj, n[0], theta0) for nj in n])
    kz = 2 * np.pi * n * cos / lam

    num = n.size
    r_if = np.empty(num - 1, dtype=complex)
    t_if = np.empty(num - 1, dtype=complex)
    for j in range(num - 1):
        r_if[j], t_if[j] = fresnel(pol, n[j], n[j + 1], cos[j], cos[j + 1])

    m_tot = np.array([[1, r_if[0]], [r_if[0], 1]], dtype=complex) / t_if[0]
    for j in range(1, num - 1):
        delta = kz[j] * d[j]
        # Guard against overflow in very thick absorbing layers.
        if delta.imag > 35:
            delta = delta.real + 35j
        prop = np.array([[np.exp(-1j * delta), 0], [0, np.exp(1j * delta)]])
        inter = np.array([[1, r_if[j]], [r_if[j], 1]], dtype=complex) / t_if[j]
        m_tot = m_tot @ prop @ inter

    r = m_tot[1, 0] / m_tot[0, 0]
    t = 1 / m_tot[0, 0]

    R = abs(r) ** 2
    if pol == "s":
        T = abs(t) ** 2 * (n[-1] * cos[-1]).real / (n[0] * cos[0]).real
    else:
        T = (
            abs(t) ** 2
            * (n[-1] * np.conj(cos[-1])).real
            / (n[0] * np.conj(cos[0])).real
        )
    return TMMResult(r=complex(r), t=complex(t), R=float(R), T=float(T))


def spectrum(pol: str, n_func, d_list, theta0: float, wavelengths):
    """Evaluate R and T over a wavelength grid.

    ``n_func(lam)`` must return the full ``n_list`` for that wavelength, which
    lets each layer be dispersive.

    Returns
    -------
    R, T : ndarray
    """
    wl = np.asarray(wavelengths, dtype=float)
    R = np.empty_like(wl)
    T = np.empty_like(wl)
    for i, lam in enumerate(wl):
        res = coh_tmm(pol, n_func(lam), d_list, theta0, lam)
        R[i], T[i] = res.R, res.T
    return R, T
