"""thinfilm: transfer-matrix modelling of optical thin-film coatings."""

from .design import (
    bragg_mirror,
    quarter_wave_stack_reflectance,
    quarter_wave_thickness,
    single_layer_ar,
)
from .materials import bk7, cauchy, fused_silica, sellmeier
from .tmm import TMMResult, coh_tmm, fresnel, spectrum

__version__ = "0.1.0"

__all__ = [
    "TMMResult",
    "coh_tmm",
    "fresnel",
    "spectrum",
    "sellmeier",
    "cauchy",
    "fused_silica",
    "bk7",
    "quarter_wave_thickness",
    "single_layer_ar",
    "bragg_mirror",
    "quarter_wave_stack_reflectance",
]
