# thinfilm-optics

![tests](https://github.com/luisrowork1992-netizen/thinfilm-optics/actions/workflows/tests.yml/badge.svg)

A small, dependency-light Python package for modelling **optical thin-film coatings** with the coherent transfer-matrix method (TMM).

## Features

- Reflectance and transmittance of arbitrary planar multilayers for **s and p polarisation** at any angle of incidence
- Handles **absorbing layers** (complex index `n + ik`) and **total internal reflection**
- Design helpers: quarter-wave layers, single-layer AR coatings, (HL)^N Bragg mirrors and their closed-form peak reflectance
- Sellmeier dispersion for fused silica (Malitson) and Schott N-BK7, plus a generic Cauchy model

## Install

```bash
git clone https://github.com/luisrowork1992-netizen/thinfilm-optics.git
cd thinfilm-optics
pip install -e .
```

## Quick start

```python
import math
import thinfilm as tf

# Uncoated glass at normal incidence: ~4 % reflection
print(tf.coh_tmm("s", [1.0, 1.5], [math.inf, math.inf], 0.0, 0.55).R)   # 0.04

# Ideal quarter-wave AR coating on n = 1.52 glass at 550 nm
n_list, d_list = tf.single_layer_ar(1.0, 1.52, lam0=0.55)
print(tf.coh_tmm("s", n_list, d_list, 0.0, 0.55).R)                    # ~0

# 8-pair TiO2/SiO2 Bragg mirror for 1064 nm
n_list, d_list = tf.bragg_mirror(1.0, 2.35, 1.46, 1.52, pairs=8, lam0=1.064)
print(tf.coh_tmm("p", n_list, d_list, 0.0, 1.064).R)                   # ~0.9987
```

Lengths can be in any unit, as long as the thicknesses and the wavelength use the same one. The Sellmeier functions expect micrometres.

## Conventions

| Item | Convention |
| --- | --- |
| Complex index | `n = n' + ik`, with `k >= 0` for absorption |
| Layer order | incident medium, layers..., exit medium |
| Outer thicknesses | `math.inf` |
| Angle | radians, measured in the incident medium |

## Notebook

[`notebooks/coating_design.ipynb`](notebooks/coating_design.ipynb) compares bare BK7, a MgF2 AR coating and Bragg mirrors with 2, 4 and 8 pairs. It also plots the s/p reflectance against angle, including Brewster's angle.

## Tests

```bash
pip install pytest
pytest
```

The tests check the model against analytic results: 4 % Fresnel reflection, a reflectionless ideal AR coating, the half-wave absentee layer, Brewster's angle, total internal reflection, energy conservation (R + T = 1) for lossless stacks, and the closed-form Bragg-mirror reflectance.

## References

- H. A. Macleod, *Thin-Film Optical Filters*, 5th ed., CRC Press.
- I. H. Malitson, "Interspecimen comparison of the refractive index of fused silica," *J. Opt. Soc. Am.* 55, 1205 (1965).
- S. J. Byrnes, "Multilayer optical calculations," arXiv:1603.02720.

## License

MIT
