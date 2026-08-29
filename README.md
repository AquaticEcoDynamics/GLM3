# GLM 3 : The General Lake Model (legacy 3.x line)

[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![GLM](https://img.shields.io/badge/GLM-3.3.5-lightgrey)](https://github.com/AquaticEcoDynamics/GLM3/tags)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

> [!WARNING]
> **This repository is an archive.**
> `GLM3` holds the legacy GLM 3.x source line and is no longer actively developed.
> Active GLM development has moved to
> [**`AquaticEcoDynamics/GLM`**](https://github.com/AquaticEcoDynamics/GLM), which carries the
> 4.x line. This repository is retained read-only so that the 3.x source, tags and commit
> history stay available for reference and for reproducing published simulations.
>
> **Compatibility notice:** the current
> [`glm-aed`](https://github.com/AquaticEcoDynamics/glm-aed) bundle is **no longer compatible
> with this version of GLM**. It builds against the GLM 4.x source line and its revised AED
> coupling interface. Do not attempt to build current `glm-aed` or `libaed-*` releases against
> this repository — use the GLM-AED releases that predate the 4.x transition if you need a
> 3.x-based coupled model.

<br>

<img src="glm.png" align="right"  width="100" >
The General Lake Model (GLM) is a water balance and one-dimensional vertical stratification
hydrodynamic model. It accounts for the effect of inflows/outflows, mixing and surface heating
and cooling, including the effect of ice cover. It is suited to longer-term investigations
ranging from seasons to decades, and for coupling with biogeochemical models to explore the
role that stratification and vertical mixing has on biogeochemical and ecological dynamics of lakes, reservoirs, ponds and wetlands.

<br>

## Which repository do I want?

| If you want to... | Go to |
|---|---|
| Run or develop the current model | [`GLM`](https://github.com/AquaticEcoDynamics/GLM) (4.x, active) |
| Get executables, examples and a quick start | [`glm-aed`](https://github.com/AquaticEcoDynamics/glm-aed) |
| Reproduce an older 3.x simulation, or consult the 3.x history | this repository |

## Contents of this archive

The final state of this line is **GLM 3.3.5**; the most recent tagged release is
[`GLM_v3.3.0`](https://github.com/AquaticEcoDynamics/GLM3/tags), alongside earlier `v3.0`,
`v3.0.0` and `v3.1.0` tags. The `patches/` directory holds the optional bubbler and restart
source patches for this line — both capabilities are now built into GLM 4.x and the patches are
not carried forward.

The build scripts (`build_glm.sh`, `build_glm.bat`) and platform notes
(`README.Windows`, `README.Macintosh`) are preserved as they were, but are unsupported: they are
known to work against the compiler and library versions current at the time of archiving, and no
fixes will be issued for newer toolchains.

## Reference

Refer to the following paper for a scientific description of the model:

Hipsey, M.R., Bruce, L.C., Boon, C., Busch, B., Carey, C.C., Hamilton, D.P., Hanson, P.C., Read, J.S., de Sousa, E., Weber, M. and Winslow, L.A., 2019. A General Lake Model (GLM 3.0) for linking with high-frequency sensor data from the Global Lake Ecological Observatory Network (GLEON). *Geoscientific Model Development*, **12**(1), pp.473-523. [https://doi.org/10.5194/gmd-12-473-2019](https://doi.org/10.5194/gmd-12-473-2019)
