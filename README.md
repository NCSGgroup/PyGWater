# PyGWater

[![tests](https://github.com/NCSGgroup/PyGWater/actions/workflows/tests.yml/badge.svg)](https://github.com/NCSGgroup/PyGWater/actions/workflows/tests.yml)

**PyGWater** is an open-source, multi-basis toolbox for modelling the Earth's elastic response to surface loading
and for estimating terrestrial water storage from geodetic deformation.

## Overview

PyGWater computes the vertical and horizontal displacement, geoid height, surface gravity, gravity at satellite
altitude and InSAR line-of-sight displacement caused by a surface mass load such as water storage. The response is
represented with **six interchangeable basis functions** that share one set of physical constants and load Love
numbers, following Yang et al. (2025). PyGWater is the loading component of the GravityPython ecosystem of the NCSG
group.

### Key features

- Basis functions: spherical harmonics, finite point-mass Green's function, point-load kernels with Kummer's
  transformation, disk-load Green's function (D-LGF), clustered disk-load Green's function (CD-LGF) and Slepian
  functions
- Observables: vertical, north and east displacement, geoid height, surface gravity (sensor below or above the load),
  gravity disturbance at altitude, InSAR line of sight
- Earth models: load Love numbers of PREM, iasp91, ak135 and their crustal variants (Wang et al., 2012) and REF 6371
  (Bevis et al., 2016), in the CE, CM or CF frame
- Spherical-harmonic analysis with selectable schemes (Driscoll-Healy weights as in SaGEA, the quadrature of
  Yang et al. (2025), least squares, or automatic choice)
- Evaluation on grids or at any points (GNSS stations, InSAR pixels), time series in one pass
- Inversion of GNSS (and InSAR) displacements for water storage with any basis: GNSS time-series cleaning and
  monthly means, regularised least squares, GCV / L-curve with a plausibility safeguard, resolution and covariance;
  closed-loop and checkerboard tests; regional masks from shapefiles

## Repository structure

```
PyGWater/
├── src/                     # the library, imported as `pygwater` (pip install -e .)
│   ├── __init__.py          #   public API: pygwater.LoveNumber, GreenFunctionBasis, RegularizedLeastSquares, ...
│   ├── config.py            #   data root (installation/pygwater_paths.cfg or PYGWATER_DATA)
│   ├── earth/               #   load Love numbers, observables (displacement, geoid, gravity)
│   ├── basis/               #   basis functions behind one interface (base.py: design_matrix, predict, ...)
│   │   ├── spherical_harmonic/  # SH analysis/synthesis (SaGEA-style), SHC/GRID data classes, SH basis
│   │   ├── green/           #     point-load and disk-load kernels; point mass, D-LGF, CD-LGF basis
│   │   └── slepian.py       #     Slepian functions: concentration problem, truncation
│   ├── forward/             #   forward model (any basis), grid convolutions (point mass, disk, clustered disk)
│   ├── inversion/           #   regularised least squares, GCV / L-curve / discrepancy, resolution, covariance
│   ├── observations/        #   GNSS time series (NGL .tenv3), synthetic networks, checkerboards
│   ├── mesh/                #   regular grids, shapefiles and regional masks
│   └── utils/               #   constants, enums, geometry and math tools
├── examples/                # short scripts showing one feature each
├── studies/                 # reproducible studies, one folder each: JGR2025_basis_functions, China, Fu2015_PNW
├── tests/                   # consistency tests (pytest)
├── doc/                     # installation guide and user guide
├── installation/            # environment files; pygwater_paths.cfg (local data-root setting)
└── pyproject.toml           # package definition

PyGWater_external_data/      # the data root, outside the code
├── input/                   # external sources: LLN/, land_mask/, mascon/, shp/, stations/, reference/
├── cache/                   # built by PyGWater and reused: CD-LGF templates, masks
├── output/                  # results, one folder per study
└── temp/                    # scratch
```

## Requirements

- Python >= 3.10 with numpy, scipy, h5py and netCDF4; geopandas and shapely for regional masks; pygmt (optional) for
  the figures of the demos. Tested with Python 3.13.
- The load Love-number files of Wang et al. (2012) (see [doc/installation.md](doc/installation.md)).

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/NCSGgroup/PyGWater.git
cd PyGWater
```

### 2. Install the Python environment and PyGWater

```bash
mamba env create -f installation/environment.yml      # or: pip install -r installation/requirements.txt
conda activate pygwater
pip install -e .                                       # makes `import pygwater` work everywhere (editable)
```

### 3. Provide the data

All data live outside the code under one data root, e.g. `PyGWater_external_data/` with the Love numbers in
`input/LLN/`. Set it in `installation/pygwater_paths.cfg` (copy `installation/pygwater_paths.cfg.example`) or with
`export PYGWATER_DATA=/path/to/PyGWater_external_data`. Caches (`cache/`), results (`output/`) and scratch files
(`temp/`) are written there too. The files and their sources are listed in
[doc/installation.md](doc/installation.md).

## Usage

After `pip install -e .`, PyGWater is imported as `pygwater` from any script or notebook; the main classes are at the
top level (`import pygwater as pgw; pgw.GreenFunctionBasis`), everything else in the subpackages:

```python
import numpy as np
from pygwater.earth.love_numbers import LoveNumber
from pygwater.basis.spherical_harmonic.data_class import GRID
from pygwater.basis.slepian import SlepianBasis
from pygwater.utils.enums import LLN_Data, PhysicalDimensions, Displacement
from pygwater.utils.geomath import GeoMathKit
from pygwater.earth import observables as Observables

# a 10 cm water load on a global 2-degree grid
res = 2.0
lat, lon = GeoMathKit.global_equal_distance_grid(res)
ewh = np.zeros((len(lat), len(lon)))
ewh[np.ix_((lat > 40) & (lat < 50), (lon > 0) & (lon < 10))] = 0.1

# spherical harmonics: analysis of the load, synthesis of any observable
lln = LoveNumber().config(lmax=90, method=LLN_Data.PREM).get_Love_number()
shc = GRID(grid=ewh, lat=lat, lon=lon).analysis(lmax=90, from_type=PhysicalDimensions.EWH,
                                                to_type=PhysicalDimensions.EWH, LLN=lln)
up = shc.synthesis_observable(res, Displacement.Vertical).value                      # m
gravity = shc.synthesis_observable(res, Displacement.Gravity, newtonian='above').value  # m/s^2

# InSAR line of sight (Sentinel-1 ascending)
east = shc.synthesis_observable(res, Displacement.East).value
north = shc.synthesis_observable(res, Displacement.North).value
los = Observables.line_of_sight(east, north, up, Observables.look_vector(39.0, -12.0))

# regional Slepian basis and design matrix for GNSS stations
lo2, la2 = np.meshgrid(*GeoMathKit.global_equal_distance_grid(0.5)[::-1])
mask = (la2 > 35) & (la2 < 70) & (lo2 > -10) & (lo2 < 40)
sb = SlepianBasis(30).configure(mask=mask, lat=la2, lon=lo2, resolution=0.5).truncate(beta_min=0.1)
lln30 = LoveNumber().config(lmax=30, method=LLN_Data.PREM).get_Love_number()
A = sb.design_matrix(np.array([45.0, 52.0]), np.array([5.0, 10.0]), lln30, variable=Displacement.Vertical)

# the same forward model with any basis (common interface of pygwater.basis)
from pygwater.basis.spherical_harmonic.sh_basis import SphericalHarmonicBasis
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.forward.model import ForwardModel

sta_lat, sta_lon = np.array([44.0, 47.5]), np.array([3.0, 6.5])          # e.g. GNSS stations
for basis in (SphericalHarmonicBasis(90), GreenFunctionBasis.from_grid(lat, lon, res, mask=ewh > 0)):
    fm = ForwardModel(basis, lln)
    d = fm.run(fm.parameterize(ewh, lat, lon), sta_lat, sta_lon,
               observables=(Displacement.Vertical, Displacement.North, Displacement.East))
```

Demonstrations and tests:

```bash
python examples/slepian_truncation.py       # Slepian truncation error over a regional mask
python -m pytest tests -q             # consistency between basis functions, frames, schemes and observables
```

All options (basis functions, observables, analysis schemes, conventions, limitations) are described in
[doc/user_guide.md](doc/user_guide.md).

## Related software (GravityPython)

| Toolbox | Scope |
|---|---|
| [PyHawk](https://github.com/NCSGgroup/PyHawk) | GRACE/GRACE-FO gravity recovery, Level-1B to Level-2 |
| [SaGEA](https://github.com/NCSGgroup/SaGEA) | Level-2 post-processing and error assessment |
| [PyGLDA](https://github.com/AAUGeodesyGroup/PyGLDA) | GRACE data assimilation into hydrological models |
| [SAGEA-fluid](https://github.com/NCSGgroup/SAGEA-fluid) | geophysical fluids, degree-1, Earth rotation |
| **PyGWater** | surface loading: deformation, gravity and water storage from geodetic observations |

## Citation

If you use PyGWater in your research, please cite:

> Yang, F., Shen, Y., Zhang, W., & Forootan, E. (2025). On the choice of basis functions for modeling Earth's elastic
> deformation due to surface loading. *Journal of Geophysical Research: Solid Earth*, 130, e2025JB031662.
> https://doi.org/10.1029/2025JB031662

## License

PyGWater is licensed under the [MIT License](LICENSE). The load Love numbers and Green's functions used as input data
(Wang et al., 2012; Bevis et al., 2016) are third-party data and are not covered by this license; please cite their
sources.

## Contact

Fan Yang — [fany@plan.aau.dk](mailto:fany@plan.aau.dk)
Geodesy Group, Aalborg University (AAU)
