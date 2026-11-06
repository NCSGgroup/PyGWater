# PyGWater — user guide

> **PyGWater**: a multi-basis toolbox for modelling the elastic response of the Earth to surface loading
> (displacement, geoid, gravity, InSAR line of sight) and for estimating terrestrial water storage from geodetic
> deformation.
> Status: development version (not yet released).
> Reference: Yang, F., Shen, Y., Zhang, W., & Forootan, E. (2025). On the choice of basis functions for modeling
> Earth's elastic deformation due to surface loading. *Journal of Geophysical Research: Solid Earth*, 130,
> e2025JB031662. https://doi.org/10.1029/2025JB031662

---

## Table of Contents

1. [Concept](#1-concept)
2. [Installation and data](#2-installation-and-data)
3. [Conventions](#3-conventions)
4. [Load Love numbers](#4-load-love-numbers)
5. [Observables](#5-observables)
6. [Basis functions](#6-basis-functions)
   - 6.1 [Spherical harmonics](#61-spherical-harmonics)
   - 6.2 [Finite point-mass Green's function (finite LGF)](#62-finite-point-mass-greens-function-finite-lgf)
   - 6.3 [Point-load kernels with Kummer's transformation](#63-point-load-kernels-with-kummers-transformation)
   - 6.4 [Disk-load Green's function (D-LGF)](#64-disk-load-greens-function-d-lgf)
   - 6.5 [Clustered disk-load Green's function (CD-LGF)](#65-clustered-disk-load-greens-function-cd-lgf)
   - 6.6 [Slepian functions](#66-slepian-functions)
   - 6.7 [Common interface and ForwardModel](#67-common-interface-and-forwardmodel)
7. [Choosing a basis function](#7-choosing-a-basis-function)
8. [InSAR line of sight](#8-insar-line-of-sight)
9. [Building blocks for inversion](#9-building-blocks-for-inversion)
10. [Regional masks](#10-regional-masks)
11. [Physical constants](#11-physical-constants)
12. [Tests and demos](#12-tests-and-demos)
13. [Known limitations](#13-known-limitations)
14. [References](#14-references)

---

## 1. Concept

Every quantity in PyGWater is a linear response of a spherically symmetric, non-rotating, elastic, isotropic
(SNREI) Earth to a thin surface load, described by the load Love numbers h_n, l_n, k_n. With the load given as
equivalent water height (EWH) E and F_n = 3 rho_w / (rho_e (2n+1)), the response is, degree by degree,

| Observable | Factor per degree | Evaluated with |
|---|---|---|
| Vertical displacement | F_n h_n | Y_nm |
| Horizontal displacement (north, east) | F_n l_n | surface gradient of Y_nm |
| Geoid height | F_n (1 + k_n) | Y_nm |
| Surface gravity (gravimeter) | (g/a) F_n [N_n + (n+1) k_n - 2 h_n] | Y_nm |
| Gravity disturbance at height H | (g/a) F_n (n+1)(1 + k_n)(a/(a+H))^(n+2) | Y_nm |

(N_n is the Newtonian term, see [Section 5](#5-observables).)

The same physics can be represented with different **basis functions**: spherical harmonics, Green's functions
(point mass, disk, clustered disks) and Slepian functions. Yang et al. (2025) show that the choice and configuration
of the basis can cause errors as large as those from the assumed Earth structure. PyGWater implements these bases
with **identical constants and Love numbers**, so they can be compared and swapped.

```
          src_EarthModel: Love numbers (LoveNumber.py) + observable table (Observables.py)
                                      |
   src_Basis: one interface (base.Basis: design_matrix, predict, load_synthesis, load_parameters)
     +--------------------------------+-------------------------------+
     |                                |                               |
 SphericalHarmonic                GreenFunction                     Slepian
 (SaGEA-style analysis/synthesis, (point-load and disk-load        (regional, truncated;
  SphericalHarmonicBasis)          kernels, GreenFunctionBasis)     SlepianBasis)
     |                                |                               |
     +--------------------------------+-------------------------------+
                                      |
   src_Forward: ForwardModel (any basis -> any observable at any point; stacked operators for inversion)
                + grid convolutions (point mass, disk, clustered disk) for large global grids
```

---

## 2. Installation and data

See [installation.md](installation.md) for the environment, the external data (Love numbers in `<data root>/LLN/`,
data root set by `PYGWATER_DATA`, default `<code root>/data`) and how to run the code.

The packages are imported as `src_<name>` from the code root:

| Package | Content |
|---|---|
| `src_auxiliary` | constants, enums, geometry and math tools, file utilities |
| `src_EarthModel` | load Love numbers (`LoveNumber.py`) and observable definitions (`Observables.py`) |
| `src_Basis` | common interface `base.py` (class `Basis`) |
| `src_Basis.SphericalHarmonic` | spherical harmonics: `DataClass.py` (SHC, GRID), `Harmonic.py`, `ConvertSHC.py`, `SHBasis.py` |
| `src_Basis.GreenFunction` | Legendre polynomials, point-load (`PointLoad.py`) and disk-load (`DiskLoad.py`) kernels, `GFBasis.py` |
| `src_Basis.Slepian` | `Slepian.py` (class `SlepianBasis`, a subclass of `SphericalHarmonicBasis`) |
| `src_Forward` | `ForwardModel.py` (any basis) and grid convolutions `PointMass.py`, `Disk.py`, `DiskCluster.py` |
| `src_LoadMesh` | `RegularGrid.py`, `ShapeAndMask.py` |
| `src_demo` | demonstration scripts; `paper_JGR2025.py` reproduces the figures of Yang et al. (2025) |

---

## 3. Conventions

- **Grids.** Global regular grids are cell-centred: latitudes `-90 + res/2 ... 90 - res/2` (south to north),
  longitudes `-180 + res/2 ... 180 - res/2`, as returned by `GeoMathKit.global_equal_distance_grid(res)`.
  Spherical-harmonic analysis expects this layout. Masks written by `src_LoadMesh` run north to south; the Slepian
  module accepts any ordering.
- **Units.** SI throughout: displacement and geoid in metres, gravity in m/s^2 (1 microGal = 1e-8 m/s^2), angles in
  degrees, areas in m^2, loads in metres of EWH unless a function says otherwise.
- **Time series.** Loads may carry a time dimension: `(ntime, nlat, nlon)` grids, `(ntime, L+1, L+1)`
  coefficients, `(Ncell, Ntime)` in the Green's-function forward models.
- **Signs.** Up, north and east are positive. Horizontal kernels (`getHorizental`) are positive *away* from the load
  (Farrell, 1972); geographic components are obtained with `Displacement.North` / `Displacement.East`. Gravity is
  positive for an increase of gravity. Line of sight is positive towards the satellite.
- **Reference frame.** The Love-number files are in the CE frame; use `LoveNumber.convert()` for CM or CF
  ([Section 4](#4-load-love-numbers)). Degree 1 matters for geocentre-sensitive comparisons such as GNSS.

---

## 4. Load Love numbers

```python
from src_EarthModel.LoveNumber import LoveNumber
from src_auxiliary.EnumClasses import LLN_Data, LLN_variable, Frame

lln = LoveNumber().config(lmax=1000, method=LLN_Data.PREM).get_Love_number()
h = lln.LLN[LLN_variable.h]          # h_0 ... h_1000 (numpy array)
lln.convert(target=Frame.CF)         # degree-1 numbers in the CF frame (also CE, CM)
```

| Option | Values | Notes |
|---|---|---|
| `config(lmax=...)` | int | Wang: <= 13 000; complete files: <= 46 000; REF: <= 40 000 |
| `config(method=...)` | `LLN_Data.Wang` (default), `PREM`, `iasp91`, `ak135`, `PREMhard`, `PREMsoft`, `iasp91hard`, `ak135hard`, `REF` | Earth model |
| `convert(target=...)` | `Frame.CE`, `Frame.CM`, `Frame.CF` | Blewitt (2003): CM: x_1 - 1; CF: x_1 - (h_1 + 2 l_1)/3, using the degree-1 numbers of the loaded model |

Data are read on `get_Love_number()` or on first access to `.LLN`. Degree 0: h_0 = -0.13205 (PREM), l_0 = k_0 = 0.

---

## 5. Observables

Observables are members of `src_auxiliary.EnumClasses.Displacement` (alias `Observables.Observable`). Their
definitions live in **one** table, `src_EarthModel/Observables.py`, which spherical harmonics and Slepian functions
read. Green's-function kernels implement each observable individually.

| `variable` | Quantity | Unit | Extra options | SH | Slepian | Point LGF | D-LGF / CD-LGF |
|---|---|---|---|---|---|---|---|
| `Vertical` | vertical displacement | m | | yes | yes | yes | yes |
| `North`, `East` | horizontal displacement | m | | yes | yes | yes | yes |
| `Horizontal` | horizontal displacement along the great circle, positive away from the load (kernel level) | m | | | | yes | yes |
| `Geoheight` | geoid height | m | | yes | yes | yes | yes |
| `Gravity` | surface gravity change at the deforming surface | m/s^2 | `newtonian` | yes | yes | yes | not implemented |
| `GravityDisturbance` | gravity disturbance at height H above the sphere (satellite, airborne) | m/s^2 | `height` (m) | yes | yes | | |

**Surface gravity** (Farrell, 1972) = Newtonian attraction of the load + attraction of the deformed Earth
((n+1) k_n) + free-air effect of the vertical displacement (-2 h_n):

| `newtonian=` | Newtonian term N_n | Use |
|---|---|---|
| `'below'` (default) | -n | sensor beneath the load (e.g. atmosphere); Farrell's and Wang et al.'s convention |
| `'above'` | n + 1 | sensor on top of the load (e.g. groundwater, soil moisture underneath) |
| `'none'` | 0 | elastic/deformation part only |

`'below'` and `'above'` differ only at the sensor, by 4 pi G rho_w E (about 83.9 microGal per metre of water).
Away from the load, Green's functions are identical for both. The factor g/a is evaluated with g = GM/a^2 of the
model Earth, so that (g/a) F_n = 4 pi G rho_w/(2n+1) exactly.

Low-level access:

```python
from src_EarthModel import Observables
Observables.love_factor(variable, lln, n, newtonian='below', height=0.0)   # per-degree factor X_n
Observables.kind(variable)                                                  # 'value' | 'north' | 'east'
Observables.sh_synthesis(load_cqlm, load_sqlm, lln, variable, lat, lon, **kwargs)
```

---

## 6. Basis functions

### 6.1 Spherical harmonics

Analysis of the **load** grid, then synthesis of any observable (synthesis only; observables are never analysed):

```python
import numpy as np
from src_Basis.SphericalHarmonic.DataClass import GRID
from src_auxiliary.EnumClasses import PhysicalDimensions, Displacement
from src_auxiliary.GeoMathKit import GeoMathKit

res = 2.0
lat, lon = GeoMathKit.global_equal_distance_grid(res)      # cell centres, south->north, -180->180
ewh = np.zeros((len(lat), len(lon)))                       # load in metres of equivalent water height
ewh[np.ix_((lat > 40) & (lat < 50), (lon > 0) & (lon < 10))] = 0.1   # a 10 cm water load in Europe

lln = LoveNumber().config(lmax=90, method=LLN_Data.PREM).get_Love_number()
shc = GRID(grid=ewh, lat=lat, lon=lon).analysis(lmax=90, from_type=PhysicalDimensions.EWH,
                                                to_type=PhysicalDimensions.EWH, LLN=lln)   # method='auto'

up = shc.synthesis_observable(res, Displacement.Vertical).value          # metre
north = shc.synthesis_observable(res, Displacement.North).value
geoid = shc.synthesis_observable(res, Displacement.Geoheight).value
grav = shc.synthesis_observable(res, Displacement.Gravity, newtonian='above').value       # m/s^2
grav_sat = shc.synthesis_observable(res, Displacement.GravityDisturbance, height=450e3).value
```

`synthesis_observable(grid_space, variable, lln=None, **kwargs)` uses the Love numbers attached at analysis time
(`LLN=` argument) unless `lln` is given; `.value` has shape `(ntime, nlat, nlon)`.

**Analysis methods** (`GRID.analysis(..., method=)`, `GRID.to_SHC(..., method=)`, `Harmonic.analysis(grid, method=)`):

| `method` | Scheme | Exact for band-limited fields when | Use |
|---|---|---|---|
| `'auto'` (default) | `'dh'` if lmax <= nlat/2 - 1, otherwise `'midpoint'` | as the chosen scheme | general use |
| `'dh'` | SaGEA scheme: longitude sums / nlon, Driscoll-Healy latitude weights | lmax <= nlat/2 - 1 (44 on 2 deg, 89 on 1 deg); warns beyond | accurate analysis within the sampling limit |
| `'midpoint'` (alias `'quadrature'`) | two-step quadrature of Sneeuw (1994), Eqs. 3-5 of Yang et al. (2025) | never exactly (zonal terms biased) | block-mean grids analysed up to Nmax = pi/dlambda; equivalent to the finite LGF (Sect. 4.2 of the paper) |
| `'least_squares'` | longitude sums, per-order least squares in latitude | lmax < nlat | point-type grid values that must be reproduced exactly |

The term of the Nyquist order (delta_{m,Nmax} in Eq. 5) is applied in all schemes. Truncation: on a grid of spacing
dlambda, choose `lmax = 180/dlambda` (Eq. 6 of the paper); higher values alias, lower values over-smooth.

The original unit-conversion route is still available: `shc.synthesis(grid_space, from_type=PhysicalDimensions.EWH,
to_type=PhysicalDimensions.VerticalDisplacement | Geoid | HorizontalDisplacementNorth | HorizontalDisplacementEast)`.

### 6.2 Finite point-mass Green's function (finite LGF)

Truncated series G(alpha) = a/M sum_{n<=N} h_n P_n(cos alpha) (Eq. 7). Equivalent to spherical harmonics if N
matches the load resolution (N = 180/res); larger N aliases (Sect. 4 of the paper).

```python
from src_Forward.PointMass import GFA_displacement, Grids_generation
from src_auxiliary.constants import EarthConstant

grids = Grids_generation.Equal_angular_distance(resolution=res)       # 'lat', 'lon', 'area' (m^2), flattened
grids['EWH'] = ewh.reshape(-1, 1) * EarthConstant.rhow                # kg/m^2 !  shape (Ncell, Ntime)
lln90 = LoveNumber().config(lmax=90, method=LLN_Data.PREM).get_Love_number()   # lmax = 180/res (Eq. 6)
fm = GFA_displacement(lln90)
fm.configure(grids=grids)
stations = {'lat': np.array([45.0, 52.0]), 'lon': np.array([5.0, 10.0])}
u = fm.evaluation(points=stations, variable=Displacement.Vertical)      # (Nstation, Ntime), metre
g = fm.evaluation(points=stations, variable=Displacement.Gravity)       # m/s^2
```

> **Note.** In the point-mass forward models `grids['EWH']` is a surface **mass density** (kg/m^2), i.e. EWH times
> `EarthConstant.rhow`.

| Class / method | Options |
|---|---|
| `GFA_displacement(lln, verbose=False)` | general: any sources, any evaluation points |
| `.configure(grids)` | `grids = {'lat', 'lon', 'area', 'EWH'}` |
| `.evaluation(points, variable)` | `points = {'lat', 'lon'}`; `variable`: `Vertical`, `North`, `East`, `Horizontal`, `Geoheight`, `Gravity` |
| `GFA_regular_grid(lln, verbose=False)` | fast version: one kernel per latitude row, shifted in longitude; sources **and** points must be the same global grid |
| `.evaluation(points, variable, resolution=2)` | `resolution` = grid spacing in degrees |

```python
from src_Forward.PointMass import GFA_regular_grid
lo2, la2 = np.meshgrid(lon, lat)
fm = GFA_regular_grid(lln90)
fm.configure(grids=grids)
u_grid = fm.evaluation(points={'lat': la2.ravel(), 'lon': lo2.ravel()},
                       variable=Displacement.Vertical, resolution=res)   # points must be the same global grid
```

### 6.3 Point-load kernels with Kummer's transformation

Kernels per kg of load as a function of angular distance (degrees). `LGF_truncation` is the plain truncated series;
`LGF_Kummer` adds the closed-form asymptotic part (Eq. 14), which converges for any N (Fig. 3 of the paper).

```python
from src_Basis.GreenFunction.PointLoad import LGF_Kummer, LGF_truncation
lln10k = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
k = LGF_Kummer(lln10k).configure(residence=np.array([0.1, 1.0, 10.0]))   # angular distance, degree
print(k.getVertical(), k.getHorizental(), k.getGeoidheight(), k.getGravity())  # per kg of load
```

| Method | Returns (per kg) | Options |
|---|---|---|
| `getVertical()` | vertical displacement, degree 0 included | `LGF_Kummer.getVertical_exclude_d0()` |
| `getHorizental()` | horizontal displacement, positive away from the load | |
| `getGeoidheight()` | geoid height (`LGF_Kummer`: degree 0 excluded) | |
| `getGravity(newtonian='below', include_degree0=True)` | surface gravity | with `include_degree0=False`, `LGF_Kummer` reproduces the gravity Green's function of Wang et al. (2012) to 0.1 % |

Distances below 1e-5 deg and above 180 - 1e-5 deg are clipped to avoid singularities. `GeoMathKit.kernel_at_points(obj,
load_lat, load_lon, points_lat, points_lon, variable)` evaluates any kernel for a load at a position, including the
north/east projection u_N = -v cos(az), u_E = -v sin(az) (az = azimuth from the point to the load).

### 6.4 Disk-load Green's function (D-LGF)

A uniform disk of thickness 1 m EWH (Eqs. 16-21). The expansion must satisfy Eq. 24:
N >= max(2 x 360/radius, 10 000) (`DiskLoad.required_lmax(radius)`); a `UserWarning` is issued otherwise.

```python
from src_Basis.GreenFunction.DiskLoad import DiskLoad
from src_Forward.Disk import GFA_displacement as DiskForward, grid2radius
print(DiskLoad.required_lmax(0.5))
radius = grid2radius(lat_center=grids['lat'], grid_size=res)[0]           # equal-area disk radius, degree
dgrid = {'lat': grids['lat'], 'lon': grids['lon'], 'radius': radius, 'EWH': ewh.reshape(-1, 1)}  # EWH in metre
lln2k = LoveNumber().config(lmax=2000, method=LLN_Data.PREM).get_Love_number()   # (Eq. 24 wants 10 000)
dfm = DiskForward(lln2k)
dfm.configure(grids=dgrid)                                              # cf=None: no far-field truncation
u_north = dfm.evaluation(points=stations, variable=Displacement.North)
```

| Class / method | Options |
|---|---|
| `DiskLoad(lln, radius, thickness=1, check_lmax=True)` | `.configure(residence)`, `.getVertical()`, `.getHorizental()`, `.getGeoidheight()`, `.release()` |
| `DiskLoad_constrain(lln, radius, thickness=1, constrain_factor=10)` | sets the response beyond `constrain_factor x radius` to **zero** (approximation) |
| `Disk.GFA_displacement(lln, verbose=False)` | `.configure(grids={'lat','lon','radius','EWH'}, cf=None)` (`EWH` in metres); `cf=K` uses `DiskLoad_constrain`. `.evaluation(points, variable)` |
| `Disk.GFA_regular_grid(lln)` | as 6.2, `.evaluation(points, variable, resolution=2)` |
| `grid2radius(lat_center, grid_size)` | equal-area disk radius of a cell (degree, metre) |

### 6.5 Clustered disk-load Green's function (CD-LGF)

Each cell is split into sub-cells, each represented by a small disk (Sect. 7 of the paper). Two routes:

**At any points (recommended):** `GreenFunctionBasis` with `kernel='cluster_disk'` (see 6.7). The D-LGF of every
distinct sub-disk radius is computed once on a dense grid of distances (`DiskLoad.DiskLoadTable`: one Legendre
series per radius, step radius/200 across the disk and its edge) and interpolated; the interpolation error is about
1e-5 of the kernel peak (vertical, geoid) and 3e-5 (horizontal), and no far-field truncation is applied. A regional
2-deg mascon grid of 25 cells with 0.1-deg sub-cells (N = 13 000) takes about 4 s for 200 stations.

```python
from src_Basis.GreenFunction.GFBasis import GreenFunctionBasis
lln13k = LoveNumber().config(lmax=13000, method=LLN_Data.PREM).get_Love_number()      # Eq. 24 for 0.1-deg sub-cells
lat, lon = GeoMathKit.global_equal_distance_grid(2.0)
lo2, la2 = np.meshgrid(lon, lat)
cells = (la2 > 40) & (la2 < 50) & (lo2 > 0) & (lo2 < 10)                                # 25 mascons of 2 deg
cd = GreenFunctionBasis.from_grid(lat, lon, 2.0, mask=cells, kernel='cluster_disk', subgrid=0.1)
sta_lat, sta_lon = np.array([44.0, 47.5, 52.0]), np.array([3.0, 6.5, 12.0])
A = cd.design_matrix(sta_lat, sta_lon, lln13k, Displacement.Vertical)                  # (N_sta, N_cell)
u = A @ cd.load_parameters(np.where(cells, 0.1, 0.0), lat, lon)                       # 10 cm load -> m
```

Supported observables: Vertical, North, East, Geoheight (no surface gravity, as for the D-LGF). `cd.sub_cells(j)`
returns the sub-disk centres and radii of cell j.

**Templates on a fixed set of points** (as used for the figures of the paper): response templates are computed once
per latitude row, saved, and reused by rolling along longitude on a regular grid.

```python
from src_Forward.DiskCluster import cluster_disk, GFA_displacement_from_template
cd = cluster_disk(lln2k)
cd.configure(grid={'lat': 45.0, 'lon': 5.0, 'delta_lat': 2.0, 'delta_lon': 2.0}, small_disk=0.5)
tpl = cd.template(points=stations, variable=Displacement.Vertical, postfix='demo', out_dir='res/LGD_template')
```

| Class / method | Options |
|---|---|
| `cluster_disk(lln, verbose=False)` | `.configure(grid={'lat','lon','delta_lat','delta_lon'}, small_disk=0.1)`; `small_disk` is the **sub-grid size** in degrees (radius ~ small_disk/sqrt(pi) at the equator). The paper recommends 0.1 (<= 0.2). |
| `.template(points, variable, postfix, out_dir=None)` | response to 1 m EWH over the cell; saved as `LGD_<postfix>.npy` (Vertical) or `LGD_<Variable>_<postfix>.npy` |
| `cluster_disk_fast` | approximates the far field (> 8 deg) by the response of ONE sub-disk (index 100, so the cell needs > 100 sub-cells) for every sub-disk; prefer `cluster_disk` or the basis above |
| `GFA_displacement_from_template(template_dir)` | `.evaluation_general(grids={'EWH','Template'}, variable)`; `.evaluation_regular_grid(grids, dim=nlat, variable)`; `EWH` in metres |

Templates for the 0.1 deg sub-grid need lmax >= ~12 800 at the equator (Eq. 24), more poleward.

### 6.6 Slepian functions

Band-limited functions concentrated in a region (Sect. 8 of the paper). Built from any 0/1 (or fractional) mask on a
regular grid; use a mask clearly finer than 180/L degrees.

```python
from src_Basis.Slepian.Slepian import SlepianBasis
L = 30
lo2, la2 = np.meshgrid(*GeoMathKit.global_equal_distance_grid(0.5)[::-1])
mask = (la2 > 35) & (la2 < 70) & (lo2 > -10) & (lo2 < 40)                 # any 0/1 (or fractional) mask
sb = SlepianBasis(L).configure(mask=mask, lat=la2, lon=lo2, resolution=0.5)
print(sb.shannon, (sb.beta > 0.1).sum())
sb.truncate(beta_min=0.1)                                                  # or Lc=..., or shannon=True
llnL = LoveNumber().config(lmax=L, method=LLN_Data.PREM).get_Love_number()
shcL = GRID(grid=ewh, lat=lat, lon=lon).analysis(lmax=L, from_type=PhysicalDimensions.EWH,
                                                 to_type=PhysicalDimensions.EWH, LLN=llnL)
c, s = shcL.get_cs2d()
u_sl = sb.forward(c, s, llnL, stations['lat'], stations['lon'], variable=Displacement.Vertical)   # (ntime, Nstation)
```

| Method | Purpose |
|---|---|
| `SlepianBasis(lmax)` | band limit L; coefficients ordered by degree, then cos m = 0..n, then sin m = 1..n |
| `.configure(mask, lat, lon, resolution=None)` | concentration matrix (Eq. 28) and eigen-decomposition; `.beta` (descending), `.G` |
| `SlepianBasis.from_mask_h5(lmax, h5file, dataset='basin')` | build from a `src_LoadMesh` mask file |
| `.shannon` | Shannon number (L+1)^2 A/(4 pi) |
| `.truncate(Lc=None, beta_min=None, shannon=False)` | number of retained functions (exactly one option) |
| `.forward(load_c, load_s, lln, lat, lon, variable, **kw)` | any observable of Section 5 at points (Eqs. 29-31) |
| `.displacement_coefficients(...)`, `.synthesis(S, lat, lon, component)` | the two steps of `forward` |
| `.design_matrix(lat, lon, lln, variable, **kw)` | inversion operator A (Section 9) |
| `.design_matrix_los(lat, lon, lln, look)` | A for InSAR line of sight |
| `.cs2d_to_vec`, `.vec_to_cs2d`, `.Ymatrix(lat, lon, component)` | helpers |

With all (L+1)^2 functions retained, the Slepian solution equals the spherical-harmonic one exactly (Eq. 30).
The concentration matrix has (L+1)^4 entries: practical up to L ~ 80-100.

### 6.7 Common interface and ForwardModel

Every basis is a linear parameterisation of the load, EWH(Omega) = sum_j x_j phi_j(Omega), with the operator that maps
x to any observable at any point: d = A x. All bases derive from `src_Basis.base.Basis`:

| Method | Returns |
|---|---|
| `n_param` | number of parameters |
| `design_matrix(lat, lon, lln, variable, **kw)` | A, shape (N, n_param), observable per metre EWH |
| `predict(x, lat, lon, lln, variable, **kw)` | A x, shape (N,) or (ntime, N) |
| `design_matrix_los(lat, lon, lln, look)` | A for InSAR line of sight |
| `load_synthesis(x, lat, lon)` | EWH [m] at points |
| `load_parameters(ewh, lat, lon)` | x of a gridded EWH load (nlat, nlon) or (ntime, nlat, nlon) |

| Basis | Parameters x | Constructor |
|---|---|---|
| `SphericalHarmonicBasis(lmax)` | SH coefficients of EWH, (L+1)^2 (ordering: degree, cosine m=0..n, sine m=1..n) | `src_Basis/SphericalHarmonic/SHBasis.py` |
| `SlepianBasis(lmax).configure(mask, lat, lon).truncate(...)` | Slepian coefficients of EWH, Lc | `src_Basis/Slepian/Slepian.py` |
| `GreenFunctionBasis.from_grid(lat, lon, res, mask=None, kernel='point_mass' \| 'disk' \| 'cluster_disk', subgrid=0.1)` | EWH per cell | `src_Basis/GreenFunction/GFBasis.py` |

`load_parameters` performs SH analysis (`method='auto'`, see 6.1) for the SH basis, SH analysis followed by projection
onto the retained Slepian functions for the Slepian basis, and samples the grid at the cell centres for the
Green's-function basis. `GreenFunctionBasis.load_synthesis` is piecewise constant (value of the containing cell,
0 outside). For the disk kernel the Love numbers must satisfy Eq. 24, and gravity is not available (6.4).

`src_Forward.ForwardModel` wraps any basis:

```python
from src_Basis.SphericalHarmonic.SHBasis import SphericalHarmonicBasis
from src_Basis.GreenFunction.GFBasis import GreenFunctionBasis
from src_Forward.ForwardModel import ForwardModel

basis = GreenFunctionBasis.from_grid(lat, lon, res, mask=ewh > 0, kernel='point_mass')   # or SphericalHarmonicBasis(90)
fm = ForwardModel(basis, lln)
x = fm.parameterize(ewh, lat, lon)
d = fm.run(x, sta_lat, sta_lon, observables=(Displacement.Vertical, Displacement.North, Displacement.East))
A = fm.design_matrix(sta_lat, sta_lon, observables=(Displacement.Vertical, 'LOS'), look=look)  # stacked rows
```

`run` returns a dict {observable: values}; `design_matrix` stacks the observables row-wise (observable by
observable). The grid convolutions in `src_Forward/PointMass.py`, `Disk.py` and `DiskCluster.py` remain the faster
choice for global grids; `DiskCluster.py` keeps the template route of the CD-LGF.

---

## 7. Choosing a basis function

Following Yang et al. (2025), first decide what a grid value represents:

| Grid value represents | Typical data | Recommended basis | Configuration |
|---|---|---|---|
| I. block mean | reanalysis, hydrological models | spherical harmonics (global); finite LGF (regional) | Nmax = 180/res |
| II. point value | in-situ gauges | infinite LGF (not implemented in PyGWater) | |
| III. uniform value per cell (mascon) | GRACE mascons, constrained GNSS inversion | CD-LGF (sub-grid 0.1-0.2 deg); D-LGF only for cells <= ~0.1 deg | N >= max(720/radius, 10 000) |

Slepian functions are a regional alternative to spherical harmonics; their truncation (Lc) adds an approximation
error of several millimetres for vertical displacement (Sect. 8.3).

---

## 8. InSAR line of sight

```python
from src_EarthModel import Observables
look = Observables.look_vector(incidence=39.0, heading=-12.0)        # Sentinel-1 ascending, right-looking
east = shc.synthesis_observable(res, Displacement.East).value
los = Observables.line_of_sight(east, north, up, look)                # positive = towards the satellite
```

| Function | Options |
|---|---|
| `look_vector(incidence, heading, right_looking=True)` | degrees; returns the ground-to-satellite unit vector (E, N, U); arrays give (3, N) |
| `line_of_sight(east, north, up, look)` | projection of a displacement field |
| `<any basis>.design_matrix_los(lat, lon, lln, look)` | LOS design matrix |

---

## 9. Building blocks for inversion

Forward models are linear, d = A x, and every basis provides A through the common interface (6.7):

- `ForwardModel(basis, lln).design_matrix(lat, lon, observables=(Vertical, North, East, 'LOS'), look=...)` stacks
  multi-component data (GNSS three components, InSAR LOS, gravity) for one parameter vector.
- Regional inversions typically use `SlepianBasis` (parameters: load Slepian coefficients) or `GreenFunctionBasis`
  (parameters: mascon EWH); map estimates back with `basis.load_synthesis(x, lat, lon)`.
- CD-LGF: `GreenFunctionBasis(..., kernel='cluster_disk', subgrid=0.1)` gives A directly at the stations; the
  templates of `src_Forward/DiskCluster.py` are pre-computed columns of A on a fixed grid.

A dedicated inversion module (regularisation, variance components, normal-equation export) is planned
(`src_Inversion`).

---

## 10. Regional masks

`src_LoadMesh/ShapeAndMask.py` (requires geopandas and shapely):

| Class / function | Purpose |
|---|---|
| `basinMask(mask_resolution, new_mask_name, save_dir)` | `.load_shp(shp)`, `.shp_to_mask_global()`, `.shp_to_mask_crop(crop_box=None)`; writes an HDF5 file with `basin`, `sub_basin_<id>`, `lat`, `lon`, `resolution` |
| `basin2grid_shp(grid=(3, 3))` | splits a region into coarse cells; `.delete_tiny_grids()`, `.delete_ocean_grid()` |
| `getMask_flow(basin_name, coarse_cell=3, small_cell=0.1, is_ocean_removed=False)` | coarse cells plus fine masks in one go |
| `box2shp(box_area)` | rectangle to shapefile |

---

## 11. Physical constants

`src_auxiliary/constants.py` is the single source of truth, used by every basis:

| Constant | Value | Note |
|---|---|---|
| `radiusm` | 6 378 136.3 m | |
| `rhoear` | 5517 kg/m^3 | |
| `Mass` | 4/3 pi a^3 rho_e | consistent with `rhoear` and `radiusm` (required for basis equivalence) |
| `rhow` | 1000 kg/m^3 | fresh water, as in Yang et al. (2025); change here only |
| `ggg` | 6.67384e-11 | |
| `GM` | G x Mass | |
| `grav` | 9.80665 m/s^2 | used for pressure conversion; gravity observables use GM/a^2 |

---

## 12. Tests and demos

```bash
python -m pytest tests -q       # consistency tests (needs data/LLN)
python src_demo/demo_slepian.py # Slepian truncation error over a regional mask
```

`tests/test_basis.py` checks the common interface: point-wise SH basis against grid synthesis, Slepian (all
functions retained) against SH, the Green's-function basis against the grid convolutions of `src_Forward`, agreement
of SH and point-mass forward models, LOS and load synthesis. `tests/test_consistency.py` checks the constants, Love-number frames, kernel shapes, north/east projection, the
spherical-harmonic schemes, Slepian equivalence for every observable, the gravity kernel against Wang et al. (2012)
and the above/below gravity relation. `src_demo/paper_JGR2025.py` reproduces the figures of Yang et al. (2025).

---

## 13. Known limitations

- The infinite LGF of Sect. 5 of the paper is not implemented (the paper's reference solution was computed with
  LoadVCD).
- Surface gravity is not available for D-LGF / CD-LGF (`DiskLoad.getGravity` raises).
- `LGF_Kummer.getGeoidheight` excludes degree 0 while the truncated and disk kernels include it.
- Tilt, strain and deflection of the vertical are not implemented.
- Spherical-harmonic analysis and the point/disk kernels store full Legendre arrays; memory grows as
  nlat x (L+1)^2 (about 375 MB at 0.5 deg, 3 GB at 0.25 deg).
- Earth-rotation excitation and geocentre motion belong to SAGEA-fluid; spectral unit conversions to SaGEA.

---

## 14. References

- Bevis, M., Melini, D., & Spada, G. (2016). On computing the geoelastic response to a disk load. *Geophysical
  Journal International*, 205(3), 1804-1812. https://doi.org/10.1093/gji/ggw115
- Blewitt, G. (2003). Self-consistency in reference frames, geocenter definition, and surface loading of the solid
  Earth. *Journal of Geophysical Research*, 108(B2), 2103. https://doi.org/10.1029/2002JB002082
- Farrell, W. E. (1972). Deformation of the Earth by surface loads. *Reviews of Geophysics*, 10(3), 761-797.
  https://doi.org/10.1029/RG010i003p00761
- Simons, F. J., Dahlen, F. A., & Wieczorek, M. A. (2006). Spatiospectral concentration on a sphere. *SIAM Review*,
  48(3), 504-536. https://doi.org/10.1137/S0036144504445765
- Sneeuw, N. (1994). Global spherical harmonic analysis by least-squares and numerical quadrature methods in
  historical perspective. *Geophysical Journal International*, 118(3), 707-716.
  https://doi.org/10.1111/j.1365-246X.1994.tb03995.x
- Wang, H., Xiang, L., Jia, L., Jiang, L., Wang, Z., Hu, B., & Gao, P. (2012). Load Love numbers and Green's functions
  for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0.
  *Computers & Geosciences*, 49, 190-199. https://doi.org/10.1016/j.cageo.2012.06.022
- Yang, F., Shen, Y., Zhang, W., & Forootan, E. (2025). On the choice of basis functions for modeling Earth's elastic
  deformation due to surface loading. *Journal of Geophysical Research: Solid Earth*, 130, e2025JB031662.
  https://doi.org/10.1029/2025JB031662
