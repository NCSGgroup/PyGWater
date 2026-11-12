# PyGWater

**Surface loading and GNSS inversion of terrestrial water storage, with a choice of basis functions.**

PyGWater models the elastic deformation of the Earth under surface water loads (GNSS displacement, geoid, gravity,
InSAR line of sight) and inverts geodetic displacements for the water storage that caused them. Spherical harmonics,
point-mass and disk-load Green's functions (D-LGF, CD-LGF) and Slepian functions share one interface, so any
computation can be repeated with another basis (Yang et al., 2025).

## Installation

```bash
git clone https://github.com/NCSGgroup/PyGWater.git
cd PyGWater
mamba env create -f installation/environment.yml
conda activate pygwater
pip install -e .
```

Data are kept outside the code, in one data root. Copy `installation/pygwater_paths.cfg.example` to
`installation/pygwater_paths.cfg` and set `data_root`; the Love numbers go in `<data root>/input/LLN/`.
See [doc/installation.md](doc/installation.md).

## Quick start

```python
import numpy as np
import pygwater as pgw

lln = pgw.LoveNumber().config(lmax=13000, method=pgw.LLN_Data.Wang).get_Love_number()

# 0.25-degree disk cells and 60 GNSS stations
r = 0.25
cells = pgw.GreenFunctionBasis.from_grid(np.arange(44 + r / 2, 48, r), np.arange(4 + r / 2, 9, r), r, kernel='disk')
rng = np.random.default_rng(1)
slat, slon = rng.uniform(44.2, 47.8, 60), rng.uniform(4.2, 8.8, 60)
A = cells.design_matrix(slat, slon, lln, pgw.Displacement.Vertical)

# forward: vertical displacement from a water load; inverse: the load back from noisy data
x_true = 0.2 * np.exp(-((cells.lat - 46) ** 2 + (cells.lon - 6.5) ** 2))     # [m EWH]
d = A @ x_true + rng.normal(0, 1e-3, slat.size)                              # [m]

solver = pgw.RegularizedLeastSquares(A, 1e-3, pgw.regularization.curvature(cells.lat, cells.lon, r))
alpha, _ = solver.choose_alpha(d, alphas=np.logspace(-6, 5, 45), method='gcv')
x = solver.solve_bounded(d, alpha, lb=0.0).x
```

More in [examples/](examples) and the reproducible [studies/](studies); all options in the
[user guide](doc/user_guide.md). Run the tests with `python -m pytest -q`.

## Citation

> Yang, F., Shen, Y., Zhang, W., & Forootan, E. (2025). On the choice of basis functions for modeling Earth's elastic
> deformation due to surface loading. *Journal of Geophysical Research: Solid Earth*, 130, e2025JB031662.
> https://doi.org/10.1029/2025JB031662

See also [CITATION.cff](CITATION.cff). PyGWater belongs to the GravityPython toolboxes of the NCSG group, with
[PyHawk](https://github.com/NCSGgroup/PyHawk), [SaGEA](https://github.com/NCSGgroup/SaGEA),
[PyGLDA](https://github.com/AAUGeodesyGroup/PyGLDA) and [SAGEA-fluid](https://github.com/NCSGgroup/SAGEA-fluid).

## License

[MIT](LICENSE). The load Love numbers (Wang et al., 2012; Bevis et al., 2016) are third-party data.
Contact: Fan Yang, Aalborg University — [fany@plan.aau.dk](mailto:fany@plan.aau.dk)
