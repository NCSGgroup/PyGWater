# PyGWater — installation guide

This guide covers the Python environment, the external data PyGWater needs, and how to run the code. Everything
needed for the environment lives in `installation/`:

| file | purpose |
|---|---|
| `environment.yml` | conda/mamba environment `pygwater` |
| `requirements.txt` | the same dependencies for `pip` |

---

## 1. Environment

PyGWater has no pinned binary stack (unlike PyGLDA v2, which is locked by WaterGAP 2.2e); any recent scientific
Python works. Tested with Python 3.13, numpy 2.5, scipy 1.18, h5py 3.16, netCDF4 1.7, geopandas 1.2 and shapely 2.2.

With conda / mamba:

```bash
git clone https://github.com/NCSGgroup/PyGWater.git
cd PyGWater
mamba env create -f installation/environment.yml
conda activate pygwater
```

With pip:

```bash
pip install -r installation/requirements.txt
```

| Dependency | Needed for |
|---|---|
| numpy, scipy, h5py, netCDF4 | everything |
| geopandas, shapely | regional masks from shapefiles (`src/mesh/shape_mask.py`) |
| pytest | the tests in `tests/` |
| pygmt | figures in `examples/` only (optional; can live in a separate environment) |

---

## 2. External data

All data PyGWater reads or writes live outside the code, under one folder, the **data root**, e.g.
`/media/user/My Book/Fan/PyGWater_external_data` (cf. `PyGLDA_v2_external_data` for PyGLDA v2). Nothing is written
inside the code folder. Set the data root once, in either of two ways:

1. a local config file (recommended; ignored by git): copy `installation/pygwater_paths.cfg.example` to
   `installation/pygwater_paths.cfg` and set

   ```
   data_root = /path/to/PyGWater_external_data
   ```

2. or the environment variable, which takes precedence over the file:

   ```bash
   export PYGWATER_DATA=/path/to/PyGWater_external_data      # e.g. in ~/.bashrc or the IDE run configuration
   ```

Without either, PyGWater stops with an error: there is no default data folder, and a data root inside the code
folder is refused. The code folder holds code only.

The data root is organised by role:

```
<data root>/
├── input/               external sources, only read by PyGWater
│   ├── LLN/             load Love numbers and Green's functions
│   ├── gnss/NGL/        DataHoldings.txt, steps.txt, tenv3/ (studies/GNSS_GRACE/download.py)
│   ├── land_mask/       land_mask_res0.1.h5 (removing ocean cells in pygwater.mesh)
│   ├── mascon/          GRACE mascon NetCDF (a symbolic link is fine) and its 2-degree extract mascon.txt
│   ├── shp/             region outlines: China/, NorthAmerica/, ...
│   ├── stations/        station lists and observations, e.g. fu2015_*_digitized.txt
│   └── reference/       solutions of other software for comparison, e.g. loadvcd.txt (LoadVCD)
├── cache/               built by PyGWater, regenerable and reused
│   ├── LGD_template*/   CD-LGF templates
│   ├── mask/            regional masks (.h5) from pygwater.mesh
│   └── shp/             generated grid shapefiles (EU/, exp/)
├── output/              results, one folder per study
│   ├── paper_JGR2025/   data of the figures of Yang et al. (2025): fig2/, fig8/, ..., BasisStudy/
│   ├── China_basis_comparison/
│   ├── inversion_China/
│   ├── inversion_PNW_Fu2015/
│   └── GNSS_GRACE/<region>/
└── temp/                scratch, safe to delete
```

In code: `FileTool.get_input_dir()`, `get_mascon_dir()`, `get_cache_dir()`, `get_output_dir()` and
`get_temp_dir()` (each takes an optional sub-path); the older demos use the variables `INPUT`, `CACHE`, `OUTPUT`
and `TEMP`. Never write a path to data directly in a script.

| File in `input/LLN/` | Content | Used by |
|---|---|---|
| `PREM-LLNs-complete.dat`, `iasp91-...`, `ak135-...`, `PREMhard-...`, `PREMsoft-...`, `iasp91hard-...`, `ak135hard-...` | load Love numbers, degrees 1 to 46 340, CE frame (Wang et al., 2012) | `LLN_Data.PREM`, ... |
| `Wang.dat` | load Love numbers to degree 13 000, CE frame | `LLN_Data.Wang` (default) |
| `REF_6371_loading_love_numbers.mat` | REF 6371 Love numbers, CM frame (Bevis et al., 2016) | `LLN_Data.REF` |
| `PREM-LGFs.dat` | tabulated PREM Green's functions (Wang et al., 2012) | one test, some demos |

These files are third-party data and are not covered by the PyGWater license; please cite their sources.

`installation/pygwater_paths.cfg` is ignored by git; it is the only local file in the code folder.

---

## 3. Installing and running PyGWater

The library is the folder `src/`; installing it in editable mode makes it importable as `pygwater` from anywhere,
while your edits in `src/` take effect immediately:

```bash
cd PyGWater
pip install -e .                    # once, in the pygwater environment
python -c "import pygwater; print(pygwater.__version__)"
```

Scripts and tests are then run from any folder:

```bash
python examples/slepian_truncation.py          # one feature per example
python studies/Fu2015_PNW/real_data.py         # a reproducible study
python -m pytest tests -q                      # all tests; those needing the Love numbers are skipped if input/LLN/ is missing
```

(The tests also run from a plain checkout without installation.) In PyCharm, select the `pygwater` interpreter after
the install; no source-root marking is needed.

See [user_guide.md](user_guide.md) for all options and usage.
