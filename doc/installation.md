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
| geopandas, shapely | regional masks from shapefiles (`src_LoadMesh/ShapeAndMask.py`) |
| pytest | the tests in `tests/` |
| pygmt | figures in `src_demo/` only (optional; can live in a separate environment) |

---

## 2. External data

The data are kept outside the code, in a folder such as `PyGWater_external_data/` (cf. `PyGLDA_v2_external_data`
for PyGLDA v2). Point PyGWater to it with the environment variable `PYGWATER_DATA` (if unset, `<code root>/data/`
is used):

```bash
export PYGWATER_DATA=/path/to/PyGWater_external_data      # e.g. in ~/.bashrc or the IDE run configuration
```

```
<data root>/
├── LLN/          inputs: load Love numbers and Green's functions
├── land_mask/    inputs: land_mask_res0.1.h5 (only for removing ocean cells in src_LoadMesh)
├── res/          written by PyGWater: CD-LGF templates, masks, shapefiles, results of the demos
└── temp/         written by PyGWater: scratch files
```

In code: `FileTool.get_data_dir()`, `FileTool.get_res_dir()`, `FileTool.get_temp_dir()`; the demos use the
variables `RES` and `TEMP`.

| File in `LLN/` | Content | Used by |
|---|---|---|
| `PREM-LLNs-complete.dat`, `iasp91-...`, `ak135-...`, `PREMhard-...`, `PREMsoft-...`, `iasp91hard-...`, `ak135hard-...` | load Love numbers, degrees 1 to 46 340, CE frame (Wang et al., 2012) | `LLN_Data.PREM`, ... |
| `Wang.dat` | load Love numbers to degree 13 000, CE frame | `LLN_Data.Wang` (default) |
| `REF_6371_loading_love_numbers.mat` | REF 6371 Love numbers, CM frame (Bevis et al., 2016) | `LLN_Data.REF` |
| `PREM-LGFs.dat` | tabulated PREM Green's functions (Wang et al., 2012) | one test, some demos |

These files are third-party data and are not covered by the PyGWater license; please cite their sources.

If the default `<code root>/data/` is used, it is ignored by git (as are `res/` and `temp/` in the code root).

---

## 3. Running the code

The packages are imported as `src_<name>` from the code root (e.g. `from src_Basis.SphericalHarmonic.DataClass import GRID`). Scripts in
`src_demo/` and the tests put the code root on the path themselves and can be run from any folder:

```bash
python src_demo/demo_slepian.py      # Slepian truncation error over a regional mask
python -m pytest tests -q            # all tests; those needing the Love numbers are skipped if LLN/ is missing
```

For your own scripts, either run them from the code root or add it to `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/PyGWater:$PYTHONPATH
```

See [user_guide.md](user_guide.md) for all options and usage.
