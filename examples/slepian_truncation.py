"""
Slepian basis demo (Sect. 8 of Yang et al., 2025): concentration spectrum (cf. Fig. 13) and truncation
(approximation) error of the Slepian forward model against spherical harmonics (cf. Fig. 14), for a regional mask
produced by pygwater.mesh.shape_mask.basinMask. Run from any folder:  python examples/slepian_truncation.py
"""
import os
import sys
from pathlib import Path


from pygwater.config import FileTool
INPUT = str(FileTool.get_input_dir()) + '/'     # <data root>/input/  (data root: installation/pygwater_paths.cfg)
CACHE = str(FileTool.get_cache_dir()) + '/'     # <data root>/cache/
OUTPUT = str(FileTool.get_output_dir()) + '/'   # <data root>/output/
TEMP = str(FileTool.get_temp_dir()) + '/'   # <data root>/temp/
import h5py
import numpy as np

import pygwater.utils.enums as Enums
from pygwater.utils.geomath import GeoMathKit
from pygwater.earth.love_numbers import LoveNumber
from pygwater.basis.slepian import SlepianBasis
from pygwater.basis.spherical_harmonic.data_class import GRID


def demo_region(mask_file=FileTool.get_cache_dir('mask/EU_res_0.5_crop.h5'),
                load_file=FileTool.get_mascon_dir('mascon.txt'), L=60):
    lln = LoveNumber().config(lmax=L, method=Enums.LLN_Data.PREM).get_Love_number()

    '''Slepian basis of the region'''
    sb = SlepianBasis.from_mask_h5(L, mask_file)
    print('Shannon number %.1f; beta>0.1: %d; beta>0.01: %d'
          % (sb.shannon, (sb.beta > 0.1).sum(), (sb.beta > 0.01).sum()))

    '''surface load: CSR-Mascon 2020-01 at 2 deg (column 3, cm), as in the paper'''
    lat, lon = GeoMathKit.global_equal_distance_grid(2)
    load = np.loadtxt(load_file)[:, 2].reshape(len(lat), len(lon)) * 0.01  # metre EWH
    shc = GRID(grid=load, lat=lat, lon=lon).analysis(lmax=L, from_type=Enums.PhysicalDimensions.EWH,
                                                     to_type=Enums.PhysicalDimensions.EWH, LLN=lln)
    c, s = shc.get_cs2d()

    '''evaluation points: cells inside the mask'''
    with h5py.File(mask_file, 'r') as f:
        inside = f['basin'][()] > 0
        plat, plon = f['lat'][()][inside], f['lon'][()][inside]

    u_sh = (sb.cs2d_to_vec(c, s) * sb._displacement_factor(lln)) @ sb.Ymatrix(plat, plon).T
    for Lc in (sb.ncoef, int((sb.beta > 0.01).sum()), int((sb.beta > 0.1).sum())):
        u = sb.truncate(Lc=Lc).forward(c, s, lln, plat, plon)
        r = (u - u_sh) * 1000
        print('Lc = %4d : Slepian - SH  rms %.3f mm, max %.3f mm' % (Lc, np.sqrt(np.mean(r ** 2)), np.abs(r).max()))

    return sb, sb.beta


if __name__ == '__main__':
    demo_region()
