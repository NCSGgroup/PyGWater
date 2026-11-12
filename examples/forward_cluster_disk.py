import os
import sys
from pathlib import Path


from pygwater.config import FileTool
INPUT = str(FileTool.get_input_dir()) + '/'     # <data root>/input/  (data root: installation/pygwater_paths.cfg)
CACHE = str(FileTool.get_cache_dir()) + '/'     # <data root>/cache/
OUTPUT = str(FileTool.get_output_dir()) + '/'   # <data root>/output/
TEMP = str(FileTool.get_temp_dir()) + '/'   # <data root>/temp/


import pygmt
import numpy as np

from pygwater.utils.enums import Displacement
from pygwater.earth.love_numbers import LoveNumber, LLN_Data, LLN_variable, Frame
from pygwater.basis.green.disk_load import DiskLoad, EarthConstant
from pygwater.forward.disk import GFA_displacement, grid2radius, GFA_regular_grid
import netCDF4 as nc
from pathlib import Path
from pygwater.utils.geomath import GeoMathKit
from pygwater.forward.disk_cluster import GFA_displacement_from_template


def demo1_convolve():
    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_displacement_from_template(template_dir=CACHE + 'LGD_template/')

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = str(FileTool.get_mascon_dir())
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]
    # mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]# unit: cm

    mascon_data = np.roll(mascon_data, shift=int(180 / res), axis=-1)

    # mascon_data = mascon_data.flatten()
    mascon_data = mascon_data.reshape((1, -1)).T
    # mascon_data = mascon_data.reshape((index+1, -1)).T

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    a = np.arange(90)
    b = np.ones(180)

    b, a = np.meshgrid(b, a)

    # sn = np.array([lat.flatten(),lon.flatten(),mascon_data.flatten()])
    #
    # np.savetxt(INPUT + 'mascon/mascon.txt', sn.T)

    grids = {
        'EWH': mascon_data * 0.01,  # cm to meter
        'Template': a.flatten()
    }

    dist = gfa.evaluation_regular_grid(grids=grids, dim=90)

    np.save(TEMP + 'cluster_disk.npy', dist)

    pass


if __name__ == '__main__':
    demo1_convolve()
