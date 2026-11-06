import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE_ROOT = HERE.parent
sys.path.insert(0, str(CODE_ROOT))   # run from any folder: the src_* packages import from the code root
os.chdir(HERE)                       # relative paths in this script ('../res', '../temp') refer to src_demo/


import pygmt
import numpy as np

from src_auxiliary.EnumClasses import Displacement
from src_EarthModel.LoveNumber import LoveNumber, LLN_Data, LLN_variable, Frame
from src_Basis.GreenFunction.DiskLoad import DiskLoad, EarthConstant
from src_Forward.Disk import GFA_displacement, grid2radius, GFA_regular_grid
import netCDF4 as nc
from pathlib import Path
from src_auxiliary.GeoMathKit import GeoMathKit
from src_Forward.DiskCluster import GFA_displacement_from_template


def demo1_convolve():
    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_displacement_from_template(template_dir='../res/LGD_template/')

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
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
    # np.savetxt('../temp/mascon.txt', sn.T)

    grids = {
        'EWH': mascon_data * 0.01,  # cm to meter
        'Template': a.flatten()
    }

    dist = gfa.evaluation_regular_grid(grids=grids, dim=90)

    np.save('../temp/cluster_disk.npy', dist)

    pass


if __name__ == '__main__':
    demo1_convolve()
