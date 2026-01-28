import sys

sys.path.append('../')
sys.path.append('../../')

import numpy as np
from src.Auxiliary.EnumClasses import Displacement
from src.ForwardModelling.PointMass_truncated import GFA_displacement, LoveNumber, LLN_Data, Grids_generation, \
    GFA_regular_grid
import netCDF4 as nc
from pathlib import Path
import pygmt
from src.Auxiliary.GeoMathKit import GeoMathKit
from src.Auxiliary.constants import EarthConstant


def demo_Green_genenral():
    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_displacement(lln=lln)

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]
    # mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]# unit: cm

    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1)

    # mascon_data = mascon_data.flatten()
    mascon_data = mascon_data.reshape((1, -1)).T
    # mascon_data = mascon_data.reshape((index+1, -1)).T

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    # sn = np.array([lat.flatten(),lon.flatten(),mascon_data.flatten()])
    #
    # np.savetxt('../temp/mascon.txt', sn.T)

    grids = Grids_generation.Equal_angular_distance(resolution=res)

    grids['EWH'] = mascon_data * 0.01 *EarthConstant.rhow # cm to meter

    gfa.configure(grids=grids)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical)

    np.save('../temp/pointload_dist.npy', dist)

    pass


def demo_Green_fast():
    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=90, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_regular_grid(lln=lln)

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]
    # mascon_data = np.array(csr['lwe_thickness'])[index, res * 2::res * 4, res * 2::res * 4]# unit: cm

    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1)

    # mascon_data = mascon_data.flatten()
    mascon_data = mascon_data.reshape((1, -1)).T
    # mascon_data = mascon_data.reshape((index+1, -1)).T

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    # sn = np.array([lat.flatten(),lon.flatten(),mascon_data.flatten()])
    #
    # np.savetxt('../temp/mascon.txt', sn.T)

    grids = Grids_generation.Equal_angular_distance(resolution=res)

    grids['EWH'] = mascon_data * 0.01* EarthConstant.rhow # cm to meter

    gfa.configure(grids=grids)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical, resolution=res)

    np.save('../temp/pointload_dist_fast.npy', dist)

    pass


if __name__ == '__main__':
    demo_Green_fast()
    # demo_Green_genenral()
