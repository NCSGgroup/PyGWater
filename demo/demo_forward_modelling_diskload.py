import sys

sys.path.append('../')
sys.path.append('../../')

import pygmt
import numpy as np

from src.Auxiliary.EnumClasses import Displacement
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.DiskLoad import DiskLoad, EarthConstant
from src.ForwardModelling.Disk import GFA_displacement, grid2radius, GFA_regular_grid
import netCDF4 as nc
from pathlib import Path
from src.Auxiliary.GeoMathKit import GeoMathKit


def demo1_mascon():
    """Use CSR-Mascon to do the tests"""
    # lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    # gfa = GFA_displacement(lln=lln)

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 2 degree"""
    res = 1
    mascon_data = np.array(csr['lwe_thickness'])[180, res * 2::res * 4, res * 2::res * 4]  # unit: cm
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)

    lon = np.roll(lon, shift=int(180 / res))
    lon[lon < 0] += 360
    region = [min(lon), max(lon), min(lat), max(lat)]

    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    pygmt.config(FONT_ANNOT='11p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-500, 500], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=mascon_data.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf50', 'yf+lEWH (cm)'])

    fig.show()

    pass


def demo2_Green():
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

    rr = grid2radius(lat_center=lat, grid_size=res)

    grids = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
        'radius': rr[0].flatten(),
        'EWH': mascon_data * 0.01  # cm to meter
    }

    gfa.configure(grids=grids, cf=1000)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical)

    np.save('../temp/100_dis2.npy', dist)

    pass


def demo2_Green_fast():
    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
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

    rr = grid2radius(lat_center=lat, grid_size=res)

    grids = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
        'radius': rr[0].flatten(),
        'EWH': mascon_data * 0.01  # cm to meter
    }

    gfa.configure(grids=grids, cf=1000)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical)

    np.save('../temp/100_dis_fast.npy', dist)

    pass

def demo2_plot():
    dist = np.load('dis.npy') * 1000

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    # lon = np.roll(lon, shift=int(180 / res))
    # lon[lon < 0] += 360
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)
    dist = np.roll(dist.reshape((int(180 / res), int(360 / res))), shift=int(-180 / res), axis=1)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    pygmt.config(FONT_ANNOT='11p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-50, 50], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=dist.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lUplift (mm)'])

    fig.show()

    pass


def demo3_SH():
    from src.SphericalHarmonics.DataClass import GRID, Enums

    ''''''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, int(res * 2)::int(res * 4),
                  int(res * 2)::int(res * 4)] * 0.01  # unit: cm

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1)

    gd = GRID(grid=mascon_data, lat=lat, lon=lon)
    sh = gd.analysis(lmax=None, from_type=Enums.PhysicalDimensions.EWH, to_type=Enums.PhysicalDimensions.EWH)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.EWH,
                        to_type=Enums.PhysicalDimensions.VerticalDisplacement)

    disp = disp.value * 1000  # convert to mm

    np.save('../temp/sh2.npy', disp)

    pass


def compare_SH_Green():
    from src.SphericalHarmonics.DataClass import GRID, Enums

    res = 2
    num = 181
    sh = np.load('../temp/sh.npy').reshape((num, -1)).T
    green = np.load('../temp/100_dis.npy') * 1000

    index = 180

    '''std'''
    std_green = np.std(green, axis=1)
    std_sh = np.std(sh, axis=1)
    std_diff = np.std(green - sh, axis=1)
    vr = (1 - std_diff / std_green) * 100
    # vr = sh[:,index]-green[:,index]
    # vr = sh[:, index]
    vr = green[:, index]
    # vr = std_green
    # vr = std_sh

    # ref = np.load('../temp/dis.npy')*1000
    # vr = ref - green[:, 180]
    #
    # sh2 = np.load('../temp/sh2.npy').flatten()
    # vr = green[:, 180]
    # convert to mm
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    sn = np.array([lat.flatten(), lon.flatten(), vr.flatten()])

    np.savetxt('../temp/GreenFunction.txt', sn.T)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    pygmt.config(FONT_ANNOT='11p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-50, 50], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=vr.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])
    #
    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lstd (mm)'])

    fig.show()

    pass


def compare_yang_shen():
    num = 181
    sh = np.load('../temp/sh.npy').reshape((num, -1)).T
    sh = np.load('../temp/sh2.npy')
    # green = np.load('../temp/100_dis2.npy') * 1000
    # green= np.load('../temp/100_dis_fast.npy') * 1000
    green = np.load('../temp/cluster_disk.npy') * 1000
    # green = np.load('../temp/pointload_dist_fast.npy') * 1000

    index = 180

    '''std'''
    # std_green = np.std(green, axis=1)
    # std_sh = np.std(sh, axis=1)
    # std_diff = np.std(green - sh, axis=1)
    # vr = (1 - std_diff / std_green) * 100
    # vr = sh[:,index]-green[:,index]
    # vr = sh[:, index]
    green = green[:, 0]
    # sh = sh[:, index]

    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    pygmt.config(FONT_ANNOT='11p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=shen.flatten() - sh.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Shen minus SH'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])
    #
    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lstd (mm)'])

    '''SH- Fan'''
    fig.shift_origin(xshift='13c')

    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=green.flatten() - sh.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) Green of Fan minus SH'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    #
    '''Green:  Fan - Shen'''
    fig.shift_origin(xshift='-13c', yshift='-9c')

    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=-green.flatten() + shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c) Green: Shen minus Fan'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    fig.show()

    pass


def cal_SH_loss():
    from src.SphericalHarmonics.DataClass import GRID, Enums
    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)
    # shen = np.load('../temp/100_dis2.npy') * 1000

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    ''''''
    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    pygmt.config(FONT_ANNOT='11p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-50, 50], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Green displacement'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    '''grid --> SH --> grid'''
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    gd = GRID(grid=shen.reshape((90, 180)), lat=lat, lon=lon)
    sh = gd.analysis(lmax=None, from_type=Enums.PhysicalDimensions.Dimensionless, to_type=Enums.PhysicalDimensions.Dimensionless)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.Dimensionless,
                        to_type=Enums.PhysicalDimensions.Dimensionless).value

    lon, lat = np.meshgrid(lon, lat)
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten(),
                          spacing=(res, res), region=region)

    fig.shift_origin(xshift='13c')
    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) SH displacement'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    ''''''
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=-disp.flatten()+shen.flatten(),
                          spacing=(res, res), region=region)

    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    fig.shift_origin(xshift='-13c',yshift='-9c')
    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c) Green minus SH displacement'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    fig.show()

    pass


if __name__ == '__main__':
    # demo1()
    # demo2_plot()
    # demo3()
    # demo3_SH()
    # compare_SH_Green()
    # demo2_Green()
    compare_yang_shen()
    # cal_SH_loss()
    # demo2_Green_fast()
