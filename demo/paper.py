import sys

sys.path.append('../')
sys.path.append('../../')

import pygmt
import numpy as np

from src.Auxiliary.EnumClasses import Displacement
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
import netCDF4 as nc
from pathlib import Path
from src.Auxiliary.GeoMathKit import GeoMathKit
from src.Auxiliary.constants import EarthConstant


def fig1():
    from src.SphericalHarmonics.DataClass import GRID, Enums

    ''''''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, int(res * 2)::int(res * 4),
                  int(res * 2)::int(res * 4)]  # unit: cm

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1)

    gd = GRID(grid=mascon_data, lat=lat, lon=lon)
    sh = gd.analysis(lmax=90, from_type=Enums.PhysicalDimensions.EWH, to_type=Enums.PhysicalDimensions.EWH)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.EWH,
                        to_type=Enums.PhysicalDimensions.EWH)

    disp = disp.value  # convert to mm

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-60, 60], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=mascon_data.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Original'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa20f', 'yf+lEWH (cm)'])
    #
    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lstd (mm)'])

    '''(b)'''
    fig.shift_origin(xshift='13c')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten() - mascon_data.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) Signal loss and Gibbs effect'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa20f', 'yf+lEWH (cm)'])

    '''(c)'''
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    gd = GRID(grid=mascon_data, lat=lat, lon=lon)
    sh = gd.analysis(lmax=120, from_type=Enums.PhysicalDimensions.EWH, to_type=Enums.PhysicalDimensions.EWH)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.EWH,
                        to_type=Enums.PhysicalDimensions.EWH)

    disp = disp.value  # convert to mm

    fig.shift_origin(xshift='-13c', yshift='-9c')

    lon, lat = np.meshgrid(lon, lat)
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten() - mascon_data.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c): Aliasing error (N=120)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa20f', 'yf+lEWH (cm)'])

    '''(d)'''
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    gd = GRID(grid=mascon_data, lat=lat, lon=lon)
    sh = gd.analysis(lmax=180, from_type=Enums.PhysicalDimensions.EWH, to_type=Enums.PhysicalDimensions.EWH)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.EWH,
                        to_type=Enums.PhysicalDimensions.EWH)

    disp = disp.value  # convert to mm

    fig.shift_origin(xshift='13c')

    lon, lat = np.meshgrid(lon, lat)
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten() - mascon_data.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(d): Aliasing error (N=180)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa20f', 'yf+lEWH (cm)'])

    fig.show()

    pass


def fig2_prepare():
    from src.ForwardModelling.PointMass_truncated import GFA_regular_grid, Grids_generation
    from src.SphericalHarmonics.DataClass import GRID, Enums

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 2 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, int(res * 2)::int(res * 4),
                  int(res * 2)::int(res * 4)]  # unit: cm

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1) * 0.01  # convert to meter

    '''SH solution'''
    gd = GRID(grid=mascon_data, lat=lat, lon=lon)
    sh = gd.analysis(lmax=90, from_type=Enums.PhysicalDimensions.EWH, to_type=Enums.PhysicalDimensions.EWH)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.EWH,
                        to_type=Enums.PhysicalDimensions.VerticalDisplacement)
    disp = disp.value
    np.save('../res/fig2/SH.npy', disp)

    '''Green function: N=180'''
    mascon_data = mascon_data.reshape((1, -1)).T
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    grids = Grids_generation.Equal_angular_distance(resolution=res)

    grids['EWH'] = mascon_data * EarthConstant.rhow

    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=360, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_regular_grid(lln=lln)
    gfa.configure(grids=grids)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical, resolution=res)

    np.save('../res/fig2/Green360.npy', dist)

    pass


def fig2():
    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    '''SH'''
    SH1 = np.load('../res/fig2/SH.npy') * 1000

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-40, 40], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=SH1.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Spherical Harmonic (N=90)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa20f', 'yf+luplift (mm)'])

    '''Green N=90'''
    Green = np.load('../res/fig2/Green90.npy') * 1000

    SH = Green.flatten() - SH1.flatten()
    fig.shift_origin(xshift='13c')

    pygmt.makecpt(cmap='haxby', series=[-5, 5], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=SH.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) Finite Green function (N=90) minus (a)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa1f', 'yf+luplift (mm)'])

    '''Green N=180'''
    Green = np.load('../res/fig2/Green180.npy') * 1000

    SH = Green.flatten() - SH1.flatten()
    fig.shift_origin(xshift='-13c', yshift='-9c')

    pygmt.makecpt(cmap='haxby', series=[-20, 20], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=SH.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c) same as (b) but with N=180'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa5f', 'yf+luplift (mm)'])

    '''Green N=180'''
    Green = np.load('../res/fig2/Green360.npy') * 1000

    SH = Green.flatten() - SH1.flatten()
    fig.shift_origin(xshift='13c')

    pygmt.makecpt(cmap='haxby', series=[-20, 20], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=SH.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(d) same as (b) but with N=360'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xa5f', 'yf+luplift (mm)'])

    fig.show()
    pass


def fig3():
    from src.GreenFunction.PointLoad import LGF_Kummer, LGF_truncation

    '''(a)'''
    from src.Auxiliary.FileTool import FileTool

    file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    theta = np.loadtxt(file, skiprows=1, usecols=0)
    theta[-1] -= 0.01
    data = {}
    lmax = 38000
    # list_LLN = [LLN_Data.PREM,LLN_Data.ak135,LLN_Data.iasp91,LLN_Data.REF]
    list_LLN = [LLN_Data.PREM, LLN_Data.PREMsoft, LLN_Data.PREMhard]
    for item in list_LLN:
        lln = LoveNumber().config(lmax=lmax, method=item).get_Love_number().convert(target=Frame.CE)
        load = LGF_Kummer(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        v = load.getHorizental() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        g = load.getGeoidheight() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = [u, v, g]

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='9p', FONT_TITLE='14p', MAP_TITLE_OFFSET=0.1)
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne+t(a) Infinite LGF (N=40000)", "xf3", "ya100f50"],
                region=[0.00001, 179, -400, 100],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item][0]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label=item.name,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 6,
            pen="2p,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)
    # fig.text(
    #     text="(a) Infinite LGF up to N=40000",
    #     font='13p',
    #     position="TL",  # Top Left
    #     justify="TL",  # Top Left
    #     offset="0.1c/-0.1c",
    # )

    """(b)"""
    fig.shift_origin(xshift='9.2c')
    for item in list_LLN:
        lln = LoveNumber().config(lmax=lmax, method=item).get_Love_number().convert(target=Frame.CE)
        load = LGF_truncation(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        # v = load.getHorizental() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        # g = load.getGeoidheight() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = [u, v, g]

    fig.basemap(frame=["WSne+t(b) Finite LGF (N=40000)", "xf3", "ya20f10"],
                region=[0.00001, 179, -100, 50],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item][0]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label=item.name,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 6,
            pen="2p,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)

    """(c)"""
    lmax = 38000
    list_LLN = [100, 1000, 40000]
    data = {}
    for item in list_LLN:
        lln = LoveNumber().config(lmax=item, method=LLN_Data.PREM).get_Love_number()
        load = LGF_Kummer(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = u

    fig.shift_origin(xshift='-9.2c', yshift='-9c')

    fig.basemap(frame=["WSne+t(c) Infinite LGF (PREM)", "xa1f3+l@[\\alpha@[", "ya10f5"],
                region=[0.00001, 179, -50, 10],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label='N=' + str(item),
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 6,
            pen="2p,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)

    """(d)"""
    lmax = 38000
    list_LLN = [100, 1000, 10000, 40000]
    data = {}
    for item in list_LLN:
        lln = LoveNumber().config(lmax=item, method=LLN_Data.PREM).get_Love_number()
        load = LGF_truncation(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = u

    fig.shift_origin(xshift='9.2c')

    fig.basemap(frame=["WSne+t(d) Finite LGF (PREM)", "xa1f3+l@[\\alpha@[", "ya10f5"],
                region=[0.00001, 179, -50, 10],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label='N=' + str(item),
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 6,
            pen="2p,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)

    fig.show()

    pass


def fig7():
    from src.SphericalHarmonics.DataClass import GRID, Enums

    sh = np.load('../res/fig2/SH.npy') * 1000
    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)
    # shen = np.load('../temp/cluster_disk.npy') * 1000

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    """(a)"""
    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-40, 40], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Infinite LGF (N=10000)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    """(b)"""
    fig.shift_origin(xshift='13c')
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=sh.flatten() - shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b): Spherical harmonic (N=90) minus (a)'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    """(c)"""
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    gd = GRID(grid=shen.reshape((90, 180)), lat=lat, lon=lon)
    sh = gd.analysis(lmax=None, from_type=Enums.PhysicalDimensions.Dimensionless,
                     to_type=Enums.PhysicalDimensions.Dimensionless)
    disp = sh.synthesis(grid_space=res, from_type=Enums.PhysicalDimensions.Dimensionless,
                        to_type=Enums.PhysicalDimensions.Dimensionless).value

    fig.shift_origin(xshift='-13c', yshift='-9c')
    pygmt.makecpt(cmap='haxby', series=[-40, 40], background='w')
    lon, lat = np.meshgrid(lon, lat)
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c): Harmonic analysis of (a) up to N=90'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    """(d)"""
    fig.shift_origin(xshift='13c')
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=disp.flatten() - shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30', 'yf15'] + ['+t(d): (c) minus (a)'],
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


def fig8_prepare():
    from src.ForwardModelling.PointMass_truncated import GFA_displacement, LoveNumber, LLN_Data, Grids_generation, \
        GFA_regular_grid

    """Use CSR-Mascon to do the tests"""
    lln = LoveNumber().config(lmax=1000, method=LLN_Data.PREM).get_Love_number()
    gfa = GFA_regular_grid(lln=lln)

    '''===============handling mascon data========================'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 1
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

    grids['EWH'] = mascon_data * 0.01 * EarthConstant.rhow  # cm to meter

    gfa.configure(grids=grids)

    point = {
        'lat': lat.flatten(),
        'lon': lon.flatten(),
    }

    dist = gfa.evaluation(points=point, variable=Displacement.Vertical, resolution=res)

    # np.save('../res/fig8/degree2.npy', dist)
    np.save('../res/fig8/degree1.npy', dist)

    pass


def fig8():
    # shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)
    # shen = np.load('../res/fig8/degree2.npy') * 1000

    """(a)"""
    shen = np.load('../res/fig8/degree1.npy') * 1000
    res = 1
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p')
    pygmt.makecpt(cmap='haxby', series=[-40000, 40000], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) Resolution at 1 degree'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    """(b)"""
    shen = np.load('../res/fig8/degree2.npy') * 1000
    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig.shift_origin(xshift='13c')
    pygmt.makecpt(cmap='haxby', series=[-40000, 40000], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) Resolution at 2 degree'],
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


def fig9():
    from src.GreenFunction.PointLoad import LGF_Kummer, LGF_truncation

    '''(a)'''
    from src.Auxiliary.FileTool import FileTool

    file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    theta = np.loadtxt(file, skiprows=1, usecols=0)
    theta[-1] -= 0.01

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='9p', FONT_TITLE='14p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    lmax = 38000
    list_LLN = [10000]
    data = {}
    for item in list_LLN:
        lln = LoveNumber().config(lmax=item, method=LLN_Data.PREM).get_Love_number()
        load = LGF_Kummer(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = u

    fig.basemap(frame=["WSne+t(a) LGF versus D-LGF, load=1 kg", "xa1f3+l@[\\alpha@[", "ya10f5"],
                region=[0.00001, 179, -50, 10],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label='Infinite LGF (N=10000)',
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 6,
            pen="2p,+z,-",
        )

    list_LLN = [10000]
    data = {}
    for item in list_LLN:
        lln = LoveNumber().config(lmax=item, method=LLN_Data.PREM).get_Love_number()
        load = LGF_truncation(lln=lln).configure(residence=theta)
        u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[item] = u

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label='Finite LGF (N=10000)',
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue + 8,
            pen="2p,+z,-",
        )

    from src.GreenFunction.DiskLoad import DiskLoad

    lmax = 10000
    lln = LoveNumber().config(lmax=lmax, method=LLN_Data.PREM).get_Love_number()
    radius = 0.05
    area = 2 * np.pi * (1 - np.cos(np.deg2rad(radius))) * EarthConstant.radiusm ** 2

    M = 1
    thickness = M / EarthConstant.rhow / area
    load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    u2 = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    xpoints = theta
    ypoints = u2
    fig.plot(
        # Set the figures frame and color as well as
        # annotations, ticks, and gridlines
        x=xpoints,
        y=ypoints,
        label='D-LGF with radius = 0.05\\260 (N=10000)',
        # Set the line thickness to "1p", the color to "blue",
        # and the style to "-", i.e. "dashed"
        cmap=True,
        zvalue=1,
        pen="2p,--",
    )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jTL", box=True)

    """(b)"""
    fig.shift_origin(xshift='9.2c')

    fig.basemap(frame=["WSne+t(b) D-LGF, load=1 kg", "xa1f3+l@[\\alpha@[", "ya10f5"],
                region=[0.00001, 179, -50, 10],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    list_radius = [3, 2, 1, 0.5, 0.1, 0.05, 0.001]
    data = {}
    for radius in list_radius:
        area = 2 * np.pi * (1 - np.cos(np.deg2rad(radius))) * EarthConstant.radiusm ** 2
        M = 1
        thickness = M / EarthConstant.rhow / area
        load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
        u2 = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
        data[radius] = u2

    for zvalue, item in enumerate(list_radius):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            label='Radius=%s\\260' % str(item),
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)

    fig.show()

    pass


def fig10():
    sh = np.load('../temp/sh2.npy')
    green = np.load('../temp/100_dis_fast.npy') * 1000
    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p', COLOR_NAN='grey')
    pygmt.makecpt(cmap='haxby', series=[-40, 40], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=green.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) D-LGF'],
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

    '''(b)'''
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    fig.shift_origin(xshift='-6.5c', yshift='-9c')

    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=-green.flatten() + sh.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) D-LGF minus Spherical harmonic'],
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
    '''(c)'''
    fig.shift_origin(xshift='13c')

    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=-green.flatten() + shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(c) D-LGF minus infinite LGF'],
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


def fig11_prepare():
    from src.GreenFunction.DiskLoad import DiskLoad

    lmax_list = np.arange(100, 40100, 100)
    radius = 0.1
    M = 1

    radius_list = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]

    for radius in radius_list:
        print(radius)
        # area = 2 * np.pi * (1 - np.cos(np.deg2rad(radius))) * EarthConstant.radiusm ** 2
        # thickness = M / EarthConstant.rhow / area

        thickness = 1

        K = 1.5
        theta = K * radius
        data = []
        for lmax in lmax_list:
            lln = LoveNumber().config(lmax=lmax, method=LLN_Data.REF).get_Love_number()
            load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
            u = load.getVertical()
            data.append(u)

        np.save('../res/fig11/far_%s.npy' % radius, np.array(data))

        K = 0.1
        theta = K * radius
        data = []
        for lmax in lmax_list:
            lln = LoveNumber().config(lmax=lmax, method=LLN_Data.REF).get_Love_number()
            load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
            u = load.getVertical()
            data.append(u)

        np.save('../res/fig11/near_%s.npy' % radius, np.array(data))

    pass


def fig11():
    from src.GreenFunction.DiskLoad import DiskLoad

    lmax_list = np.arange(100, 40100, 100)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=0.1)
    # pygmt.config(FONT_ANNOT='13p')
    pygmt.config(FONT_LABEL='11p', FONT_ANNOT_PRIMARY='8p', FONT_TITLE='14p', FONT_ANNOT='10p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    radius_list = [0.05, 0.1, 0.5, 1.0, 2.0, 3.0]

    n = -1
    label = {
        0: '(a) radius=0.05\\260',
        1: '(b)',
        2: '(c) radius=0.1\\260',
        3: '(d)',
        4: '(e) radius=0.5\\260',
        5: '(f)',
        6: '(g) radius=1.0\\260',
        7: '(h)',
        8: '(i) radius=2.0\\260',
        9: '(j)',
        10: '(k) radius=3.0\\260',
        11: '(l)',
    }
    for radius in radius_list:
        n += 1
        # radius =0.1
        ypoints = np.load('../res/fig11/near_%s.npy' % radius) * 1000
        xpoints = lmax_list

        if radius == 3.0:
            fig.basemap(frame=["WSne",
                               "xa1pf3+l N@-max@-", "yaf+l mm"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")
        elif radius == 0.05:
            fig.basemap(frame=["WSne+t Near-field (at 0.1\\327radius)",
                               "xf3", "yaf+l mm"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")
        else:
            fig.basemap(frame=["WSne",
                               "xf3", "yaf+l mm"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")

        zvalue = 9

        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )
        fig.text(
            text=label[n],
            font='10p',
            position="TL",  # Top Left
            justify="TL",  # Top Left
            offset="0.1c/-0.1c",
        )

        n += 1
        nROT = 360 / radius
        fig.plot(x=[nROT, nROT], y=[-200, 200], pen='1p,Grey,--')

        fig.plot(x=[2 * nROT, 2 * nROT], y=[-200, 200], pen='1p,Red,--')

        fig.plot(x=[10000, 10000], y=[-200, 200], pen='1p,Black,--')

        fig.shift_origin(xshift='9.2c')
        ypoints = np.load('../res/fig11/far_%s.npy' % radius) * 1000
        xpoints = lmax_list

        if radius == 3.0:
            fig.basemap(frame=["WSne",
                               "xa1pf3+l N@-max@-", "yaf"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")
        elif radius == 0.05:
            fig.basemap(frame=["WSne+t Far-field (at 1.5\\327radius)",
                               "xf3", "yaf"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")

        else:
            fig.basemap(frame=["WSne",
                               "xf3", "yaf"], region=[100, 40000, min(ypoints), max(ypoints)],
                        # Set a logarithmic transformation on the x-axis
                        projection="X8cl/2c")

        zvalue = 6

        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )

        nROT = 360 / radius
        fig.plot(x=[nROT, nROT], y=[-200, 200], pen='1p,Grey,--')

        fig.plot(x=[2 * nROT, 2 * nROT], y=[-200, 200], pen='1p,Red,--')

        fig.plot(x=[10000, 10000], y=[-200, 200], pen='1p,Black,--')

        fig.text(
            text=label[n],
            font='10p',
            position="TL",  # Top Left
            justify="TL",  # Top Left
            offset="0.1c/-0.1c",
        )

        fig.shift_origin(xshift='-9.2c', yshift='-2.5c')

    fig.show()

    pass


def fig14():
    from src.GreenFunction.DiskLoad import DiskLoad
    from src.ForwardModelling.Disk import grid2radius_type2

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='14p', FONT_ANNOT_PRIMARY='12p', FONT_TITLE='9p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    """CD-LGF"""
    lat0, lon0 = 1, -179
    lat1 = lat0 * np.ones(1000)
    lon1 = lon0 + np.linspace(start=0.0001, stop=10, num=1000)
    X_angle1 = GeoMathKit.angular_distance(lat0, lon0, lat1, lon1)
    CD1 = np.load('../res/fig14/LGD_fig14_res_0.05.npy') * 1000
    CD2 = np.load('../res/fig14/LGD_fig14_res_0.1.npy') * 1000
    CD3 = np.load('../res/fig14/LGD_fig14_res_0.2.npy') * 1000
    CD4 = np.load('../res/fig14/LGD_fig14_res_0.5.npy') * 1000

    """D-LGF"""
    LLN_list = [LLN_Data.PREM, LLN_Data.PREMhard]
    U, V, G = {}, {}, {}
    for llndata in LLN_list:
        lln = LoveNumber().config(lmax=10000, method=llndata).get_Love_number()

        # file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
        # theta = np.loadtxt(file, skiprows=1, usecols=0)
        # theta = 0.001100220044009
        radius = grid2radius_type2(lat_center=1, grid_size=2)[0]
        thickness = 1
        theta = X_angle1
        load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
        u = load.getVertical()
        U[llndata] = u * 1000

    ypoints = {
        'D-LGF (2.0\\260)': U[LLN_Data.PREM],
        'CD-LGF (0.05\\260)': CD1,
        'CD-LGF (0.1\\260)': CD2,
        'CD-LGF (0.2\\260)': CD3,
        'CD-LGF (0.5\\260)': CD4,
    }

    fig.basemap(frame=["WSne",
                       "xaf+u\\260+l@[\\alpha@[", "yaf+l mm"], region=[0, 2, -15, -1],
                # Set a logarithmic transformation on the x-axis
                projection="X7.5c/8c")

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=X_angle1,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend(position="jBR", box=True)

    fig.text(
        text="(a)",
        font='17p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    """(b)"""
    ypoints = {
        'D-LGF (PREM versus PREMhard)': U[LLN_Data.PREM] - U[LLN_Data.PREMhard],
        'CD-LGF (0.05\\260) versus D-LGF (2.0\\260)': CD1 - U[LLN_Data.PREM],
        'CD-LGF (0.1\\260) versus CD-LGF (0.05\\260)': CD2 - CD1,
        'CD-LGF (0.2\\260) versus CD-LGF (0.05\\260)': CD3 - CD1,
        # 'CD-LGF (0.5\\260)': CD4,
    }
    fig.shift_origin(xshift='8.5c')
    fig.basemap(frame=["WSne",
                       "xaf+u\\260+l@[\\alpha@[", "yaf"], region=[0, 2, -2, 3],
                # Set a logarithmic transformation on the x-axis
                projection="X7.5c/8c")

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=X_angle1,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend(position="jBR", box=True)

    fig.text(
        text="(b)",
        font='17p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    fig.show()
    pass


def fig15():
    """(a)"""
    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)
    green2 = np.load('../res/fig15/cluster_disk_0.2.npy') * 1000
    green1 = np.load('../res/fig15/cluster_disk_0.1.npy') * 1000

    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.1)
    pygmt.config(FONT_ANNOT='13p')
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=green2.flatten() - shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(a) CD-LGF (0.2\\260) minus infinite LGF'],
        dpi=100,
        projection='Q12c',
        region=region,
        # interpolation='n'
        # cmap = 'batlow',
    )

    fig.coast(shorelines="1/0.5p", region=region, projection="Q12c")

    # fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+lVR (%)'])

    fig.colorbar(position="JBC+w8.5c/0.25c+h+o-1.2c/0.6c", frame=['xaf', 'yf+luplift (mm)'])

    """(b)"""
    res = 2
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]
    lon, lat = np.meshgrid(lon, lat)

    fig.shift_origin(xshift='13c')
    pygmt.makecpt(cmap='haxby', series=[-10, 10], background='w')
    grid1 = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=green1.flatten() - shen.flatten(),
                          spacing=(res, res), region=region)

    fig.grdimage(
        grid=grid1,
        cmap=True,
        # frame=['WSne']+ ['xf10', 'yf1'] + ['+t(d) Weak > 0.3'],
        frame=['xf30',
               'yf15'] + ['+t(b) CD-LGF (0.1\\260) minus infinite LGF'],
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


def save_data_to_repository():
    from src.SphericalHarmonics.DataClass import GRID, Enums
    import h5py

    '''create a data'''
    bsd = h5py.File(name='/media/user/My Book/Fan/BasisStudy/res.hdf5', mode='w')

    '''input CSR-Mascon data'''
    filename = 'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc'

    dir_in = '/media/user/Backup Plus/GRACE/Mascon/CSR'
    csr = nc.Dataset(Path(dir_in) / filename)

    """upscale to 1 degree"""
    res = 2
    index = 180
    mascon_data = np.array(csr['lwe_thickness'])[index, int(res * 2)::int(res * 4),
                  int(res * 2)::int(res * 4)]  # unit: cm

    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    region = [min(lon), max(lon), min(lat), max(lat)]

    mascon_data = np.roll(mascon_data, shift=int(-180 / res), axis=-1)*10 # convert to mm

    bsd.create_dataset(name='lat', data=lat)
    bsd.create_dataset(name='lon', data=lon)
    bsd.create_dataset(name='Load_mascon_ewh', data=mascon_data)

    dt = h5py.special_dtype(vlen=str)
    bsd.create_dataset(name='description', dtype=dt, data=['The unit is mm.'])

    '''spherical harmonic res'''
    SH1 = np.load('../res/fig2/SH.npy')[0] * 1000
    bsd.create_dataset(name='Uplift_spherical_harmonic', data=SH1)

    '''Finite LGF'''
    Green = np.load('../res/fig2/Green90.npy') * 1000
    Green = Green.reshape((90,180))
    bsd.create_dataset(name='Uplift_finite_LGF', data=Green)

    '''Infinite LGF'''
    shen = np.loadtxt('../temp/loadvcd.txt', skiprows=1, usecols=2)
    shen = shen.reshape((90, 180))
    bsd.create_dataset(name='Uplift_infinite_LGF', data=shen)

    '''D-LGF'''
    d_LGF = np.load('../temp/100_dis_fast.npy') * 1000
    d_LGF = d_LGF.reshape((90, 180))
    bsd.create_dataset(name='Uplift_DLGF', data=d_LGF)

    '''CD-LGF'''
    cd_LGF = np.load('../res/fig15/cluster_disk_0.1.npy') * 1000
    cd_LGF = cd_LGF.reshape((90, 180))
    bsd.create_dataset(name='Uplift_CDLGF', data=cd_LGF)

    pass


if __name__ == '__main__':
    # fig10()
    # fig1()
    fig2_prepare()
    # fig2()
    # fig8_prepare()
    # save_data_to_repository()
    # fig7()
