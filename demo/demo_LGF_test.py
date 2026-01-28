import pygmt
import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.PointLoad import LGF_Kummer, EarthConstant, LGF_truncation
from src.GreenFunction.DiskLoad import DiskLoad, EarthConstant
from src.Auxiliary.FileTool import FileTool


def compare():
    lln = LoveNumber().config(lmax=90, method=LLN_Data.PREM).get_Love_number()

    file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    theta = np.loadtxt(file, skiprows=1, usecols=0)
    theta[-1] -= 0.01

    '''point load'''
    M = 1 #kg
    # load = LGF(lln=lln).configure(residence=theta)
    load = LGF_Kummer(lln=lln).configure(residence=theta)
    # u1 = load.getVertical() * (10 ** 12) * EarthConstant.radiusm/40
    # v1 = load.getHorizental() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    # g1 = load.getGeoidheight() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)

    # radius = 0.5
    # area = 2*np.pi*(1-np.cos(np.deg2rad(radius)))*EarthConstant.radiusm**2
    # thickness = M/EarthConstant.rhow/area
    # load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    # u1 = load.getVertical()* (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)

    load = LGF_truncation(lln=lln).configure(residence=theta)
    u1 = load.getVertical() * (10 ** 12) * EarthConstant.radiusm/40

    '''disk load'''
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    radius = 2
    area = 2*np.pi*(1-np.cos(np.deg2rad(radius)))*EarthConstant.radiusm**2
    thickness = M/EarthConstant.rhow/area
    load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    u2 = load.getVertical()* (10 ** 12) * EarthConstant.radiusm/40
    # v2 = load.getHorizental()* (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    # g2 = load.getGeoidheight()* (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)




    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='8p', FONT_TITLE='9p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne+tDisk radius: %s degree; Load thickness: %s meter" % (radius, thickness),
                       "xa10f5", "ya10f5"], region=[-1, 10, -50, 5],
                # Set a logarithmic transformation on the x-axis
                projection="X8c/8c")
    xpoints = theta

    zvalue = 1

    fig.plot(
        # Set the figures frame and color as well as
        # annotations, ticks, and gridlines
        x=xpoints,
        y=u1,
        label='point load',
        # Set the line thickness to "1p", the color to "blue",
        # and the style to "-", i.e. "dashed"
        cmap=True,
        zvalue=zvalue,
        pen="thick,+z,-",
    )

    zvalue += 1

    fig.plot(
        # Set the figures frame and color as well as
        # annotations, ticks, and gridlines
        x=xpoints,
        y=u2,
        label='disk load',
        # Set the line thickness to "1p", the color to "blue",
        # and the style to "-", i.e. "dashed"
        cmap=True,
        zvalue=zvalue,
        pen="thick,+z,-",
    )
    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend()


    fig.shift_origin(xshift='10c')

    fig.basemap(frame=["WSne+tDisk radius: %s degree; Load thickness: %s meter" % (radius, thickness),
                       "xa0.5f0.5", "ya1f1"], region=[0, 4, 0, 4],
                # Set a logarithmic transformation on the x-axis
                projection="X8c/8c")
    xpoints = theta / radius

    zvalue = 1

    u3 = u2/u1
    u3[xpoints<1] = u1[xpoints<1]/u2[xpoints<1]

    fig.plot(
        # Set the figures frame and color as well as
        # annotations, ticks, and gridlines
        x=xpoints,
        y=u3,
        label='point load',
        # Set the line thickness to "1p", the color to "blue",
        # and the style to "-", i.e. "dashed"
        cmap=True,
        zvalue=zvalue,
        pen="thick,+z,-",
    )


    fig.show()

    pass


if __name__ == '__main__':
    compare()