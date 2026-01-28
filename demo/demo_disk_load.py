import pygmt
import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.DiskLoad import DiskLoad, EarthConstant


def demo1():
    from src.Auxiliary.FileTool import FileTool
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.PREM).get_Love_number()

    # file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    # theta = np.loadtxt(file, skiprows=1, usecols=0)
    # theta = 0.001100220044009
    radius = 0.10
    thickness = 1
    theta = np.linspace(0, radius * 25, 5000)
    load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    u = load.getVertical()
    v = load.getHorizental()
    g = load.getGeoidheight()

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='8p', FONT_TITLE='9p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne+tDisk radius: %s degree; Load thickness: %s meter" % (radius, thickness),
                       "xa0.5f0.5+l@[\\theta/\\alpha@[", "ya0.5f0.5+l mm"], region=[0, 25, -2.5, 1],
                # Set a logarithmic transformation on the x-axis
                projection="X8c/8c")
    ypoints = {'U': u * 1000, 'V': v * 1000, 'G': g * 1000}
    xpoints = theta / radius

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="thick,+z,-",
        )
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend()

    fig.show()
    pass


def demo2():
    """
    To calculate the truncation error
    :return:
    """

    from src.Auxiliary.FileTool import FileTool
    lmax_list = np.arange(100, 40100, 100)
    radius = 0.1
    thickness = 1
    K= 0.1
    theta = K * radius

    data_key = {'U': [], 'V': [], 'G': []}

    for lmax in lmax_list:
        print(lmax)
        lln = LoveNumber().config(lmax=lmax, method=LLN_Data.REF).get_Love_number()
        load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
        u = load.getVertical()
        v = load.getHorizental()
        g = load.getGeoidheight()
        data_key['U'].append(u)
        data_key['V'].append(v)
        data_key['G'].append(g)

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='15p', FONT_ANNOT_PRIMARY='8p', FONT_TITLE='7p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne+tTruncation error: Loading response computed at %s*radius as a function of l@-max@-"%K,
                       "xa1pf3+ll@-max@-", "ya1f0.5+l mm"], region=[100, 40000, -3, 1.5],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    ypoints= {}
    for key, value in data_key.items():
        ypoints[key] = np.array(value)*1000

    xpoints = lmax_list

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="2p,+z,-",
        )
        fig.plot(x=[100, 40000], y=[value[-1], value[-1]], cmap=True,
            zvalue=zvalue, pen='1p,--,+z')
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend(position='jCR',box=True)

    nROT = 360/radius
    fig.plot(x=[nROT, nROT], y=[-200, 200], pen='1p,Black,--')

    fig.plot(x=[2*nROT, 2*nROT], y=[-200, 200],  pen='1p,Black,--')

    fig.show()

    pass


def demo3():
    from src.Auxiliary.FileTool import FileTool

    LLN_list = [LLN_Data.PREM,LLN_Data.PREMhard,LLN_Data.PREMsoft]
    U,V,G={},{},{}
    for llndata in LLN_list:
        lln = LoveNumber().config(lmax=40000, method=llndata).get_Love_number()

        # file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
        # theta = np.loadtxt(file, skiprows=1, usecols=0)
        # theta = 0.001100220044009
        radius = 0.1
        thickness = 1
        theta = np.linspace(0, radius * 5, 5000)
        load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
        u = load.getVertical()
        v = load.getHorizental()
        g = load.getGeoidheight()
        U[llndata] = u*1000
        V[llndata] = v*1000
        G[llndata] = g*1000


    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='8p', FONT_TITLE='9p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne+tDisk radius: %s degree; Load thickness: %s meter" % (radius, thickness),
                       "xa0.5f0.5+l@[\\theta/\\alpha@[", "ya0.1f0.1+l mm"], region=[0, 2.5, -0.4, 0.2],
                # Set a logarithmic transformation on the x-axis
                projection="X8c/8c")
    ypoints = {'U': U[LLN_Data.PREM]-U[LLN_Data.PREMhard],
               'V': V[LLN_Data.PREM]-V[LLN_Data.PREMhard],
               'G': G[LLN_Data.PREM]-G[LLN_Data.PREMhard]}
    xpoints = theta / radius

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="thick,+z,-",
        )
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend()

    fig.text(
        text="(a) PREM - PREM_hard",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    ''''''
    fig.shift_origin(xshift='9c')

    fig.basemap(frame=["WSne+tDisk radius: %s degree; Load thickness: %s meter" % (radius, thickness),
                       "xa0.5f0.5+l@[\\theta/\\alpha@[", "ya0.1f0.1"], region=[0, 2.5, -0.3, 0.4],
                # Set a logarithmic transformation on the x-axis
                projection="X8c/8c")
    ypoints = {'U': U[LLN_Data.PREM] - U[LLN_Data.PREMsoft],
               'V': V[LLN_Data.PREM] - V[LLN_Data.PREMsoft],
               'G': G[LLN_Data.PREM] - G[LLN_Data.PREMsoft]}
    xpoints = theta / radius

    zvalue = 4
    for key, value in ypoints.items():
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=value,
            label=key,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="thick,+z,-",
        )
        zvalue += 1

    with pygmt.config(FONT_ANNOT_PRIMARY='10p'):
        fig.legend()

    fig.text(
        text="(b) PREM - PREM_soft",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    fig.show()
    pass


def demo4():
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.PREM).get_Love_number()

    # file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    # theta = np.loadtxt(file, skiprows=1, usecols=0)
    # theta = 0.001100220044009
    radius = 0.10
    thickness = 1
    theta = np.linspace(0, radius * 5, 5000)
    load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    u1 = load.getVertical()
    v1 = load.getHorizental()
    g1 = load.getGeoidheight()

    radius = 0.10
    thickness = 2
    theta = np.linspace(0, radius * 5, 5000)
    load = DiskLoad(lln=lln, radius=radius, thickness=thickness).configure(residence=theta)
    u2 = load.getVertical()
    v2 = load.getHorizental()
    g2 = load.getGeoidheight()

    pass

if __name__ == '__main__':
    demo3()
