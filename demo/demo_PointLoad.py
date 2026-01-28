import pygmt
import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.PointLoad import LGF_Kummer, EarthConstant


def demo1():
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
    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='8p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    fig.basemap(frame=["WSne", "xa1f3+l@[\\theta@[", "ya100f50+l@[u@["], region=[0.00001, 179, -400, 100],
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
            zvalue=zvalue,
            pen="thick,+z,-",
        )

    with pygmt.config(FONT_ANNOT_PRIMARY='9p'):
        fig.legend(position="jBR", box=True)
    fig.text(
        text="(a)",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    fig.shift_origin(xshift='9.7c')
    fig.basemap(frame=["WSne", "xa1f3+l@[\\theta@[", "ya40f20+l@[v@["], region=[0.00001, 179, -120, 40],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item][1]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            # label=item.name,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="thick,+z,-",
        )
    fig.text(
        text="(b)",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    fig.shift_origin(xshift='9.7c')
    fig.basemap(frame=["WSne", "xa1f3+l@[\\theta@[", "ya5f2.5+l@[N@["], region=[0.00001, 179, -15, 10],
                # Set a logarithmic transformation on the x-axis
                projection="X8cl/8c")

    for zvalue, item in enumerate(list_LLN):
        # if zvalue>1:break
        xpoints = theta
        ypoints = data[item][2]
        fig.plot(
            # Set the figures frame and color as well as
            # annotations, ticks, and gridlines
            x=xpoints,
            y=ypoints,
            # label=item.name,
            # Set the line thickness to "1p", the color to "blue",
            # and the style to "-", i.e. "dashed"
            cmap=True,
            zvalue=zvalue,
            pen="thick,2p,+z,-",
        )

    fig.text(
        text="(c)",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )

    fig.show()
    pass


if __name__ == '__main__':
    demo1()
