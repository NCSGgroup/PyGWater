import pygmt
import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame


def demo():
    data={}
    lmax= 10000
    for item in LLN_Data:
        data[item] = LoveNumber().config(lmax=lmax, method=item).get_Love_number().convert(target=Frame.CM)
        # data[item] = LoveNumber().config(lmax=lmax, method=item).get_Love_number()

    fig = pygmt.Figure()

    pygmt.config(FONT_LABEL='20p', FONT_ANNOT_PRIMARY='13p')
    pygmt.makecpt(cmap="categorical", series=[0, 9, 1], color_model="+c0-9")

    # fig.basemap(frame=["WSne+gbisque", "xa1f3", "ya2f1"],region=[1, 11000, -8, 0],
    #         # Set a logarithmic transformation on the x-axis
    #         projection="X10cl/10c")

    fig.basemap(frame=["WSne", "xa1f3+l@[n@[", "ya2f1+l@[h_n@["],region=[1, 11000, -8, 0],
            # Set a logarithmic transformation on the x-axis
            projection="X8cl/8c")

    for zvalue, item in enumerate(LLN_Data):
        # if zvalue>1:break
        xpoints = np.arange(1, lmax)
        ypoints = data[item].LLN[LLN_variable.h][xpoints]
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
        fig.legend()
    fig.text(
        text="(a)",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )
    # Plot square root values as points on the line
    # Style of points is 0.3 cm squares, color fill is "red" with a "black" outline
    # Points are not clipped if they go off the figure

    fig.shift_origin(xshift='9.7c')
    fig.basemap(frame=["WSne", "xa1f3+l@[n@[", "ya0.4f0.2+l@[nl_n@["],region=[1, 11000, 0, 2],
            # Set a logarithmic transformation on the x-axis
            projection="X8cl/8c")

    for zvalue, item in enumerate(LLN_Data):
        # if zvalue>1:break
        xpoints = np.arange(1, lmax)
        ypoints = data[item].LLN[LLN_variable.l][xpoints]*xpoints
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
    fig.basemap(frame=["WSne", "xa1f3+l@[n@[", "ya1f0.5+l@[nk_n@["],region=[1, 11000, -3, 0],
            # Set a logarithmic transformation on the x-axis
            projection="X8cl/8c")

    for zvalue, item in enumerate(LLN_Data):
        # if zvalue>1:break
        xpoints = np.arange(1, lmax)
        ypoints = data[item].LLN[LLN_variable.k][xpoints]*xpoints
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
        text="(c)",
        font='15p',
        position="TL",  # Top Left
        justify="TL",  # Top Left
        offset="0.1c/-0.1c",
    )


    fig.show()

    pass

if __name__ == '__main__':
    demo()