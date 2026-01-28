from enum import Enum


class LLN_Data(Enum):
    PREM = 1
    REF = 2
    Wang = 3
    iasp91 = 4
    ak135 = 5
    iasp91hard = 6
    ak135hard = 7
    PREMhard = 8
    PREMsoft = 9


class LLN_variable(Enum):
    h = 1
    k = 2
    l = 3


class Frame(Enum):
    CM = 0
    CF = 1
    CE = 2


class Displacement(Enum):
    Vertical = 0
    Horizontal = 1
    Geoheight = 2
    Gravity = 3


class PhysicalDimensions(Enum):
    Dimensionless = 0
    EWH = 1
    Pressure = 2
    Density = 3
    Geoid = 4
    Gravity = 5
    HorizontalDisplacementEast = 6
    HorizontalDisplacementNorth = 7
    VerticalDisplacement = 8


def match_string(name, obj, ignore_case=False):
    obj_list = list(obj)
    names = [obj_list[i].name for i in range(len(obj_list))]

    if ignore_case:
        names = [names[i].lower() for i in range(len(names))]
        name = name.lower()

    assert name in names

    return obj_list[names.index(name)]