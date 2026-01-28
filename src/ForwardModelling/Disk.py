import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.DiskLoad import DiskLoad_constrain
from src.Auxiliary.EnumClasses import Displacement
from src.Auxiliary.GeoMathKit import GeoMathKit


class GFA_displacement:
    """
    Green Function Approach (GFA): based on disk load green function;
    """

    def __init__(self, lln: LoveNumber):
        self._grids = None
        self._DL = None
        self._lln = lln

    def configure(self, grids: dict, cf=15):
        """
        define the unit uniform disk loads (thickness = 1 meter), which is dependent on the radius.
        Avoid repeating definition of the disk by telling the radius.
        :param grids: ('lat', 'lon', 'radius', 'EWH') of the grid; radius [degree], lat, lon [degree], EWH [meter].
        EWH could be actually a matrix containing the dimension of time to speed up computation.
        :return:
        """

        Unit_Disk_loads = {}
        for rr in grids['radius']:
            if rr not in Unit_Disk_loads.keys():
                Unit_Disk_loads[rr] = DiskLoad_constrain(lln=self._lln, radius=rr, thickness=1, constrain_factor=cf)

        self._DL = Unit_Disk_loads
        self._grids = grids
        return self

    def evaluation(self, points: dict, variable=Displacement.Vertical):
        """
        Evaluation of desired variable at specified points
        :param points: ('lat', 'lon')
        :param variable:
        :return:
        """
        grids = self._grids
        dl = self._DL
        dis = np.zeros_like(grids['EWH'])
        for id, rr in enumerate(list(grids['radius'])):
            print(id)
            theta = GeoMathKit.angular_distance(grids['lat'][id], grids['lon'][id], points['lat'],
                                                points['lon'])

            '''To avoid singularity'''
            theta[theta == 0] += 1e-6
            theta[theta == 180] += 1e-6

            '''calculate the displacement'''
            dis += self._getFunc(DL=dl[rr].configure(residence=theta), variable=variable)[:, None] @ grids['EWH'][id][
                                                                                                     None, :]

            '''release memory'''
            dl[rr].release()

        return dis

    def _getFunc(self, variable: Displacement, DL: DiskLoad_constrain):
        if variable == Displacement.Vertical:
            return DL.getVertical()
        elif variable == Displacement.Horizontal:
            return DL.getHorizental()
        elif variable == Displacement.Geoheight:
            return DL.getGeoidheight()


class GFA_regular_grid(GFA_displacement):

    def __init__(self, lln: LoveNumber):
        super().__init__(lln)

    def evaluation(self, points: dict, variable=Displacement.Vertical, resolution=2):
        """
        Evaluation of desired variable at specified points
        :param resolution:
        :param points: ('lat', 'lon')
        :param variable:
        :return:
        """
        grids = self._grids
        dl = self._DL
        dis = np.zeros_like(grids['EWH'])

        Num = int(180 / resolution)
        lat0 = -1000
        for id, rr in enumerate(list(grids['radius'])):
            print(id)
            if grids['lat'][id] != lat0:
                lat0 = grids['lat'][id]
                theta = GeoMathKit.angular_distance(grids['lat'][id], grids['lon'][id], points['lat'],
                                                    points['lon'])
                '''To avoid singularity'''
                theta[theta == 0] += 1e-6
                theta[theta == 180] += 1e-6

                '''calculate the displacement'''
                temp = self._getFunc(DL=dl[rr].configure(residence=theta), variable=variable)

            else:
                temp = np.roll(temp.reshape((Num, -1)), shift=1, axis=1)
                temp = temp.flatten()

            dis += temp[:, None] @ grids['EWH'][id][None, :]

            '''release memory'''
            dl[rr].release()

        return dis


def grid2radius(lat_center, grid_size):
    """
    calculate the radius of disk load that equals to the area of pixel
    :param lat_center: center of the pixel, latitude, [degree]
    :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
    :return: theta_radius [degree], length_radius [m]
    """
    from src.Auxiliary.constants import EarthConstant
    example = 30  # longitude, but indeed the choice could be arbitrary.

    a = GeoMathKit.angular_distance(point1_lat=lat_center + grid_size / 2, point1_lon=example,
                                    point2_lat=lat_center - grid_size / 2, point2_lon=example)

    b = GeoMathKit.angular_distance(point1_lat=lat_center, point1_lon=example, point2_lat=lat_center,
                                    point2_lon=example + grid_size)

    # print(np.deg2rad(a) * np.deg2rad(b))
    r = np.sqrt(np.deg2rad(a) * np.deg2rad(b) / np.pi)

    return np.rad2deg(r), EarthConstant.radiusm * r


def grid2radius_type2(lat_center, grid_size):
    """
    calculate the radius of disk load that equals to the area of pixel
    :param lat_center: center of the pixel, latitude, [degree]
    :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
    :return: theta_radius [degree], length_radius [m]
    """
    from src.Auxiliary.constants import EarthConstant

    area = np.cos(np.deg2rad(lat_center)) * np.deg2rad(grid_size) ** 2

    # print(area)
    x = 1 - area / (2 * np.pi)

    r = np.arccos(x)

    return np.rad2deg(r), EarthConstant.radiusm * r


def Mass(r, thickness):
    """
    calculate the radius of disk load that equals to the area of pixel
    :param lat_center: center of the pixel, latitude, [degree]
    :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
    :return: theta_radius [degree], length_radius [m]
    """
    from src.Auxiliary.constants import EarthConstant

    M1 = 2 * np.pi * EarthConstant.radiusm ** 2 * (
                1 - np.cos(r / EarthConstant.radiusm)) * EarthConstant.rhow * thickness

    M2 = np.pi * r ** 2 * EarthConstant.rhow * thickness / (1e12)

    return M1, M2


def demo1():
    a, b = grid2radius(lat_center=np.array([80]), grid_size=0.05)
    print(a)
    print(b)

    a, b = grid2radius_type2(lat_center=np.array([89]), grid_size=2)
    print(a)
    print(b)
    #
    print(Mass(r=b, thickness=0.068))
    pass


if __name__ == '__main__':
    demo1()
