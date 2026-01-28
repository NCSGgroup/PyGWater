import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.PointLoad import LGF_truncation
from src.Auxiliary.EnumClasses import Displacement
from src.Auxiliary.GeoMathKit import GeoMathKit
from src.Auxiliary.constants import EarthConstant


class GFA_displacement:
    """
    Green Function Approach (GFA): based on truncated point mass load green function;
    !! The maximal degree of LLN should be consistent with the resolution of load!!!
    For example, for a load at 2 degree, the lmax of LLN should be 180/2 = 90; otherwise, it would introduce large errors.
    """

    def __init__(self, lln: LoveNumber):
        self._grids = None
        # self._PL = LGF(lln=lln)
        self._PL = LGF_truncation(lln=lln)
        self._lln = lln

    def configure(self, grids: dict):
        """
        define the unit uniform disk loads (thickness = 1 meter), which is dependent on the radius.
        Avoid repeating definition of the disk by telling the radius.
        :param grids: ('lat', 'lon', 'area', 'EWH') of the grid; area [m**2], lat, lon [degree], EWH [meter].
        EWH could be actually a matrix containing the dimension of time to speed up computation.
        :return:
        """
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
        pl = self._PL
        dis = np.zeros_like(grids['EWH'])
        for id, rr in enumerate(list(grids['area'])):
            print(id)
            theta = GeoMathKit.angular_distance(grids['lat'][id], grids['lon'][id], points['lat'],
                                                points['lon'])
            '''To avoid singularity'''
            # theta[theta == 0] += 1e-6
            # theta[theta == 180] += 1e-6

            '''calculate the displacement'''
            a= self._getFunc(PL=pl.configure(residence=theta), variable=variable)
            temp = (a * rr)[:, None]
            dis += temp @ grids['EWH'][id][None, :]

        return dis

    def _getFunc(self, variable: Displacement, PL: LGF_truncation):
        if variable == Displacement.Vertical:
            return PL.getVertical()
        elif variable == Displacement.Horizontal:
            return PL.getHorizental()
        elif variable == Displacement.Geoheight:
            return PL.getGeoidheight()


class GFA_regular_grid(GFA_displacement):
    """
    This fast version only works when the grid network is defined as equal angular distance grid.
    """

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
        pl = self._PL
        dis = np.zeros_like(grids['EWH'])

        Num = int(180 / resolution)
        lat0 = -1000
        for id, rr in enumerate(list(grids['area'])):
            print(id)
            if grids['lat'][id] != lat0:
                lat0 = grids['lat'][id]
                theta = GeoMathKit.angular_distance(grids['lat'][id], grids['lon'][id], points['lat'],
                                                    points['lon'])
                '''To avoid singularity'''
                # theta[theta == 0] += 1e-6
                # theta[theta == 180] += 1e-6

                '''calculate the displacement'''
                a = self._getFunc(PL=pl.configure(residence=theta), variable=variable)
                temp = a * rr
            else:
                temp = np.roll(temp.reshape((Num, -1)), shift=1, axis=1)
                temp = temp.flatten()

            dis += temp[:, None] @ grids['EWH'][id][None, :]

        return dis


class Grids_generation:

    @staticmethod
    def Equal_angular_distance(resolution=0.5, Earth_radius=EarthConstant.radiusm):
        """

        :param Earth_radius: [meter]
        :param resolution: [degree]
        :return: grid, dict
        """

        lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=resolution)
        lon, lat = np.meshgrid(lon, lat)

        area = np.cos(np.deg2rad(lat)) * np.deg2rad(resolution) ** 2 * Earth_radius ** 2

        grids = {
            'lat': lat.flatten(),
            'lon': lon.flatten(),
            'area': area.flatten()
        }

        return grids
