import numpy as np
from src_EarthModel.LoveNumber import LoveNumber, LLN_Data, LLN_variable, Frame
from src_Basis.GreenFunction.DiskLoad import DiskLoad, DiskLoad_constrain
from src_auxiliary.EnumClasses import Displacement
from src_auxiliary.GeoMathKit import GeoMathKit
from src_LoadMesh.RegularGrid import regular


class GFA_displacement:
    """
    Green Function Approach (GFA): based on disk load green function;
    """

    def __init__(self, lln: LoveNumber, verbose=False):
        self._grids = None
        self._DL = None
        self._lln = lln
        self._verbose = verbose

    def configure(self, grids: dict, cf=None):
        """
        define the unit uniform disk loads (thickness = 1 meter), which is dependent on the radius.
        Avoid repeating definition of the disk by telling the radius.
        :param grids: ('lat', 'lon', 'radius', 'EWH') of the grid; radius [degree], lat, lon [degree], EWH [meter].
        EWH could be actually a matrix containing the dimension of time to speed up computation.
        :param cf: None (default): full D-LGF, no truncation (as in Yang et al., 2025).
                   A number K: the response beyond K * radius is set to zero (approximation, use with care;
                   the previous default K=15 silently removed the far-field deformation).
        :return:
        """

        Unit_Disk_loads = {}
        for rr in grids['radius']:
            if rr not in Unit_Disk_loads.keys():
                if cf is None:
                    Unit_Disk_loads[rr] = DiskLoad(lln=self._lln, radius=rr, thickness=1)
                else:
                    Unit_Disk_loads[rr] = DiskLoad_constrain(lln=self._lln, radius=rr, thickness=1,
                                                             constrain_factor=cf)

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
        '''output at the evaluation points (any number, e.g. GNSS stations), one column per epoch'''
        dis = np.zeros((np.size(points['lat']), np.shape(grids['EWH'])[1]))
        for id, rr in enumerate(list(grids['radius'])):
            if self._verbose:
                print(id)
            '''calculate the displacement: Vertical, Geoheight, Horizontal (away from load), North or East'''
            a = GeoMathKit.kernel_at_points(dl[rr], grids['lat'][id], grids['lon'][id], points['lat'], points['lon'],
                                             variable)
            dis += a[:, None] @ grids['EWH'][id][None, :]

            '''release memory'''
            dl[rr].release()

        return dis



class GFA_regular_grid(GFA_displacement):

    def __init__(self, lln: LoveNumber, verbose=False):
        super().__init__(lln, verbose=verbose)

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
            if self._verbose:
                print(id)
            if grids['lat'][id] != lat0:
                lat0 = grids['lat'][id]
                '''calculate the displacement: Vertical, Geoheight, Horizontal, North or East (roll-invariant)'''
                temp = GeoMathKit.kernel_at_points(dl[rr], grids['lat'][id], grids['lon'][id], points['lat'], points['lon'],
                                                    variable)

            else:
                temp = np.roll(temp.reshape((Num, -1)), shift=1, axis=1)
                temp = temp.flatten()

            dis += temp[:, None] @ grids['EWH'][id][None, :]

            '''release memory'''
            dl[rr].release()

        return dis


def grid2radius(lat_center, grid_size):
    """Equal-area disk radius of a cell; see src_LoadMesh.RegularGrid.regular.grid2radius (kept here for compatibility)."""
    return regular.grid2radius(lat_center=lat_center, grid_size=grid_size)


def grid2radius_type2(lat_center, grid_size):
    """Spherical-cap version; see src_LoadMesh.RegularGrid.regular.grid2radius_type2 (kept here for compatibility)."""
    return regular.grid2radius_type2(lat_center=lat_center, grid_size=grid_size)


def Mass(r, thickness):
    """
    calculate the radius of disk load that equals to the area of pixel
    :param lat_center: center of the pixel, latitude, [degree]
    :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
    :return: theta_radius [degree], length_radius [m]
    """
    from src_auxiliary.constants import EarthConstant

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
