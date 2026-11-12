import numpy as np
from pygwater.earth.love_numbers import LoveNumber, LLN_Data, LLN_variable, Frame
from pygwater.basis.green.point_load import LGF_truncation
from pygwater.utils.enums import Displacement
from pygwater.utils.geomath import GeoMathKit
from pygwater.utils.constants import EarthConstant


class GFA_displacement:
    """
    Green Function Approach (GFA): based on truncated point mass load green function;
    !! The maximal degree of LLN should be consistent with the resolution of load!!!
    For example, for a load at 2 degree, the lmax of LLN should be 180/2 = 90; otherwise, it would introduce large errors.
    """

    def __init__(self, lln: LoveNumber, verbose=False):
        self._grids = None
        # self._PL = LGF(lln=lln)
        self._PL = LGF_truncation(lln=lln)
        self._lln = lln
        self._verbose = verbose

    def configure(self, grids: dict):
        """
        :param grids: ('lat', 'lon', 'area', 'EWH') of the grid; area [m**2], lat, lon [degree].
        NOTE: despite its name, grids['EWH'] must be the surface MASS DENSITY [kg/m**2], i.e. EWH [m] * rho_w
        (see Eq. 8 of Yang et al., 2025), shape (Ncell, Ntime) to process a time series in one pass.
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
        '''output at the evaluation points (any number, e.g. GNSS stations), one column per epoch'''
        dis = np.zeros((np.size(points['lat']), np.shape(grids['EWH'])[1]))
        for id, rr in enumerate(list(grids['area'])):
            if self._verbose:
                print(id)
            '''calculate the displacement: Vertical, Geoheight, Horizontal (away from load), North or East'''
            a = GeoMathKit.kernel_at_points(pl, grids['lat'][id], grids['lon'][id], points['lat'], points['lon'],
                                             variable)
            temp = (a * rr)[:, None]
            dis += temp @ grids['EWH'][id][None, :]

        return dis



class GFA_regular_grid(GFA_displacement):
    """
    This fast version only works when the grid network is defined as equal angular distance grid.
    """

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
        pl = self._PL
        dis = np.zeros_like(grids['EWH'])

        Num = int(180 / resolution)
        lat0 = -1000
        for id, rr in enumerate(list(grids['area'])):
            if self._verbose:
                print(id)
            if grids['lat'][id] != lat0:
                lat0 = grids['lat'][id]
                '''calculate the displacement: Vertical, Geoheight, Horizontal, North or East. North/East
                components are invariant under a rotation about the polar axis, so rolling in longitude is valid.'''
                a = GeoMathKit.kernel_at_points(pl, grids['lat'][id], grids['lon'][id], points['lat'], points['lon'],
                                                 variable)
                temp = a * rr
            else:
                temp = np.roll(temp.reshape((Num, -1)), shift=1, axis=1)
                temp = temp.flatten()

            dis += temp[:, None] @ grids['EWH'][id][None, :]

        return dis


class Grids_generation:
    """Kept for compatibility; see pygwater.mesh.regular_grid.regular.global_grid_for_pointmass."""

    @staticmethod
    def Equal_angular_distance(resolution=0.5, Earth_radius=EarthConstant.radiusm):
        """
        :param resolution: [degree]
        :param Earth_radius: [meter]
        :return: dict with flattened 'lat', 'lon' [degree] and 'area' [m^2] of a global cell-centred grid
        """
        from pygwater.mesh.regular_grid import regular
        return regular.global_grid_for_pointmass(resolution=resolution, Earth_radius=Earth_radius)
