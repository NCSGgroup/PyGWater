import copy
import datetime
import pathlib
import warnings

import h5py
import netCDF4
import numpy as np

from src_auxiliary.MathTool import MathTool
from src_Basis.SphericalHarmonic.CoreGRID import CoreGRID
from src_Basis.SphericalHarmonic.CoreSHC import CoreSHC
import src_auxiliary.EnumClasses as Enums
from src_Basis.SphericalHarmonic.ConvertSHC import ConvertSHC
from src_Basis.SphericalHarmonic.Harmonic import Harmonic
from src_EarthModel.LoveNumber import LoveNumber, LLN_variable


class SHC(CoreSHC):
    def __init__(self, c, s=None, LLN: LoveNumber = None):
        """
        :param LLN: LoveNumber object used for unit conversion. If None, Wang et al. (2012) PREM Love numbers are
                    read lazily (only when a conversion is requested), with lmax >= max(1000, lmax of the SHC).
                    Previously the default argument LoveNumber() was evaluated once at import time.
        """
        super().__init__(c, s)
        self._LLN = None if LLN is None else LLN.LLN
        pass

    def _get_LLN(self, lmax):
        if self._LLN is None or len(self._LLN[LLN_variable.k]) < lmax + 1:
            self._LLN = LoveNumber().config(lmax=max(1000, lmax)).get_Love_number().LLN
        return self._LLN

    def __add__(self, other):
        assert issubclass(type(other), CoreSHC)

        out = SHC(self.value + other.value)
        out._LLN = self._LLN
        return out

    def __sub__(self, other):
        assert issubclass(type(other), CoreSHC)

        out = SHC(self.value - other.value)
        out._LLN = self._LLN
        return out

    def get_degree_rms(self):
        cqlm, sqlm = self.get_cs2d()
        return MathTool.get_degree_rms(cqlm, sqlm)

    def get_degree_rss(self):
        cqlm, sqlm = self.get_cs2d()
        return MathTool.get_degree_rss(cqlm, sqlm)

    def get_std(self):
        cs_std = np.std(self.value, axis=0)
        out = SHC(cs_std)
        out._LLN = self._LLN
        return out

    def convert_type(self, from_type=None, to_type=None):
        types = list(Enums.PhysicalDimensions)
        types_string = [i.name.lower() for i in types]
        types += types_string

        if from_type is None:
            from_type = Enums.PhysicalDimensions.Dimensionless
        if to_type is None:
            to_type = Enums.PhysicalDimensions.Dimensionless

        assert (from_type.lower() if type(
            from_type) is str else from_type) in types, f"from_type must be one of {types}"
        assert (to_type.lower() if type(
            to_type) is str else to_type) in types, f"to_type must be one of {types}"

        if from_type is None:
            from_type = Enums.PhysicalDimensions.Dimensionless
        if to_type is None:
            to_type = Enums.PhysicalDimensions.Dimensionless

        if type(from_type) is str:
            from_type = Enums.match_string(from_type, Enums.PhysicalDimensions, ignore_case=True)
        if type(to_type) is str:
            to_type = Enums.match_string(to_type, Enums.PhysicalDimensions, ignore_case=True)
        lmax = self.get_lmax()

        LLN = self._get_LLN(lmax)
        ln_k = LLN[LLN_variable.k][0:(lmax + 1)]
        ln_h = LLN[LLN_variable.h][0:(lmax + 1)]
        ln_l = LLN[LLN_variable.l][0:(lmax + 1)]

        convert = ConvertSHC()
        convert.configuration.set_Love_number(ln_k=ln_k, ln_h=ln_h, ln_l=ln_l).set_input_type(from_type).set_output_type(to_type)

        self.value = convert.apply_to(self.value)
        return self

    def to_grid(self, grid_space=None, special_type: Enums.PhysicalDimensions = None):
        """pure synthesis"""

        if grid_space is None:
            grid_space = int(180 / self.get_lmax())
        assert special_type in (
            None,
            Enums.PhysicalDimensions.HorizontalDisplacementEast,
            Enums.PhysicalDimensions.HorizontalDisplacementNorth,
        )

        lat, lon = MathTool.get_global_lat_lon_range(grid_space)

        lmax = self.get_lmax()
        har = Harmonic(lat, lon, lmax, option=1)

        cqlm, sqlm = self.get_cs2d()
        grid_data = har.synthesis(cqlm, sqlm, special_type=special_type)
        grid = GRID(grid_data, lat, lon, option=1)

        return grid

    def synthesis(self, grid_space, from_type: Enums.PhysicalDimensions = None,
                  to_type: Enums.PhysicalDimensions = None):

        shc_copy = copy.deepcopy(self)

        special_type = to_type if to_type in (
            Enums.PhysicalDimensions.HorizontalDisplacementNorth,
            Enums.PhysicalDimensions.HorizontalDisplacementEast) else None

        shc_copy.convert_type(from_type=from_type, to_type=to_type)
        grid = shc_copy.to_grid(grid_space=grid_space, special_type=special_type)

        return grid


    def synthesis_observable(self, grid_space, variable, lln: LoveNumber = None, **kwargs):
        """
        Synthesis of any loading observable (src_EarthModel/Observables.py) on a global grid, with THIS SHC being the
        load in EWH [m]. Synthesis only: no analysis of the observable is involved.
        :param grid_space: grid spacing [degree]
        :param variable: EnumClasses.Displacement: Vertical, North, East, Geoheight, Gravity, GravityDisturbance
        :param lln: LoveNumber object; default: the one attached to this SHC (or Wang PREM)
        :param kwargs: newtonian='below'|'above'|'none' (Gravity), height=[m] (GravityDisturbance)
        :return: GRID (m, or m/s^2 for gravity)
        """
        from src_EarthModel import Observables
        lmax = self.get_lmax()
        if lln is None:
            class _L:  # minimal LoveNumber-like wrapper around the attached Love numbers
                pass
            lln = _L()
            lln.LLN = self._get_LLN(lmax)
        lat, lon = MathTool.get_global_lat_lon_range(grid_space)
        cqlm, sqlm = self.get_cs2d()
        data = Observables.sh_synthesis(cqlm, sqlm, lln, variable, lat, lon, **kwargs)
        return GRID(data, lat, lon, option=1)


class GRID(CoreGRID):
    def __init__(self, grid, lat, lon, option=1):
        super().__init__(grid, lat, lon, option)

    def to_SHC(self, lmax=None, special_type: Enums.PhysicalDimensions = None, LLN: LoveNumber = None,
               method='auto'):
        """
        :param LLN: LoveNumber object attached to the returned SHC (used for later unit conversion), e.g. the same
                    Earth model as used by the Green's function bases. Default: Wang et al. (2012) PREM.
        :param method: 'auto' (default: 'dh' within the sampling limit lmax <= nlat/2-1, else 'midpoint'),
                       'dh' (Driscoll-Healy weights as in SaGEA), 'midpoint' (scheme used for the figures of
                       Yang et al., 2025) or 'least_squares', see Harmonic.analysis.
        """
        grid_space = self.get_grid_space()

        assert special_type in (
            None,
            Enums.PhysicalDimensions.HorizontalDisplacementEast,
            Enums.PhysicalDimensions.HorizontalDisplacementNorth,
        )

        if special_type in (
                Enums.PhysicalDimensions.HorizontalDisplacementEast,
                Enums.PhysicalDimensions.HorizontalDisplacementNorth):
            assert False, "Horizontal Displacement is not supported yet."

        if lmax is None:
            lmax = int(180 / grid_space)

        lat, lon = MathTool.get_global_lat_lon_range(grid_space)

        har = Harmonic(lat, lon, lmax, option=1)

        grid_data = self.value
        cqlm, sqlm = har.analysis(grid_data, special_type=special_type, method=method)
        shc = SHC(cqlm, sqlm, LLN=LLN)

        return shc

    def analysis(self, lmax=None, from_type: Enums.PhysicalDimensions = None,
                 to_type: Enums.PhysicalDimensions = None, LLN: LoveNumber = None, method='auto'):

        grid_copy = copy.deepcopy(self)

        shc = grid_copy.to_SHC(lmax=lmax, LLN=LLN, method=method)
        shc.convert_type(from_type=from_type, to_type=to_type)

        return shc


if __name__ == '__main__':
    pass
