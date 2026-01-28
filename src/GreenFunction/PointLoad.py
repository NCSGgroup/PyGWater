import numpy as np

from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.Legendre import Legendre_polynomial
from src.Auxiliary.constants import EarthConstant

"""
REF:Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0
"""


class LGF_truncation:
    """
    This is theoretically right but practical wrong because of limited truncation. Therefore, a compensation should be applied, see LGF
    """

    def __init__(self, lln: LoveNumber):
        self._residence = None
        self._lln = lln

    def configure(self, residence):
        '''

        :param residence: deg, angular distance between two nodes
        :return:
        '''

        residence[residence < 1e-5] = 1e-5  # to avoid singularity
        self._residence = residence
        self._legendre()
        return self

    def _legendre(self):
        self._Pl = Legendre_polynomial(alpha=self._residence, lmax=self._lln.lmax)
        pass

    def getVertical(self):
        item = self._Pl * self._lln.LLN[LLN_variable.h]

        return np.sum(item, axis=-1) * EarthConstant.radiusm / EarthConstant.Mass

    def getHorizental(self):
        n = np.arange(self._lln.lmax + 1)

        sinTheta = np.sin(np.deg2rad(self._residence))
        cosTheta = np.cos(np.deg2rad(self._residence))

        Plb = np.roll(self._Pl, 1, axis=-1)
        differential = n / sinTheta * (cosTheta * self._Pl - Plb)

        return differential[:, 1:] @ self._lln.LLN[LLN_variable.l][1:,
                                     None] * EarthConstant.radiusm / EarthConstant.Mass

    def getGeoidheight(self):
        return self._Pl @ (self._lln.LLN[LLN_variable.k][:, None] + 1) * EarthConstant.radiusm / EarthConstant.Mass


class LGF_Kummer(LGF_truncation):
    """
    LGF using Kummer's transformation (1st order approximation) to compensate the loss caused by truncation.
    """

    def __init__(self, lln: LoveNumber):
        super().__init__(lln)
        pass

    def getHorizental(self):
        """
        Exclude degree-0
        :return:
        """
        N = self._lln.lmax
        l_inf = self._lln.LLN[LLN_variable.l][-1] * N
        n = np.arange(N + 1)

        sinTheta = np.sin(np.deg2rad(self._residence))

        sinTheta[sinTheta < 1e-5] = 1e-5  # to avoid singularity
        cosTheta = np.cos(np.deg2rad(self._residence))
        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))
        cosTheta2 = np.cos(np.deg2rad(self._residence / 2))
        item2 = cosTheta2 * (1 + 2 * sinTheta2) / (2 * sinTheta2 * (1 + sinTheta2))

        if len(np.shape(sinTheta)) != 0:
            sinTheta = sinTheta[:, None]
            cosTheta = cosTheta[:, None]

        Plb = np.roll(self._Pl, 1, axis=-1)
        differential = 1 / sinTheta * (cosTheta * self._Pl - Plb) * n

        n[0] = 1  # to avoid singularity. degree-0 does not participate in the computation.
        item = differential * (n * self._lln.LLN[LLN_variable.l] - l_inf) / n

        return (np.sum(item[:, 1:], axis=-1) -
                l_inf * item2) * EarthConstant.radiusm / EarthConstant.Mass

    def getVertical(self):
        """
        Include degree-0
        :return:
        """
        h_inf = self._lln.LLN[LLN_variable.h][-1]

        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))

        item = self._Pl * (self._lln.LLN[LLN_variable.h] - h_inf)

        return (np.sum(item, axis=-1) + h_inf * (1 / (2 * sinTheta2))) * EarthConstant.radiusm / EarthConstant.Mass

    def getVertical_exclude_d0(self):
        """
        Exclude degree-0
        :return:
        """
        h_inf = self._lln.LLN[LLN_variable.h][-1]

        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))

        item = self._Pl * (self._lln.LLN[LLN_variable.h] - h_inf)

        return (np.sum(item, axis=-1) + h_inf * (1 / (2 * sinTheta2) - 1)) * EarthConstant.radiusm / EarthConstant.Mass

    def getGeoidheight(self):
        """
        Exclude degree-0.
        :return:
        """

        N = self._lln.lmax
        k_inf = self._lln.LLN[LLN_variable.k][-1] * N
        n = np.arange(N + 1)

        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))
        cosTheta2 = np.cos(np.deg2rad(self._residence / 2))
        item2 = 1 / (2 * sinTheta2) - 1 - k_inf * np.log(sinTheta2 + sinTheta2 ** 2)

        # if len(np.shape(sinTheta)) != 0:
        #     sinTheta = sinTheta[:, None]
        #     cosTheta = cosTheta[:, None]
        n[0] = 1  # to avoid singularity. degree-0 does not participate in the computation.
        item = self._Pl * (self._lln.LLN[LLN_variable.k] - k_inf / n)

        return (np.sum(item[:, 1:], axis=-1) + item2) * EarthConstant.radiusm / EarthConstant.Mass


def demo():
    from src.Auxiliary.FileTool import FileTool
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.PREM).get_Love_number()

    file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
    theta = np.loadtxt(file, skiprows=1, usecols=0)

    # load = PointLoad_theory(lln=lln).configure(residence=theta)
    load = LGF_Kummer(lln=lln).configure(residence=theta)
    uss = load.getVertical()
    u = load.getVertical() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    v = load.getHorizental() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    g = load.getGeoidheight() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)

    # u0 = load.getVertical_exclude_d0() * (10 ** 12) * EarthConstant.radiusm * np.deg2rad(theta)
    pass


if __name__ == '__main__':
    demo()
