"""
KERNEL of a single point load: response K(psi) at angular distance psi [deg] (paper Eqs. 8, 12).
    LGF_truncation : finite LGF (truncated at the degree of the Love numbers)
    LGF_Kummer     : point-load LGF with Kummer's transformation
The basis built from these kernels (many cells, any stations, design matrix) is gf_basis.py.
"""
import numpy as np

from pygwater.earth.love_numbers import LoveNumber, LLN_Data, LLN_variable, Frame
from pygwater.basis.green.legendre import Legendre_polynomial
from pygwater.utils.constants import EarthConstant

"""
REF:Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0
"""


class LGF_truncation:
    """
    This is theoretically right but practical wrong because of limited truncation. Therefore, a compensation should be applied, see LGF
    """

    def __init__(self, lln: LoveNumber):
        self._residence = None
        self._scalar = False
        self._lln = lln

    def configure(self, residence):
        '''

        :param residence: deg, angular distance between two nodes (scalar or array)
        :return:
        '''
        self._scalar = np.ndim(residence) == 0
        residence = np.atleast_1d(np.array(residence, dtype=float))  # copy: the caller's array is not modified
        residence[residence < 1e-5] = 1e-5  # to avoid singularity
        residence[residence > 180 - 1e-5] = 180 - 1e-5  # antipode: 1/sin(theta) in the horizontal kernel
        self._residence = residence
        self._legendre()
        return self

    def _out(self, x):
        """all getters return shape (N,) for array input and a scalar for scalar input"""
        x = np.asarray(x).reshape(-1)
        return x[0] if self._scalar else x

    def _legendre(self):
        self._Pl = Legendre_polynomial(alpha=self._residence, lmax=self._lln.lmax)
        pass

    def getVertical(self):
        """Eq. 7 of Yang et al. (2025)"""
        item = self._Pl * self._lln.LLN[LLN_variable.h]

        return self._out(np.sum(item, axis=-1) * EarthConstant.radiusm / EarthConstant.Mass)

    def getHorizental(self):
        n = np.arange(self._lln.lmax + 1)

        sinTheta = np.sin(np.deg2rad(self._residence))[:, None]  # broadcast over degrees (was missing)
        cosTheta = np.cos(np.deg2rad(self._residence))[:, None]

        Plb = np.roll(self._Pl, 1, axis=-1)
        differential = n / sinTheta * (cosTheta * self._Pl - Plb)

        return self._out(differential[:, 1:] @ self._lln.LLN[LLN_variable.l][1:]
                         * EarthConstant.radiusm / EarthConstant.Mass)

    def getGeoidheight(self):
        return self._out(self._Pl @ (self._lln.LLN[LLN_variable.k] + 1) * EarthConstant.radiusm / EarthConstant.Mass)

    def getGravity(self, newtonian='below', include_degree0=True):
        """
        Surface gravity change per kg of load [m/s^2/kg], positive = increase (Farrell, 1972):
            g/M * sum_n   (g/M = G/a^2, consistent with M; the Newtonian part is then exactly G/chord^2) [newton_n + (n+1) k_n - 2 h_n] P_n
        newton_n = -n ('below'), n+1 ('above') or 0 ('none'), see src/earth/observables.py. Truncated series: as
        for the other finite LGFs, it is equivalent to the spherical-harmonic solution only if lmax matches the
        resolution of the load (Sect. 4 of Yang et al., 2025).
        """
        n = np.arange(self._lln.lmax + 1)
        h, k = self._lln.LLN[LLN_variable.h], self._lln.LLN[LLN_variable.k]
        newton = {'below': -n, 'above': n + 1, 'none': 0 * n}[newtonian]
        c = (newton + (n + 1) * k - 2 * h).astype(float)
        if not include_degree0:
            c[0] = 0.
        return self._out(self._Pl @ c * EarthConstant.ggg / EarthConstant.radiusm ** 2)


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

        sinTheta = sinTheta[:, None]
        cosTheta = cosTheta[:, None]

        Plb = np.roll(self._Pl, 1, axis=-1)
        differential = 1 / sinTheta * (cosTheta * self._Pl - Plb) * n

        n[0] = 1  # to avoid singularity. degree-0 does not participate in the computation.
        item = differential * (n * self._lln.LLN[LLN_variable.l] - l_inf) / n

        return self._out((np.sum(item[:, 1:], axis=-1) -
                          l_inf * item2) * EarthConstant.radiusm / EarthConstant.Mass)

    def getVertical(self):
        """
        Eq. 14 of Yang et al. (2025), include degree-0.
        Note: h_inf is approximated by h at the truncation degree (lmax), i.e. lim n->inf h_n = h_inf (Eq. 12);
        likewise n*l_n -> l_inf and n*k_n -> k_inf for the horizontal and geoid kernels.
        :return:
        """
        h_inf = self._lln.LLN[LLN_variable.h][-1]

        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))

        item = self._Pl * (self._lln.LLN[LLN_variable.h] - h_inf)

        return self._out((np.sum(item, axis=-1) + h_inf * (1 / (2 * sinTheta2)))
                         * EarthConstant.radiusm / EarthConstant.Mass)

    def getVertical_exclude_d0(self):
        """
        Exclude degree-0
        :return:
        """
        h_inf = self._lln.LLN[LLN_variable.h][-1]

        sinTheta2 = np.sin(np.deg2rad(self._residence / 2))

        item = self._Pl * (self._lln.LLN[LLN_variable.h] - h_inf)

        return self._out((np.sum(item, axis=-1) + h_inf * (1 / (2 * sinTheta2) - 1))
                         * EarthConstant.radiusm / EarthConstant.Mass)

    def getGravity(self, newtonian='below', include_degree0=True):
        """
        Surface gravity change per kg of load [m/s^2/kg], positive = increase, with Kummer's transformation
        (Farrell, 1972; Wang et al., 2012). Degree 0 included by default (as getVertical); the tabulated gravity
        Green's function of Wang et al. (2012, PREM-LGFs.dat) excludes it - with include_degree0=False this
        kernel reproduces that table to a constant 0.4% (different g, M, a constants) at all distances.
            elastic   : sum_n e_n P_n, e_n = (n+1) k_n - 2 h_n -> e_inf = k_inf - 2 h_inf,
                        = sum_n (e_n - e_inf) P_n + e_inf / (2 sin(theta/2))
            Newtonian : sum_n (-n) P_n = 1 / (4 sin(theta/2))   (theta > 0)
        Away from the load (theta > 0) the Newtonian kernel is the same for a sensor below or above the load
        (they differ only by the local Bouguer term at theta = 0), so newtonian='below' and 'above' are identical
        here; 'none' gives the elastic part only. Matches the gravity Green's function of Wang et al. (2012).
        """
        N = self._lln.lmax
        n = np.arange(N + 1)
        h, k = self._lln.LLN[LLN_variable.h], self._lln.LLN[LLN_variable.k]
        h_inf, k_inf = h[-1], k[-1] * N  # lim h_n = h_inf, lim n k_n = k_inf (Eq. 12, corrected)
        e = (n + 1) * k - 2 * h
        e_inf = k_inf - 2 * h_inf

        s2 = np.sin(np.deg2rad(self._residence / 2))
        elastic = self._Pl @ (e - e_inf) + e_inf / (2 * s2)
        if not include_degree0:
            elastic = elastic - e[0]
        newton = 0. if newtonian == 'none' else 1 / (4 * s2)
        return self._out((elastic + newton) * EarthConstant.ggg / EarthConstant.radiusm ** 2)

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

        return self._out((np.sum(item[:, 1:], axis=-1) + item2) * EarthConstant.radiusm / EarthConstant.Mass)


def demo():
    from pygwater.config import FileTool
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.PREM).get_Love_number()

    file = FileTool.get_input_dir('LLN/PREM-LGFs.dat')
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
