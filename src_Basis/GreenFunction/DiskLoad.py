"""
KERNEL of a single disk load (1 m EWH): response K(psi) at angular distance psi [deg] (paper Eq. 17).
    DiskLoad          : D-LGF of one radius
    DiskLoad_constrain: D-LGF set to zero beyond K * radius (approximation)
    DiskLoadTable     : D-LGF of many radii, tabulated in psi and interpolated (used for the CD-LGF)
The basis built from these kernels (many cells, any stations, design matrix) is GFBasis.py.
"""
import warnings

import numpy as np

from src_EarthModel.LoveNumber import LoveNumber, LLN_Data, LLN_variable, Frame
from src_Basis.GreenFunction.Legendre import Legendre_polynomial
from src_auxiliary.constants import EarthConstant


class DiskLoad:

    @staticmethod
    def required_lmax(radius):
        """
        Minimum Love-number / Legendre expansion for a disk of angular radius `radius` [degree], Eq. 24 of
        Yang et al. (2025): N = max(2 * 360 / radius, 10000).
        """
        return int(np.ceil(max(2 * 360.0 / np.min(radius), 10000)))

    def __init__(self, lln: LoveNumber, radius, thickness=1, check_lmax=True):
        """

        :param lln:
        :param radius: radius of the disk, unit: degree
        :param thickness: thickness of the disk in equivalent water height, unit: meter
        :param check_lmax: warn if lln.lmax violates the truncation rule of Eq. 24 (Yang et al., 2025)
        """
        self._residence = None
        self._scalar = False
        self._lln = lln

        if check_lmax and lln.lmax < self.required_lmax(radius):
            warnings.warn('D-LGF radius %.4f deg: lmax=%d is below the rule of Eq. 24 (N >= %d); '
                          'expect over-smoothing (artificial error) in the near field.'
                          % (np.min(radius), lln.lmax, self.required_lmax(radius)), stacklevel=2)

        '''compute sigma'''
        leg = Legendre_polynomial(alpha=radius, lmax=self._lln.lmax + 1)
        leg_back_2 = np.roll(leg, shift=2, axis=-1)

        sigma_n = leg_back_2 - leg
        sigma_n = np.delete(sigma_n, obj=1, axis=-1)
        sigma_n[0] = 1 - np.cos(np.deg2rad(radius))
        sigma_n *= 0.5 * thickness * EarthConstant.rhow
        self._sigma = sigma_n
        pass

    def configure(self, residence):
        """

        :param residence: deg, angular distance between two nodes (scalar or array)
        :return:
        """
        self._scalar = np.ndim(residence) == 0
        '''avoid the 1/sin(theta) singularity of the horizontal kernel at theta = 0 and 180 deg'''
        self._residence = np.clip(np.atleast_1d(np.array(residence, dtype=float)), 1e-6, 180 - 1e-6)
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
        n = np.arange(self._lln.lmax + 1)
        item = self._Pl * self._lln.LLN[LLN_variable.h] * self._sigma / (1 + 2 * n)

        return self._out(np.sum(item, axis=-1) * 3 / EarthConstant.rhoear)

    def getHorizental(self):
        """
        It could be generating NaN if theta (residence) is zero. Not the case for vertical displacement and geoid height
        :return:
        """
        n = np.arange(self._lln.lmax + 1)

        sinTheta = np.sin(np.deg2rad(self._residence))
        cosTheta = np.cos(np.deg2rad(self._residence))

        Plb = np.roll(self._Pl, 1, axis=-1)

        sinTheta = sinTheta[:, None]
        cosTheta = cosTheta[:, None]

        differential = 1 / sinTheta * (cosTheta * self._Pl - Plb) * n

        item = differential * self._lln.LLN[LLN_variable.l] * self._sigma / (1 + 2 * n)

        return self._out(np.sum(item, axis=-1) * 3 / EarthConstant.rhoear)

    def getGeoidheight(self):
        n = np.arange(self._lln.lmax + 1)
        item = self._Pl * (1 + self._lln.LLN[LLN_variable.k]) * self._sigma / (1 + 2 * n)

        return self._out(np.sum(item, axis=-1) * 3 / EarthConstant.rhoear)

    def getGravity(self, newtonian='below'):
        raise NotImplementedError(
            'Surface gravity is not implemented for the disk-load Green\'s function: the Newtonian term converges '
            'too slowly in the Legendre series and needs a closed form for the disk attraction. Use spherical '
            'harmonics, Slepian functions or the point-load LGF_Kummer.getGravity().')

    def release(self):
        """
        This is for reducing memory! Very important for global convolution
        :return:
        """
        self._residence = None
        self._Pl = None
        pass

class DiskLoad_constrain(DiskLoad):
    """
    When the angular distance is above K*radius, the impact of this load could be considered as negligible.
    """

    def __init__(self, lln: LoveNumber, radius, thickness=1, constrain_factor=10, check_lmax=True):
        """
        WARNING: the response beyond constrain_factor * radius is set to ZERO, not approximated. This truncates
        the long-wavelength (far-field) deformation and is NOT used in Yang et al. (2025), whose D-LGF results
        were computed without truncation (constrain_factor=1000).
        """
        super().__init__(lln, radius, thickness, check_lmax=check_lmax)
        self._residence_full = None
        self._threshold = constrain_factor * radius
        pass

    def configure(self, residence):
        """

        :param residence: deg, angular distance between two nodes
        :return:
        """
        self._scalar = False
        residence = np.clip(np.atleast_1d(np.array(residence, dtype=float)), 1e-6, 180 - 1e-6)
        self._residence_full = residence
        self._residence = residence[residence < self._threshold]
        self._legendre()

        return self

    def getVertical(self):
        a = super().getVertical()
        b = np.zeros_like(self._residence_full)
        b[self._residence_full < self._threshold] = a

        return b

    def getHorizental(self):
        a = super().getHorizental()
        b = np.zeros_like(self._residence_full)
        b[self._residence_full < self._threshold] = a

        return b

    def getGeoidheight(self):
        a=super().getGeoidheight()
        b = np.zeros_like(self._residence_full)
        b[self._residence_full < self._threshold] = a

        return b

    def release(self):
        """
        This is for reducing memory! Very important for global convolution
        :return:
        """
        super().release()
        self._residence_full = None
        pass


class DiskLoadTable:
    """
    D-LGF (Eq. 17 of Yang et al., 2025) of MANY disk radii, tabulated once on a dense grid of angular distances and
    evaluated by linear interpolation. This is what makes the clustered disk-load Green's function (CD-LGF) practical
    at arbitrary points: a 2-deg cell split into 0.1-deg sub-cells holds 400 disks, and evaluating each disk with
    its own Legendre series (DiskLoad) at every point costs O(N_disk * N_point * lmax).

    The Legendre polynomials P_n(cos psi) depend only on the distance psi, so for all radii r_k at once

        K_k(psi) = 3 / rho_e * sum_n  X_n * sigma_n(r_k) / (2n + 1) * P_n(cos psi)       (table = P @ C)

    with X_n = h_n (Vertical), 1 + k_n (Geoheight) and, for the horizontal kernel, P_n replaced by
    n (cos psi P_n - P_{n-1}) / sin psi and X_n = l_n (identical to DiskLoad.getHorizental).

    Distance grid: step r_min / 200 out to 10 r_max (disk interior and edge), min(0.0005 deg, r_min / 50) out to
    max(3 deg, 30 r_max), 0.005 deg to 20 deg and 0.05 deg to 180 deg. Interpolation error against the direct
    series (r = 0.056 deg, N = 13 000): below 1e-5 of the kernel peak. No far-field truncation is applied.
    Memory: (number of distances) x (number of radii) per variable, e.g. ~15 000 x 900 x 8 B = 110 MB for the
    distinct sub-cell radii of a global 0.1-deg partition.
    """

    def __init__(self, lln: LoveNumber, radii, thickness=1, check_lmax=True, chunk=500):
        """
        :param radii: disk radii [degree], any shape; duplicates are merged
        :param thickness: EWH thickness of the disks [m]
        :param chunk: number of distances per Legendre block (memory: chunk x lmax floats)
        """
        self._lln = lln
        self.radii = np.unique(np.round(np.atleast_1d(np.asarray(radii, dtype=float)), 10))
        self._chunk = chunk
        if check_lmax and lln.lmax < DiskLoad.required_lmax(self.radii):
            warnings.warn('CD-LGF smallest radius %.4f deg: lmax=%d is below the rule of Eq. 24 (N >= %d); '
                          'expect over-smoothing (artificial error) in the near field.'
                          % (self.radii.min(), lln.lmax, DiskLoad.required_lmax(self.radii)), stacklevel=2)

        # sigma_n for all radii at once (same expression as DiskLoad), stored as (lmax+1, n_radius)
        leg = Legendre_polynomial(alpha=self.radii, lmax=lln.lmax + 1)
        sigma = np.roll(leg, shift=2, axis=-1) - leg
        sigma = np.delete(sigma, obj=1, axis=-1)
        sigma[:, 0] = 1 - np.cos(np.deg2rad(self.radii))
        self._sigma = (sigma * 0.5 * thickness * EarthConstant.rhow).T

        rmin, rmax = self.radii.min(), self.radii.max()
        edge_end, near_end = 10 * rmax, max(3.0, 30 * rmax)
        psi = np.concatenate([np.arange(0, edge_end, rmin / 200),
                              np.arange(edge_end, near_end, min(5e-4, rmin / 50)), np.arange(near_end, 20, 5e-3),
                              np.arange(max(near_end, 20), 180, 5e-2), [180.0]])
        self.psi = np.unique(np.clip(psi, 1e-6, 180 - 1e-6))
        self._tables = {}

    def _build(self, kind):
        n = np.arange(self._lln.lmax + 1)
        L = self._lln.LLN
        X = {'vertical': L[LLN_variable.h], 'geoid': 1 + L[LLN_variable.k], 'horizontal': L[LLN_variable.l]}[kind]
        C = (X * 3 / EarthConstant.rhoear / (2 * n + 1))[:, None] * self._sigma  # (lmax+1, n_radius)
        out = np.empty((self.psi.size, self.radii.size))
        for i in range(0, self.psi.size, self._chunk):
            psi = self.psi[i:i + self._chunk]
            P = Legendre_polynomial(alpha=psi, lmax=self._lln.lmax)
            if kind == 'horizontal':
                cos_t, sin_t = np.cos(np.deg2rad(psi))[:, None], np.sin(np.deg2rad(psi))[:, None]
                P = n * (cos_t * P - np.roll(P, 1, axis=-1)) / sin_t
            out[i:i + self._chunk] = P @ C
        return out

    def table(self, kind):
        """kernel table (n_psi, n_radius); kind = 'vertical', 'geoid' or 'horizontal' (away from the load)"""
        if kind not in self._tables:
            self._tables[kind] = self._build(kind)
        return self._tables[kind]

    def radius_index(self, radius):
        """index of each radius in the table"""
        r = np.round(np.asarray(radius, dtype=float), 10)
        idx = np.searchsorted(self.radii, r)
        if not np.all(self.radii[np.minimum(idx, self.radii.size - 1)] == r):
            raise ValueError('radius not tabulated')
        return idx

    def evaluate(self, kind, radius, psi):
        """kernel of the disk of the given radius (scalar) at distance(s) psi [degree], any shape"""
        k = int(self.radius_index(radius))
        return np.interp(np.clip(psi, 1e-6, 180 - 1e-6), self.psi, self.table(kind)[:, k])


def demo():
    from src_auxiliary.FileTool import FileTool
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.REF).get_Love_number()

    # file = FileTool.get_data_dir('LLN/PREM-LGFs.dat')
    # theta = np.loadtxt(file, skiprows=1, usecols=0)
    # theta = 0.001100220044009
    theta = np.array([0.001000200040008, 0.001100220044009])
    load = DiskLoad(lln=lln, radius=0.10).configure(residence=theta)
    u = load.getVertical()
    v = load.getHorizental()
    g = load.getGeoidheight()
    print(u)
    print(v)
    print(g)

    load = DiskLoad_constrain(lln=lln, radius=0.10, constrain_factor=0.01).configure(residence=theta)
    u = load.getVertical()
    v = load.getHorizental()
    g = load.getGeoidheight()
    print(u)
    print(v)
    print(g)
    pass


if __name__ == '__main__':
    demo()
