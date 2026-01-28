import numpy as np

from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.Legendre import Legendre_polynomial
from src.Auxiliary.constants import EarthConstant


class DiskLoad:

    def __init__(self, lln: LoveNumber, radius, thickness=1):
        """

        :param lln:
        :param radius: radius of the disk, unit: degree
        """
        self._residence = None
        self._lln = lln

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

        :param residence: deg, angular distance between two nodes
        :return:
        """
        self._residence = residence
        self._legendre()

        return self

    def _legendre(self):
        self._Pl = Legendre_polynomial(alpha=self._residence, lmax=self._lln.lmax)
        pass

    def getVertical(self):
        n = np.arange(self._lln.lmax + 1)
        item = self._Pl * self._lln.LLN[LLN_variable.h] * self._sigma / (1 + 2 * n)

        return np.sum(item, axis=-1) * 3 / EarthConstant.rhoear

    def getHorizental(self):
        """
        It could be generating NaN if theta (residence) is zero. Not the case for vertical displacement and geoid height
        :return:
        """
        n = np.arange(self._lln.lmax + 1)

        sinTheta = np.sin(np.deg2rad(self._residence))
        cosTheta = np.cos(np.deg2rad(self._residence))

        Plb = np.roll(self._Pl, 1, axis=-1)

        if len(np.shape(sinTheta)) != 0:
            sinTheta = sinTheta[:, None]
            cosTheta = cosTheta[:, None]

        differential = 1 / sinTheta * (cosTheta * self._Pl - Plb) * n

        item = differential * self._lln.LLN[LLN_variable.l] * self._sigma / (1 + 2 * n)

        return np.sum(item, axis=-1) * 3 / EarthConstant.rhoear

    def getGeoidheight(self):
        n = np.arange(self._lln.lmax + 1)
        item = self._Pl * (1 + self._lln.LLN[LLN_variable.k]) * self._sigma / (1 + 2 * n)

        return np.sum(item, axis=-1) * 3 / EarthConstant.rhoear

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

    def __init__(self, lln: LoveNumber, radius, thickness=1, constrain_factor=10):
        super().__init__(lln, radius, thickness)
        self._residence_full = None
        self._threshold = constrain_factor * radius
        pass

    def configure(self, residence):
        """

        :param residence: deg, angular distance between two nodes
        :return:
        """
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

def demo():
    from src.Auxiliary.FileTool import FileTool
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.REF).get_Love_number()

    # file = FileTool.get_project_dir() / 'data' / 'LLN' / 'PREM-LGFs.dat'
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
