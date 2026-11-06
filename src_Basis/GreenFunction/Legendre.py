"""Legendre polynomials P_n(cos psi), used by the kernels in PointLoad.py and DiskLoad.py."""
import pickle

import numpy as np


def Legendre_polynomial(alpha, lmax: int):
    """
    get legendre polynomial up to degree/order lmax in Lat. This is traditional unnormalized Legendre function.
    :param alpha: ndarray, relative angle between two spherical nodes, unit[degree]
    :param lmax: int, max degree
    :return: 3d-ndarray, indexes stand for (len(alpha), degree l)
    """
    z = np.cos(np.deg2rad(alpha))

    P0 = np.ones_like(z)
    P1 = np.ones_like(z) * z
    P = [P0, P1]
    for l in range(2, lmax + 1):
        Pl = (2 * l - 1) * z * P[-1] / l - (l - 1) * P[-2] / l
        P.append(Pl)

    return np.array(P).T


def demo():
    P = Legendre_polynomial(alpha=np.ones(5000)*60, lmax=10000)
    pass

if __name__ == '__main__':
    demo()