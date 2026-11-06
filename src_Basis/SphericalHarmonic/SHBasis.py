"""
Spherical harmonics as a Basis (src_Basis/base.py): x = 4pi-normalised SH coefficients of the load EWH [m] up to
degree L, in vector form. This is the point-wise counterpart of the grid synthesis in Harmonic.py / DataClass.py
and the parent class of the Slepian basis.

        d(points) = Y(points) @ diag(F_n X_n) @ x        (Eq. 2 of Yang et al., 2025, generalised to any observable)

with F_n X_n the Love-number factor of the observable (src_EarthModel/Observables.love_factor).

Coefficient ordering
--------------------
Vectors of length (L+1)^2 are ordered by degree n, then cosine terms m = 0..n, then sine terms m = 1..n.
Use cs2d_to_vec / vec_to_cs2d to convert from/to the (cqlm, sqlm) arrays used by Harmonic/SHC.
"""

import numpy as np

from src_auxiliary.MathTool import MathTool
from src_auxiliary.EnumClasses import Displacement
from src_EarthModel import Observables
from src_Basis.base import Basis


class SphericalHarmonicBasis(Basis):
    name = 'spherical_harmonics'

    def __init__(self, lmax: int):
        """
        :param lmax: maximum degree L
        """
        self.lmax = int(lmax)
        L = self.lmax

        '''index map of the coefficient vector'''
        n_list, m_list, cs_list = [], [], []
        for n in range(L + 1):
            for m in range(n + 1):
                n_list.append(n), m_list.append(m), cs_list.append(0)  # cosine
            for m in range(1, n + 1):
                n_list.append(n), m_list.append(m), cs_list.append(1)  # sine
        self.idx_n = np.array(n_list)
        self.idx_m = np.array(m_list)
        self.idx_cs = np.array(cs_list)
        self.ncoef = len(n_list)
        assert self.ncoef == (L + 1) ** 2

    @property
    def n_param(self):
        return self.ncoef

    # ------------------------------------------------------------------------------------------------------------
    # spherical harmonics in vector form
    # ------------------------------------------------------------------------------------------------------------
    def cs2d_to_vec(self, cqlm, sqlm):
        """(cqlm, sqlm) arrays of shape (L'+1, L'+1) or (ntime, L'+1, L'+1), L' >= L  ->  vector(s) (..., (L+1)^2)"""
        cqlm, sqlm = np.asarray(cqlm), np.asarray(sqlm)
        c = cqlm[..., self.idx_n, self.idx_m]
        s = sqlm[..., self.idx_n, self.idx_m]
        return np.where(self.idx_cs == 0, c, s)

    def vec_to_cs2d(self, vec):
        """vector(s) (..., (L+1)^2) -> (cqlm, sqlm) arrays (..., L+1, L+1)"""
        vec = np.asarray(vec)
        shape = vec.shape[:-1] + (self.lmax + 1, self.lmax + 1)
        cqlm, sqlm = np.zeros(shape), np.zeros(shape)
        ic, is_ = self.idx_cs == 0, self.idx_cs == 1
        cqlm[..., self.idx_n[ic], self.idx_m[ic]] = vec[..., ic]
        sqlm[..., self.idx_n[is_], self.idx_m[is_]] = vec[..., is_]
        return cqlm, sqlm

    def Ymatrix(self, lat, lon, component='value'):
        """
        Real 4pi-normalised spherical harmonics (or their surface gradient) at points.
        :param lat: geographic latitude [degree], 1d array (N,)
        :param lon: longitude [degree], 1d array (N,)
        :param component: 'value' -> Y;  'north' -> -dY/dtheta;  'east' -> 1/sin(theta) dY/dlambda
        :return: (N, (L+1)^2)
        """
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        mlon = np.deg2rad(lon)[:, None] * self.idx_m[None, :]
        cos_m, sin_m = np.cos(mlon), np.sin(mlon)
        is_cos = self.idx_cs == 0

        '''Legendre functions only for the distinct latitudes (grids share latitudes), then broadcast'''
        ulat, inv = np.unique(lat, return_inverse=True)

        if component == 'value':
            P = MathTool.get_Legendre(ulat, self.lmax, option=1)[:, self.idx_n, self.idx_m][inv]
            return P * np.where(is_cos, cos_m, sin_m)
        elif component == 'north':
            dP = MathTool.get_Legendre_derivative(ulat, self.lmax, option=1)[:, self.idx_n, self.idx_m][inv]
            return -dP * np.where(is_cos, cos_m, sin_m)
        elif component == 'east':
            P = MathTool.get_Legendre(ulat, self.lmax, option=1)[:, self.idx_n, self.idx_m][inv]
            sin_theta = np.cos(np.deg2rad(lat))[:, None]  # sin(colatitude)
            return P * self.idx_m[None, :] * np.where(is_cos, -sin_m, cos_m) / sin_theta
        else:
            raise ValueError(component)

    def _blocks(self, n, max_bytes=2e8):
        """split n points into blocks so that a (block x ncoef) float64 matrix stays below max_bytes"""
        size = max(1, int(max_bytes / (8 * self.ncoef)))
        return [slice(i, min(i + size, n)) for i in range(0, n, size)]

    @staticmethod
    def _component(variable):
        """'value', 'north' or 'east', from the observable table in src_EarthModel/Observables.py"""
        return Observables.kind(variable)

    def _displacement_factor(self, lln, variable=Displacement.Vertical, **kwargs):
        """
        Per-coefficient factor from load (EWH) to the observable (Eq. 2 / Eq. 30), taken from the observable table
        (src_EarthModel/Observables.love_factor), e.g. 3 rho_w / rho_e * h_n / (2n+1) for Vertical.
        :param kwargs: newtonian=... (Gravity), height=... (GravityDisturbance)
        """
        return Observables.love_factor(variable, lln, self.idx_n, **kwargs)

    # ------------------------------------------------------------------------------------------------------------
    # Basis interface
    # ------------------------------------------------------------------------------------------------------------
    def _coefficient_map(self):
        """Columns of the parameterisation in terms of SH coefficients; identity here, G[:, :Lc] for Slepian."""
        return None

    def design_matrix(self, lat, lon, lln, variable=Displacement.Vertical, **kwargs):
        """
        A (N x n_param): observable at points per unit parameter of the load (EWH, metre).
        :param variable: any observable of src_EarthModel/Observables.py; kwargs newtonian=..., height=...
        """
        lat, lon = np.atleast_1d(lat), np.atleast_1d(lon)
        f = self._displacement_factor(lln, variable, **kwargs)[None, :]
        Gc = self._coefficient_map()
        comp = self._component(variable)
        blocks = []
        for b in self._blocks(len(lat)):
            Yf = self.Ymatrix(lat[b], lon[b], comp) * f
            blocks.append(Yf if Gc is None else Yf @ Gc)
        return np.concatenate(blocks, axis=0)

    def synthesis_vec(self, vec, lat, lon, component='value'):
        """Evaluate SH coefficient vector(s) vec (..., (L+1)^2) or their surface gradient at points."""
        vec = np.asarray(vec)
        lat, lon = np.atleast_1d(lat), np.atleast_1d(lon)
        return np.concatenate([vec @ self.Ymatrix(lat[b], lon[b], component).T for b in self._blocks(len(lat))],
                              axis=-1)

    def load_synthesis(self, x, lat, lon):
        """EWH [m] at points from the parameters x."""
        Gc = self._coefficient_map()
        x = np.asarray(x)
        return self.synthesis_vec(x if Gc is None else x @ Gc.T, lat, lon)

    def sh_coefficients(self, ewh, lat, lon, method='auto'):
        """SH coefficient vector(s) of a gridded EWH load [m] (analysis as in Harmonic.analysis, SaGEA-style)."""
        from src_Basis.SphericalHarmonic.Harmonic import Harmonic
        ewh = np.asarray(ewh, dtype=float)
        cqlm, sqlm = Harmonic(np.asarray(lat), np.asarray(lon), self.lmax, option=1).analysis(ewh, method=method)
        vec = self.cs2d_to_vec(cqlm, sqlm)
        return vec[0] if ewh.ndim == 2 and vec.ndim == 2 else vec

    def load_parameters(self, ewh, lat, lon, method='auto'):
        """x of a gridded EWH load [m]: spherical-harmonic analysis ('auto', 'dh', 'midpoint', 'least_squares')."""
        return self.sh_coefficients(ewh, lat, lon, method=method)
