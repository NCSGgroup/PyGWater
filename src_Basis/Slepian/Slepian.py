"""
Slepian basis for regional modelling of load-induced deformation (Sect. 8 of Yang et al., 2025, JGR Solid Earth).

Theory
------
Band-limited (degree <= L) real, 4pi-normalised spherical harmonics Y_k, k = 1..(L+1)^2.
Concentration matrix of a region R (Simons et al., 2006; Albertella et al., 1999):

        D_kk' = 1/(4 pi) * int_R  Y_k Y_k'  dOmega

D is symmetric with eigenvalues 0 <= beta_i <= 1 (Eq. 28). Its eigenvectors g_i (columns of G) are the
spherical-harmonic coefficients of the Slepian functions (Eq. 27). sum(beta) = trace(D) is the Shannon number
N = (L+1)^2 * A_R / (4 pi).

Forward modelling (Eqs. 29-31). With the load given as EWH spherical-harmonic coefficients E_k (metre), the
vertical displacement coefficients are (Eq. 2)

        U_k = 3 rho_w / rho_e * h_n / (2n + 1) * E_k                                   (n = degree of k)

and its Slepian coefficients and truncated representation are

        S_i = sum_k U_k g_ik          (Eq. 31, here in metres, i.e. the factor a of Eq. 26 is absorbed)
        U(Omega) ~ sum_{i<=Lc} S_i g_i(Omega)                                           (Eq. 29)

With Lc = (L+1)^2 the Slepian solution equals the spherical-harmonic one exactly (Eq. 30); truncation to Lc
introduces the "approximation" error discussed in Sect. 8.3.

Other observables (argument `variable`, EnumClasses.Displacement) use the same scalar machinery with another
Love-number factor (Farrell, 1972; Wang et al., 2012):
    Vertical  : 3 rho_w / rho_e * h_n / (2n+1)            evaluated with Y
    Geoheight : 3 rho_w / rho_e * (1 + k_n) / (2n+1)      evaluated with Y
    North     : 3 rho_w / rho_e * l_n / (2n+1)            evaluated with -dY/dtheta            (theta: colatitude)
    East      : 3 rho_w / rho_e * l_n / (2n+1)            evaluated with  1/sin(theta) dY/dlambda
For North/East, S_i are the Slepian coefficients of the horizontal displacement POTENTIAL, whose surface gradient
is the horizontal displacement vector.

Inversion
---------
If the regional LOAD is parameterised by Slepian coefficients x (EWH(Omega) = sum_i x_i g_i(Omega)), displacements
at arbitrary points (e.g. GNSS stations) are  u = A x  with

        A = Y(points) @ diag(3 rho_w / rho_e * h_n / (2n+1)) @ G[:, :Lc]

see SlepianBasis.design_matrix(). This is the regional inversion operator used later by PyGWater.

Coefficient ordering
--------------------
Vectors of length (L+1)^2 are ordered by degree n, then cosine terms m = 0..n, then sine terms m = 1..n.
Use SlepianBasis.cs2d_to_vec (inherited from SphericalHarmonicBasis) / vec_to_cs2d to convert from/to the (cqlm, sqlm) arrays used by Harmonic/SHC.
"""

import h5py
import numpy as np

from src_auxiliary.EnumClasses import Displacement
from src_EarthModel.LoveNumber import LoveNumber
from src_Basis.SphericalHarmonic.SHBasis import SphericalHarmonicBasis


class SlepianBasis(SphericalHarmonicBasis):
    name = 'slepian'

    def __init__(self, lmax: int):
        """
        :param lmax: band limit L of the Slepian functions (spherical-harmonic degree)
        """
        super().__init__(lmax)
        self.D = None  # concentration matrix
        self.beta = None  # eigenvalues, descending
        self.G = None  # eigenvectors (columns), SH coefficients of Slepian functions
        self.Lc = None  # number of retained Slepian functions
        self.area = None  # area of the region / (4 pi), i.e. fraction of the sphere

    @property
    def n_param(self):
        return self.Lc

    def _coefficient_map(self):
        return self.G[:, :self.Lc]

    # ------------------------------------------------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------------------------------------------------
    def configure(self, mask, lat, lon, resolution=None):
        """
        Build the concentration matrix D (Eq. 28) from a regional mask on a regular, cell-centred lat/lon grid,
        e.g. the output of LoadMesh.ShapeAndMask.basinMask ('basin', 'lat', 'lon', 'resolution').

        The integral over R is evaluated cell by cell: D = 1/(4pi) sum_cells w Y Y^T dOmega, with dOmega the exact
        cell area cos(lat) dlat dlon. Use a mask resolution clearly finer than 180/L degrees (e.g. <= 0.5 deg
        for L = 60) so that the region boundary, not the quadrature, controls D.

        :param mask: 2d array, 1/True inside the region (fractional weights in [0, 1] are allowed)
        :param lat: 2d (meshgrid) or 1d latitudes of cell centres [degree]
        :param lon: 2d (meshgrid) or 1d longitudes of cell centres [degree]
        :param resolution: grid spacing [degree]; inferred from lat/lon if None
        :return: self
        """
        mask = np.asarray(mask, dtype=float)
        lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
        if lat.ndim == 1 and lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        assert lat.shape == mask.shape == lon.shape

        if resolution is None:
            resolution = np.abs(np.diff(np.unique(lat))).min()
        res = np.deg2rad(resolution)

        inside = mask > 0
        if not inside.any():
            raise ValueError('The mask is empty.')

        rows = np.unique(np.nonzero(inside)[0])
        D = np.zeros((self.ncoef, self.ncoef))
        area = 0.
        for r in rows:
            sel = inside[r]
            w = mask[r, sel] * np.cos(np.deg2rad(lat[r, sel])) * res * res / (4 * np.pi)
            Y = self.Ymatrix(lat[r, sel], lon[r, sel])
            D += (Y * w[:, None]).T @ Y
            area += w.sum()

        self.D = 0.5 * (D + D.T)
        self.area = area

        beta, G = np.linalg.eigh(self.D)
        order = np.argsort(beta)[::-1]
        self.beta, self.G = beta[order], G[:, order]
        self.Lc = self.ncoef
        return self

    @classmethod
    def from_mask_h5(cls, lmax, h5file, dataset='basin'):
        """Convenience: build from a mask file written by LoadMesh.ShapeAndMask.basinMask."""
        with h5py.File(h5file, 'r') as f:
            mask, lat, lon = f[dataset][()], f['lat'][()], f['lon'][()]
            res = float(f['resolution'][()]) if 'resolution' in f else None
        return cls(lmax).configure(mask=mask, lat=lat, lon=lon, resolution=res)

    @property
    def shannon(self):
        """Shannon number N = (L+1)^2 A / (4 pi) = trace(D)"""
        return float(np.trace(self.D))

    def truncate(self, Lc=None, beta_min=None, shannon=False):
        """
        Choose the number of retained Slepian functions Lc (Eq. 29). Exactly one option:
        :param Lc: explicit number
        :param beta_min: keep functions with concentration beta > beta_min, e.g. 0.1 (Jiang et al., 2021) or
                         0.01 (Shen et al., 2024)
        :param shannon: keep round(Shannon number) functions (Simons et al., 2006)
        :return: self
        """
        assert sum(x is not None and x is not False for x in (Lc, beta_min, shannon)) == 1
        if Lc is not None:
            self.Lc = int(Lc)
        elif beta_min is not None:
            self.Lc = int(np.sum(self.beta > beta_min))
        else:
            self.Lc = int(round(self.shannon))
        self.Lc = max(1, min(self.Lc, self.ncoef))
        return self

    # ------------------------------------------------------------------------------------------------------------
    # forward modelling (Eqs. 29-31)
    # ------------------------------------------------------------------------------------------------------------
    def slepian_coefficients(self, vec):
        """Project an SH coefficient vector (or matrix (ntime, ncoef)) onto the retained Slepian functions."""
        return np.asarray(vec) @ self.G[:, :self.Lc]

    def displacement_coefficients(self, load_cqlm, load_sqlm, lln: LoveNumber, variable=Displacement.Vertical,
                                  **kwargs):
        """
        Eq. 31: Slepian coefficients S_i of the observable [m], given the LOAD as EWH spherical-harmonic
        coefficients [m] (e.g. from GRID.analysis(lmax=..., from_type=EWH, to_type=EWH)).
        For North/East the coefficients describe the horizontal displacement potential (same for both).
        :return: S, shape (Lc,) or (ntime, Lc)
        """
        U = self.cs2d_to_vec(load_cqlm, load_sqlm) * self._displacement_factor(lln, variable, **kwargs)
        return self.slepian_coefficients(U)

    def synthesis(self, S, lat, lon, component='value'):
        """
        Eq. 29: evaluate sum_i S_i g_i (component='value') or its surface gradient ('north', 'east') at points.
        :param S: (Lc,) or (ntime, Lc)
        :param lat, lon: 1d arrays of points [degree]
        :return: (N,) or (ntime, N)
        """
        return self.synthesis_vec(np.asarray(S) @ self.G[:, :self.Lc].T, lat, lon, component)

    def forward(self, load_cqlm, load_sqlm, lln: LoveNumber, lat, lon, variable=Displacement.Vertical, **kwargs):
        """
        Observable at points (lat, lon) from the truncated Slepian representation (m, or m/s^2 for gravity).
        :param variable: any observable of src_EarthModel/Observables.py: Vertical (default), North, East, Geoheight,
                         Gravity (kwarg newtonian='below'|'above'|'none'), GravityDisturbance (kwarg height=[m])
        """
        S = self.displacement_coefficients(load_cqlm, load_sqlm, lln, variable, **kwargs)
        return self.synthesis(S, lat, lon, self._component(variable))

    # ------------------------------------------------------------------------------------------------------------
    # inversion operator: design_matrix, design_matrix_los, predict and load_synthesis are inherited from
    # SphericalHarmonicBasis / Basis with the coefficient map G[:, :Lc], i.e.
    #     A = Y(points) @ diag(F_n X_n) @ G[:, :Lc],   EWH(Omega) = sum_i x_i g_i(Omega)
    # ------------------------------------------------------------------------------------------------------------
    def load_parameters(self, ewh, lat, lon, method='auto'):
        """x of a gridded EWH load [m]: SH analysis, then projection onto the retained Slepian functions."""
        return self.slepian_coefficients(self.sh_coefficients(ewh, lat, lon, method=method))


def demo():
    """Spherical cap check: Shannon number of a cap of radius 30 deg is (L+1)^2 (1 - cos 30)/2."""
    from src_auxiliary.GeoMathKit import GeoMathKit
    L, res, cap = 20, 0.5, 30.
    lat, lon = GeoMathKit.global_equal_distance_grid(res)
    lon, lat = np.meshgrid(lon, lat)
    sb = SlepianBasis(L).configure(mask=(lat > 90 - cap), lat=lat, lon=lon, resolution=res)
    print('Shannon number:', sb.shannon, ' theory:', (L + 1) ** 2 * (1 - np.cos(np.deg2rad(cap))) / 2)
    print('largest eigenvalues:', np.round(sb.beta[:5], 4))


if __name__ == '__main__':
    demo()
