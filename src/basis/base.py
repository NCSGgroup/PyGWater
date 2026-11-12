"""
Common interface of all basis functions in PyGWater.

Every basis is a LINEAR PARAMETERISATION of the surface load (equivalent water height, EWH [m]),

        EWH(Omega) = sum_j x_j phi_j(Omega),

together with the linear operator that maps the parameters x to any observable at arbitrary points (GNSS stations,
InSAR pixels, grid nodes):

        d = A(points, observable, Earth model) @ x.

The same three methods therefore serve forward modelling (x known -> d) and inversion (d observed -> x):

    design_matrix(lat, lon, lln, variable, **kw)  ->  A, shape (N_points, n_param)
    predict(x, lat, lon, lln, variable, **kw)      ->  d = A x, shape (N_points,) or (ntime, N_points)
    load_synthesis(x, lat, lon)                    ->  EWH [m] at points
    load_parameters(ewh, lat, lon)                 ->  x from a gridded EWH load (forward-modelling input)

Implementations
---------------
    pygwater.basis.spherical_harmonic.sh_basis.SphericalHarmonicBasis   x = SH coefficients of EWH (global, degree <= L)
    pygwater.basis.slepian.SlepianBasis                       x = Slepian coefficients of EWH (regional)
    pygwater.basis.green.gf_basis.GreenFunctionBasis           x = EWH of point masses or disks (cells)

Observables (argument `variable`, pygwater.utils.enums.Displacement) are defined once in
src/earth/observables.py; InSAR line of sight is built from East, North and Vertical (design_matrix_los).
"""

from abc import ABC, abstractmethod

import numpy as np

from pygwater.utils.enums import Displacement


class Basis(ABC):
    """Abstract base class: a linear load parameterisation and its response operator."""

    name = 'basis'

    @property
    @abstractmethod
    def n_param(self) -> int:
        """Number of parameters x."""

    @abstractmethod
    def design_matrix(self, lat, lon, lln, variable=Displacement.Vertical, **kwargs) -> np.ndarray:
        """
        A (N x n_param): observable at the points per unit parameter (unit: observable per metre EWH).
        :param lat, lon: points [degree], 1d arrays of equal length
        :param lln: LoveNumber object (Earth model and reference frame)
        :param variable: observable, enums.Displacement (Vertical, North, East, Geoheight, Gravity, ...)
        :param kwargs: observable options, e.g. newtonian='above' (Gravity), height=500e3 (GravityDisturbance)
        """

    @abstractmethod
    def load_synthesis(self, x, lat, lon) -> np.ndarray:
        """EWH [m] at the points from the parameters x, shape (n_param,) or (ntime, n_param)."""

    @abstractmethod
    def load_parameters(self, ewh, lat, lon) -> np.ndarray:
        """
        Parameters x of a gridded EWH load [m].
        :param ewh: (nlat, nlon) or (ntime, nlat, nlon) on the regular cell-centred grid (lat, lon)
        :param lat, lon: 1d grid coordinates [degree]
        :return: (n_param,) or (ntime, n_param)
        """

    # ------------------------------------------------------------------------------------------------------------
    # generic operators built on design_matrix
    # ------------------------------------------------------------------------------------------------------------
    def predict(self, x, lat, lon, lln, variable=Displacement.Vertical, **kwargs) -> np.ndarray:
        """d = A x at the points; x of shape (n_param,) or (ntime, n_param) -> (N,) or (ntime, N)."""
        return np.asarray(x) @ self.design_matrix(lat, lon, lln, variable, **kwargs).T

    def design_matrix_los(self, lat, lon, lln, look) -> np.ndarray:
        """
        A (N x n_param) for InSAR line-of-sight displacement (positive towards the satellite).
        :param look: ground-to-satellite unit vector(s) (east, north, up), (3,) or (3, N); see Observables.look_vector
        """
        look = np.asarray(look, dtype=float)
        le, ln_, lu = (look[i][:, None] if look.ndim == 2 else look[i] for i in range(3))
        return (le * self.design_matrix(lat, lon, lln, Displacement.East)
                + ln_ * self.design_matrix(lat, lon, lln, Displacement.North)
                + lu * self.design_matrix(lat, lon, lln, Displacement.Vertical))

    def __repr__(self):
        return f'{type(self).__name__}(n_param={self.n_param})'
