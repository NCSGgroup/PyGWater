"""
Basis-independent forward model: load (EWH) -> observables at points, for any Basis of pygwater.basis.

    from pygwater.forward.model import ForwardModel
    fm = ForwardModel(basis, lln)                       # basis: SphericalHarmonicBasis, SlepianBasis, GreenFunctionBasis
    x = fm.parameterize(ewh, lat, lon)                  # gridded EWH [m] (nlat, nlon) or (ntime, nlat, nlon) -> x
    out = fm.run(x, sta_lat, sta_lon, observables=(Displacement.Vertical, Displacement.North, Displacement.East))
    A = fm.design_matrix(sta_lat, sta_lon, observables=(...))     # stacked operator for inversion

The convolution classes in disk.py, point_mass.py and disk_cluster.py remain available for large global grids and the
cluster-disk templates; ForwardModel is the common entry point used by the flow control and the inversion.
"""

import numpy as np

from pygwater.utils.enums import Displacement
from pygwater.basis.base import Basis


class ForwardModel:

    def __init__(self, basis: Basis, lln):
        """
        :param basis: any pygwater.basis.base.Basis
        :param lln: LoveNumber object (Earth model, frame, maximum degree)
        """
        self.basis = basis
        self.lln = lln

    def parameterize(self, ewh, lat, lon, **kwargs):
        """basis parameters x of a gridded EWH load [m] (see Basis.load_parameters)"""
        return self.basis.load_parameters(ewh, lat, lon, **kwargs)

    def design_matrix(self, lat, lon, observables=(Displacement.Vertical,), look=None, **kwargs):
        """
        Operator stacked over observables, rows ordered observable by observable: shape (n_obs * N, n_param).
        :param observables: sequence of enums.Displacement, or the string 'LOS' (needs look)
        :param look: InSAR ground-to-satellite unit vector(s), see Observables.look_vector
        :param kwargs: observable options (newtonian=..., height=...)
        """
        blocks = []
        for obs in observables:
            if isinstance(obs, str) and obs.upper() == 'LOS':
                blocks.append(self.basis.design_matrix_los(lat, lon, self.lln, look))
            else:
                blocks.append(self.basis.design_matrix(lat, lon, self.lln, obs, **kwargs))
        return np.concatenate(blocks, axis=0)

    def run(self, x, lat, lon, observables=(Displacement.Vertical,), look=None, **kwargs):
        """
        Observables at the points from the parameters x ((n_param,) or (ntime, n_param)).
        :return: dict {observable: (N,) or (ntime, N)}; units m (displacement, geoid) or m/s^2 (gravity)
        """
        x = np.asarray(x)
        out = {}
        for obs in observables:
            A = self.design_matrix(lat, lon, (obs,), look=look, **kwargs)
            out[obs] = x @ A.T
        return out

    def load(self, x, lat, lon):
        """EWH [m] represented by the parameters at the points"""
        return self.basis.load_synthesis(x, lat, lon)
