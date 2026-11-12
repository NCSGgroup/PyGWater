"""
PyGWater: elastic loading deformation of the Earth (GNSS displacement, geoid, gravity) from surface water loads,
and the inversion of GNSS observations for those loads, with a choice of basis functions (Yang et al., 2025).

The code lives in the folder src/ of the repository and is imported as `pygwater` once installed
(`pip install -e .` in the repository folder, see doc/installation.md). The main classes are available directly:

    import pygwater as pgw
    lln = pgw.LoveNumber().config(lmax=13000, method=pgw.LLN_Data.Wang).get_Love_number()
    basis = pgw.GreenFunctionBasis.from_grid(lat, lon, 0.25, kernel='cluster_disk')
    A = basis.design_matrix(station_lat, station_lon, lln, pgw.Displacement.Vertical)
    sol = pgw.RegularizedLeastSquares(A, sigma, R).solve(d, alpha)

Subpackages: earth (Love numbers, observables), basis (spherical harmonics, Green's functions, Slepian),
forward, inversion, observations (GNSS, synthetic), mesh (grids, masks), utils; config holds the data root.
"""
__version__ = '0.9.0'

from pygwater.config import FileTool
from pygwater.utils.enums import Displacement, LLN_Data, LLN_variable, Frame
from pygwater.utils.constants import EarthConstant
from pygwater.earth.love_numbers import LoveNumber
from pygwater.basis.base import Basis
from pygwater.basis.spherical_harmonic.sh_basis import SphericalHarmonicBasis
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.basis.slepian import SlepianBasis
from pygwater.forward.model import ForwardModel
from pygwater.inversion.solver import RegularizedLeastSquares
from pygwater.inversion import regularization
from pygwater.preprocessing import GNSSNetwork, GNSSDataset, read_tenv3, read_ngl_steps
from pygwater.inversion.workflow import MonthlyInversion, InversionResult
from pygwater.observations.grace import Mascon

__all__ = ['FileTool', 'Displacement', 'LLN_Data', 'LLN_variable', 'Frame', 'EarthConstant', 'LoveNumber', 'Basis',
           'SphericalHarmonicBasis', 'GreenFunctionBasis', 'SlepianBasis', 'ForwardModel', 'RegularizedLeastSquares',
           'regularization', 'GNSSNetwork', 'GNSSDataset', 'read_tenv3', 'read_ngl_steps', 'MonthlyInversion', 'InversionResult', 'Mascon',
           '__version__']
