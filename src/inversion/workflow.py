"""
Monthly inversion of a GNSS network for the surface load, with any basis:

    from pygwater.inversion.workflow import MonthlyInversion
    inv = MonthlyInversion(basis, lln, dataset, components=('U',), R=regularization.curvature(lat, lon, res))
    alpha = inv.choose_alpha(method='gcv')                    # one alpha for all months (median of monthly choices)
    res = inv.run(alpha, lb=0.0)                              # InversionResult: x per month
    s = res.regional_series(lat, lon, inside)                 # area-weighted mean EWH over a region [m]

The design matrix of every station is computed once; each month uses the rows of the stations with data in that
month. Observation weights are, by default, one sigma per station and component (the median of its monthly sigmas),
so months sharing the same stations share one decomposition of the normal matrix and the whole series is fast.
"""
from dataclasses import dataclass, field

import numpy as np

from pygwater.utils.enums import Displacement
from pygwater.utils.constants import EarthConstant
from pygwater.inversion.solver import RegularizedLeastSquares

_OBS = {'U': Displacement.Vertical, 'N': Displacement.North, 'E': Displacement.East}


class MonthlyInversion:

    def __init__(self, basis, lln, dataset, components=None, R=None, weights='station', sigma_floor=1e-4):
        """
        :param basis: any pygwater Basis (GreenFunctionBasis, SphericalHarmonicBasis, SlepianBasis)
        :param lln: LoveNumber object (choose the frame consistent with the GNSS series, e.g. CF for seasonal loads)
        :param dataset: GNSSDataset (pygwater.preprocessing)
        :param components: subset of dataset.components, e.g. ('U',) or ('U', 'N', 'E')
        :param R: regularisation matrix (pygwater.inversion.regularization); default identity
        :param weights: 'station' (one sigma per station and component, the median over months; fast) or 'month'
                        (the sigma of each monthly value)
        :param sigma_floor: lower bound of the sigmas [m]
        """
        self.basis, self.lln, self.ds = basis, lln, dataset
        self.components = tuple(components or dataset.components)
        self.R = R
        self.weights = weights
        self.sigma_floor = sigma_floor
        self.A = {c: basis.design_matrix(dataset.lat, dataset.lon, lln, _OBS[c]) for c in self.components}
        self.station_sigma = {c: np.maximum(np.nanmedian(dataset.sigma[c], axis=0), sigma_floor)
                              for c in self.components}
        self._solvers = {}

    # --------------------------------------------------------------------------------------------------------------
    def _month(self, month):
        """rows for one month: A (n_comp * N, n), d, sigma, station mask"""
        i = self.ds.months.index(month)
        ok = np.all([np.isfinite(self.ds.values[c][i]) for c in self.components], axis=0)
        A = np.vstack([self.A[c][ok] for c in self.components])
        d = np.concatenate([self.ds.values[c][i, ok] for c in self.components])
        if self.weights == 'station':
            sg = np.concatenate([self.station_sigma[c][ok] for c in self.components])
        else:
            sg = np.concatenate([np.maximum(self.ds.sigma[c][i, ok], self.sigma_floor) for c in self.components])
        return A, d, sg, ok

    def solver(self, month):
        A, d, sg, ok = self._month(month)
        key = ok.tobytes() if self.weights == 'station' else month
        if key not in self._solvers:
            if len(self._solvers) > 64:
                self._solvers.clear()
            self._solvers[key] = RegularizedLeastSquares(A, sg, self.R)
        return self._solvers[key], d, ok

    def months(self, min_stations=1):
        return [m for m in self.ds.months if self._month(m)[3].sum() >= min_stations]

    def choose_alpha(self, method='gcv', alphas=None, months=None, min_stations=5, **kwargs):
        """
        One regularisation parameter for the whole series: the median of the choices of the individual months
        (more stable than a single month). See RegularizedLeastSquares.choose_alpha for the methods.
        :return: alpha; the monthly choices are kept in self.alpha_by_month
        """
        alphas = np.logspace(-6, 5, 45) if alphas is None else alphas
        months = months or self.months(min_stations)
        self.alpha_by_month = {}
        for m in months:
            sv, d, ok = self.solver(m)
            self.alpha_by_month[m] = sv.choose_alpha(d, alphas=alphas, method=method, **kwargs)[0]
        if not self.alpha_by_month:
            raise ValueError('no month with at least %d stations' % min_stations)
        return float(np.exp(np.median(np.log(list(self.alpha_by_month.values())))))

    def run(self, alpha, months=None, lb=None, ub=None, min_stations=5):
        """
        :param alpha: one value for all months, or 'per_month' to use self.alpha_by_month (run choose_alpha first)
        :param lb, ub: bounds on the parameters (e.g. lb=0 for a non-negative load); bounded solves are slower
        :return: InversionResult
        """
        months = months or self.months(min_stations)
        X, alph, nst, fit, dof = [], [], [], [], []
        for m in months:
            sv, d, ok = self.solver(m)
            a = self.alpha_by_month[m] if alpha == 'per_month' else alpha
            if lb is None and ub is None:
                sol = sv.solve(d, a)
            else:
                sol = sv.solve_bounded(d, a, lb=-np.inf if lb is None else lb, ub=np.inf if ub is None else ub)
            X.append(sol.x), alph.append(a), nst.append(int(ok.sum())), dof.append(sol.dof)
            fit.append(float(np.sqrt(np.mean((sv.A @ sol.x - d) ** 2))))
        return InversionResult(basis=self.basis, months=list(months), x=np.array(X), alpha=np.array(alph),
                               n_stations=np.array(nst), fit_rms=np.array(fit), dof=np.array(dof),
                               reference=self.ds.reference, components=self.components)


@dataclass
class InversionResult:
    basis: object
    months: list
    x: np.ndarray                   # (n_month, n_param)
    alpha: np.ndarray
    n_stations: np.ndarray
    fit_rms: np.ndarray             # RMS of the data misfit [m]
    dof: np.ndarray                 # effective number of parameters
    reference: tuple
    components: tuple
    meta: dict = field(default_factory=dict)

    def ewh(self, lat, lon):
        """EWH anomaly [m] at points, (n_month, N)"""
        return self.basis.load_synthesis(self.x, np.asarray(lat), np.asarray(lon))

    def regional_series(self, lat, lon, inside=None):
        """
        Area-weighted mean EWH [m] over a region, per month, evaluated on a regular grid.
        :param lat, lon: 1d cell-centred grid covering the region (e.g. 0.25 deg)
        :param inside: optional mask (nlat, nlon) or function (lat2d, lon2d) -> bool; default the whole grid
        """
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        w = np.cos(np.deg2rad(la2[m]))
        return self.ewh(la2[m], lo2[m]) @ w / w.sum()

    def regional_mass(self, lat, lon, inside=None):
        """Mass anomaly [Gt] over the region (regional_series x area x water density)"""
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        dlat, dlon = np.deg2rad(abs(lat[1] - lat[0])), np.deg2rad(abs(lon[1] - lon[0]))
        area = np.sum(np.cos(np.deg2rad(la2[m]))) * dlat * dlon * EarthConstant.radiusm ** 2
        return self.regional_series(lat, lon, m) * area * EarthConstant.rhow / 1e12

    def save(self, path):
        np.savez(path, months=np.array(self.months), x=self.x, alpha=self.alpha, n_stations=self.n_stations,
                 fit_rms=self.fit_rms, dof=self.dof, reference=np.array(self.reference),
                 components=np.array(self.components))
