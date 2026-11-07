"""
Regularised weighted least squares for d = A x + e,  e ~ N(0, diag(sigma^2)):

        x_alpha = argmin  || W (A x - d) ||^2 + alpha  x^T R x ,      W = diag(1 / sigma)
                = (N + alpha R)^-1  A^T W^2 d ,                      N = A^T W^2 A

Parameter choice (choose_alpha): generalised cross-validation (GCV, Golub et al., 1979), weighted GCV (WGCV, Chung
et al., 2008; omega > 1 guards against the spurious minimum at tiny alpha that appears when the data contain
modelling errors) or the corner of the L-curve (Hansen, 1992). A plausibility test (`plausible`, e.g. a physical
bound on the recovered EWH) removes candidate alphas whose solutions are not physical before the choice; it guards
against the GCV failure seen when the inversion basis does not match the data-generating physics. Quality
measures for a chosen alpha:

    resolution matrix   Rm = (N + alpha R)^-1 N                       (x_alpha = Rm x_true + noise term)
    noise covariance    Cx = (N + alpha R)^-1 N (N + alpha R)^-1      (formal, propagated from sigma)
    effective dof       trace(Rm)

All solves use a Cholesky factorisation of the (small, n_param x n_param) regularised normal matrix, so many
epochs (columns of d) are inverted in one pass.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve


@dataclass
class Solution:
    x: np.ndarray            # (n,) or (ntime, n)
    alpha: float
    residual_norm: float     # || W (A x - d) ||  (first epoch if several)
    seminorm: float          # sqrt(x^T R x)
    dof: float               # trace of the resolution matrix
    sigma0: float            # a-posteriori unit weight


class RegularizedLeastSquares:

    def __init__(self, A, sigma, R=None):
        """
        :param A: design matrix (N_obs x n), e.g. ForwardModel.design_matrix(...) or basis.design_matrix(...)
        :param sigma: observation standard deviations (N_obs,) or scalar, same units as d
        :param R: regularisation matrix (n x n), see Regularization.py; default identity
        """
        self.A = np.asarray(A, dtype=float)
        self.w = 1.0 / np.broadcast_to(np.asarray(sigma, dtype=float), (self.A.shape[0],))
        self.Aw = self.A * self.w[:, None]
        self.N = self.Aw.T @ self.Aw
        self.R = np.eye(self.A.shape[1]) if R is None else np.asarray(R, dtype=float)
        # scale alpha relative to the data: alpha = 1 balances the traces of N and R
        self.scale = np.trace(self.N) / max(np.trace(self.R), 1e-300)

    def _factor(self, alpha):
        return cho_factor(self.N + alpha * self.scale * self.R + 1e-12 * self.scale * np.eye(self.N.shape[0]))

    def solve(self, d, alpha):
        """
        :param d: observations (N_obs,) or (ntime, N_obs)
        :param alpha: regularisation parameter (relative: the matrix used is alpha * trace(N)/trace(R) * R)
        :return: Solution
        """
        d = np.asarray(d, dtype=float)
        F = self._factor(alpha)
        b = (d * self.w) @ self.Aw  # (n,) or (ntime, n)
        x = cho_solve(F, b.T).T
        r = (x @ self.A.T - d) * self.w
        r0, x0 = (r if r.ndim == 1 else r[0]), (x if x.ndim == 1 else x[0])
        dof = float(np.trace(cho_solve(F, self.N)))
        sigma0 = float(np.sqrt(np.sum(r0 ** 2) / max(len(r0) - dof, 1)))
        return Solution(x=x, alpha=alpha, residual_norm=float(np.linalg.norm(r0)),
                        seminorm=float(np.sqrt(max(x0 @ self.R @ x0, 0))), dof=dof, sigma0=sigma0)

    def choose_alpha(self, d, alphas=None, method='gcv', omega=1.0, plausible=None):
        """
        :param d: observations (N_obs,); for a time series pass one representative epoch or the stacked mean
        :param alphas: candidates (relative scale); default 10^-6 ... 10^2
        :param method: 'gcv' (minimum of the GCV function), 'wgcv' (weighted GCV with weight omega, e.g. 1.2) or
                       'lcurve' (maximum curvature of the L-curve)
        :param omega: WGCV weight; GCV(alpha) = m ||r||^2 / (m - omega * dof)^2  ('gcv' uses omega = 1)
        :param plausible: optional function x -> bool; alphas whose solution fails it are not eligible (e.g.
                          lambda x: np.sqrt(np.mean(basis.load_synthesis(x, lat, lon) ** 2)) < 1.0  for EWH < 1 m)
        :return: (best alpha, dict with 'alpha', 'residual', 'seminorm', 'gcv')
        """
        alphas = np.logspace(-6, 2, 33) if alphas is None else np.asarray(alphas, dtype=float)
        d = np.asarray(d, dtype=float)
        m = len(d)
        res, semi, gcv, ok = [], [], [], []
        for a in alphas:
            s = self.solve(d, a)
            ok.append(True if plausible is None else bool(plausible(s.x)))
            res.append(s.residual_norm), semi.append(s.seminorm)
            om = omega if method == 'wgcv' else 1.0
            gcv.append(m * s.residual_norm ** 2 / max(m - om * s.dof, 1e-9) ** 2 if m - om * s.dof > 0 else np.inf)
        res, semi, gcv, ok = np.array(res), np.array(semi), np.array(gcv), np.array(ok)
        if not ok.any():
            raise ValueError('no alpha gives a plausible solution; widen the alpha range or relax the test')
        if method in ('gcv', 'wgcv'):
            best = alphas[int(np.argmin(np.where(ok, gcv, np.inf)))]
        elif method == 'lcurve':
            lr, ls = np.log(res), np.log(np.maximum(semi, 1e-300))
            t = np.log(alphas)
            d1r, d1s = np.gradient(lr, t), np.gradient(ls, t)
            d2r, d2s = np.gradient(d1r, t), np.gradient(d1s, t)
            kappa = (d1r * d2s - d2r * d1s) / np.maximum((d1r ** 2 + d1s ** 2) ** 1.5, 1e-300)
            kappa = np.where(ok, kappa, -np.inf)
            best = alphas[int(np.argmax(kappa[1:-1])) + 1]
        else:
            raise ValueError(method)
        return float(best), {'alpha': alphas, 'residual': res, 'seminorm': semi, 'gcv': gcv, 'plausible': ok}

    def resolution_matrix(self, alpha):
        """Rm = (N + alpha R)^-1 N; row k shows how the true parameters blend into estimate k"""
        return cho_solve(self._factor(alpha), self.N)

    def covariance(self, alpha):
        """formal noise covariance of x_alpha, (N + alpha R)^-1 N (N + alpha R)^-1"""
        F = self._factor(alpha)
        M = cho_solve(F, self.N)
        return cho_solve(F, M.T).T
