"""Tests of src_Inversion (solver, regularisation) and src_OBS (synthetic observations)."""
import numpy as np

from src_auxiliary.GeoMathKit import GeoMathKit
from src_Inversion.Solver import RegularizedLeastSquares
from src_Inversion.Regularization import difference_operator, smoothness, identity
from src_OBS.Synthetic import synthetic_network, checkerboard, add_noise


def test_solver_recovers_overdetermined_problem():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(300, 40)); x = rng.normal(size=40)
    d = A @ x + rng.normal(size=300) * 0.01
    s = RegularizedLeastSquares(A, 0.01)
    alpha, curve = s.choose_alpha(d)
    sol = s.solve(d, alpha)
    assert np.abs(sol.x - x).max() < 0.01 and abs(sol.sigma0 - 1) < 0.2 and sol.dof > 39
    # several epochs in one pass
    X = s.solve(np.vstack([d, 2 * d]), alpha).x
    np.testing.assert_allclose(X[1], 2 * X[0], rtol=1e-10)
    # resolution matrix -> identity, covariance -> (A^T W A)^-1 for small alpha
    np.testing.assert_allclose(s.resolution_matrix(1e-10), np.eye(40), atol=1e-6)
    np.testing.assert_allclose(s.covariance(1e-10), np.linalg.inv(s.N), rtol=1e-4)
    assert s.choose_alpha(d, method='lcurve')[0] > 0 and s.choose_alpha(d, method='wgcv', omega=1.2)[0] > 0


def test_difference_operator_wraps_and_annihilates_constants():
    lat = np.array([1., 1., 3., 3.]); lon = np.array([179., -179., 179., -179.])
    D = difference_operator(lat, lon, 2.0)
    assert D.shape == (4, 4)
    np.testing.assert_allclose(D @ np.ones(4), 0, atol=1e-12)
    R = smoothness(lat, lon, 2.0)
    np.testing.assert_allclose(R, R.T) and np.testing.assert_allclose(R @ np.ones(4), 0, atol=1e-12)
    assert identity(3).shape == (3, 3)


def test_synthetic_network_and_checkerboard():
    lat, lon = synthetic_network(lambda a, o: (a > 20) & (a < 40) & (o > 100) & (o < 120), n_target=40,
                                 min_dist=2.0, bbox=(15, 45, 95, 125), seed=3)
    assert len(lat) == 40 and ((lat > 20) & (lat < 40)).all()
    d = np.array([GeoMathKit.angular_distance(lat[i], lon[i], np.delete(lat, i), np.delete(lon, i)).min()
                  for i in range(len(lat))])
    assert d.min() >= 2.0
    c = checkerboard(np.array([1., 1., 5.]), np.array([1., 5., 1.]), 4.0, 0.1)
    assert c[0] == -c[1] and c[0] == -c[2] and abs(c[0]) == 0.1
    assert np.std(add_noise(np.zeros(20000), 2.0) - 0) > 1.9
