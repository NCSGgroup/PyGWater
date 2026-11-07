"""
Synthetic GNSS observations for closed-loop and checkerboard tests.

    synthetic_network(inside, n_target, min_dist, seed)   quasi-uniform station positions inside a region
    add_noise(d, sigma, seed)                              white Gaussian noise
    checkerboard(lat, lon, block, amplitude)               +/- amplitude in blocks of `block` degrees
"""
import numpy as np

from src_auxiliary.GeoMathKit import GeoMathKit


def synthetic_network(inside, n_target=260, min_dist=1.5, bbox=(-90, 90, -180, 180), seed=1, max_tries=200000):
    """
    Poisson-disk (dart-throwing) sampling: random points inside the region, at least `min_dist` degrees apart.
    :param inside: function (lat, lon) -> bool array, True inside the region
    :param n_target: stop after this many stations
    :param min_dist: minimum angular distance between stations [degree]
    :param bbox: (lat_min, lat_max, lon_min, lon_max) to draw from
    :return: lat, lon of the stations [degree]
    """
    rng = np.random.default_rng(seed)
    la0, la1, lo0, lo1 = bbox
    out_lat, out_lon = [], []
    tries = 0
    while len(out_lat) < n_target and tries < max_tries:
        n = 2000
        z = rng.uniform(np.sin(np.deg2rad(la0)), np.sin(np.deg2rad(la1)), n)  # uniform on the sphere
        cand_lat, cand_lon = np.rad2deg(np.arcsin(z)), rng.uniform(lo0, lo1, n)
        ok = inside(cand_lat, cand_lon)
        for a, o in zip(cand_lat[ok], cand_lon[ok]):
            tries += 1
            if out_lat and GeoMathKit.angular_distance(a, o, np.array(out_lat), np.array(out_lon)).min() < min_dist:
                continue
            out_lat.append(a), out_lon.append(o)
            if len(out_lat) >= n_target:
                break
        tries += n - ok.sum()
    return np.array(out_lat), np.array(out_lon)


def add_noise(d, sigma, seed=2):
    """d + N(0, sigma^2), sigma scalar or array broadcastable to d"""
    rng = np.random.default_rng(seed)
    return d + rng.standard_normal(np.shape(d)) * sigma


def checkerboard(lat, lon, block, amplitude=0.1, offset=(0.0, 0.0)):
    """
    +/- amplitude in square blocks of `block` degrees (EWH [m] if amplitude in metres).
    :param lat, lon: cell centres [degree] (any shape)
    """
    i = np.floor((np.asarray(lat) - offset[0] + 90) / block).astype(int)
    j = np.floor((np.asarray(lon) - offset[1] + 180) / block).astype(int)
    return amplitude * np.where((i + j) % 2 == 0, 1.0, -1.0)
