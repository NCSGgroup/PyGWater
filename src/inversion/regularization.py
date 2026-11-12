"""
Regularisation matrices R for the inversion  min ||W (A x - d)||^2 + alpha * x^T R x.

    identity(n)                     : Tikhonov (zeroth order), x^T x; natural for Slepian / SH coefficients, whose
                                      orthonormality makes x^T x proportional to the mean-square EWH
    difference_operator(lat, lon, res) / smoothness(...)
                                    : first differences between neighbouring cells of a regular grid (4-neighbour
                                      graph, longitude wrap-around), R = D^T D = graph Laplacian; penalises roughness
                                      of mascon (point-mass, disk, cluster-disk) solutions
    area_weighted_identity(lat)     : x^T diag(cos lat) x, the mean-square EWH over the cells
"""
import numpy as np


def identity(n):
    """R = I (n x n)"""
    return np.eye(int(n))


def area_weighted_identity(lat):
    """R = diag(cos(lat)), i.e. x^T R x is proportional to the area-weighted mean square of the cell values"""
    return np.diag(np.cos(np.deg2rad(np.asarray(lat, dtype=float))))


def difference_operator(lat, lon, resolution, weighted=True):
    """
    First differences between neighbouring cells (east-west and north-south) of a set of cells on a regular,
    cell-centred grid. Only pairs where both cells are parameters are used, so the domain edge is free.
    :param lat, lon: cell centres [degree], 1d arrays (n,)
    :param resolution: grid spacing [degree]
    :param weighted: scale each difference by sqrt(cos(lat)) so that every edge counts by its length
    :return: D, shape (n_edges, n)
    """
    lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    key = {(int(round(a / resolution * 2)), int(round(((o + 180) % 360) / resolution * 2))): k
           for k, (a, o) in enumerate(zip(lat, lon))}
    step = 2  # keys are in half-resolution units
    nlon = int(round(360 / resolution)) * 2
    rows = []
    for (i, j), k in key.items():
        for di, dj in ((0, step), (step, 0)):  # east and north neighbour
            nb = key.get((i + di, (j + dj) % nlon))
            if nb is not None:
                w = np.sqrt(np.cos(np.deg2rad((lat[k] + lat[nb]) / 2))) if weighted else 1.0
                rows.append((k, nb, w))
    D = np.zeros((len(rows), lat.size))
    for e, (k, nb, w) in enumerate(rows):
        D[e, k], D[e, nb] = -w, w
    return D


def smoothness(lat, lon, resolution, weighted=True):
    """R = D^T D (graph Laplacian of the cell grid)"""
    D = difference_operator(lat, lon, resolution, weighted)
    return D.T @ D


def laplacian_operator(lat, lon, resolution):
    """
    Discrete Laplacian on the cell grid (Harris & Segall, 1987; used by Argus et al., 2014 and Fu et al., 2015):
    (L x)_i = sum over the k_i available 4-neighbours of (x_j - x_i), i.e. a free (Neumann) edge at the domain
    boundary. The penalty is ||L x||^2, i.e. R = L^T L (second-order smoothness), unlike smoothness(), which
    penalises first differences.
    :return: L, shape (n, n)
    """
    D = difference_operator(lat, lon, resolution, weighted=False)
    return -(D.T @ D)  # graph Laplacian with the sign of the continuous operator


def curvature(lat, lon, resolution):
    """R = L^T L with L = laplacian_operator (second-order smoothness, as in Fu et al., 2015)"""
    L = laplacian_operator(lat, lon, resolution)
    return L.T @ L
