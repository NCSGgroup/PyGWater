"""
BASIS built from the kernels of point_load.py and disk_load.py: many cells, any stations, design matrix.
One class, three methods of the paper selected with `kernel`:

    finite LGF : GreenFunctionBasis.from_grid(lat, lon, res, mask, kernel='point_mass')
    D-LGF      : GreenFunctionBasis.from_grid(lat, lon, res, mask, kernel='disk')
    CD-LGF     : GreenFunctionBasis.from_grid(lat, lon, res, mask, kernel='cluster_disk', subgrid=0.1)

Green's functions as a Basis (src/basis/base.py): the load is a set of cells (mascons), each carrying a uniform EWH
x_j [m], and the response is the convolution with the load Green's function (Sect. 3-6 of Yang et al., 2025):

    kernel='point_mass' : finite (truncated) point-mass LGF, Eq. 8;  cell mass = rho_w * area_j * x_j
                          (load assumption I/II: the cell mass is concentrated at its centre)
    kernel='disk'       : disk-load LGF (D-LGF), Eq. 17; a disk of radius r_j and thickness x_j per cell
                          (load assumption III: uniform load over the disk)
    kernel='cluster_disk': clustered disk-load LGF (CD-LGF, Sect. 7): each cell (size = resolution) is split into
                          sub-cells of size `subgrid` (e.g. 0.1 deg), each an equal-area disk of thickness x_j; the
                          column of cell j is the sum of its sub-disk responses (load assumption III, recommended
                          for mascons). Kernels are evaluated from a DiskLoadTable (one Legendre series per distinct
                          sub-disk radius, interpolated in distance), so any points can be used.

        d(points) = A x,     A_ij = K(psi_ij; observable) * s_j,

with psi_ij the angular distance between point i and cell j, K the kernel of the observable and s_j = rho_w area_j
(point mass) or 1 (disk, cluster_disk: the kernels are already per metre thickness).

The Love numbers must reach a degree consistent with the cell size: L >= 180 / resolution for the point mass
(Sect. 4) and N >= max(720 / r, 10000) for the disk and the sub-disks of the cluster (Eq. 24, see
DiskLoad.required_lmax; about 12 800 for a 0.1-deg sub-grid at the equator).
"""

import numpy as np

from pygwater.utils.constants import EarthConstant
from pygwater.utils.enums import Displacement
from pygwater.utils.geomath import GeoMathKit
from pygwater.basis.base import Basis
from pygwater.basis.green.point_load import LGF_truncation
from pygwater.basis.green.disk_load import DiskLoad, DiskLoadTable


class GreenFunctionBasis(Basis):
    name = 'green_function'

    def __init__(self, lat, lon, kernel='point_mass', area=None, radius=None, resolution=None, subgrid=None,
                 check_lmax=True, tabulate=True):
        """
        :param lat, lon: cell centres [degree], 1d arrays (n_cell,)
        :param kernel: 'point_mass', 'disk' or 'cluster_disk'
        :param area: cell areas [m^2] (point_mass)
        :param radius: disk radii [degree] (disk)
        :param resolution: grid spacing [degree] of a regular grid (optional for point_mass/disk, where it is used
                           by load_synthesis to leave points outside every cell at zero; required for cluster_disk,
                           where it is the cell size that is partitioned)
        :param subgrid: sub-cell size [degree] of cluster_disk (e.g. 0.1; 0.1-0.2 recommended by Yang et al., 2025);
                        resolution / subgrid must be an integer
        :param check_lmax: warn if the Love numbers are too short for the disk radius (Eq. 24)
        :param tabulate: disk kernel only: True (default) evaluates the D-LGF from a DiskLoadTable (one Legendre
                         series per distinct radius, interpolated, ~1e-5 of the peak; fast for many cells); False
                         evaluates the series of DiskLoad directly at every point (exact, slow for large N)
        """
        self.lat = np.atleast_1d(np.asarray(lat, dtype=float))
        self.lon = np.atleast_1d(np.asarray(lon, dtype=float))
        assert self.lat.shape == self.lon.shape and self.lat.ndim == 1
        if kernel not in ('point_mass', 'disk', 'cluster_disk'):
            raise ValueError(kernel)
        self.kernel = kernel
        self.area, self.radius, self.subgrid = None, None, None
        if kernel == 'cluster_disk':
            if resolution is None or subgrid is None:
                raise ValueError('kernel="cluster_disk" needs the cell size (resolution) and the sub-grid size')
            ns = resolution / subgrid
            if abs(ns - round(ns)) > 1e-9:
                raise ValueError('resolution / subgrid must be an integer')
            self.subgrid = float(subgrid)
            self._ns = int(round(ns))
        elif kernel == 'point_mass':
            if area is None:
                raise ValueError('kernel="point_mass" needs the cell areas [m^2]')
            self.area = np.broadcast_to(np.asarray(area, dtype=float), self.lat.shape).copy()
        else:
            if radius is None:
                raise ValueError('kernel="disk" needs the disk radii [degree]')
            self.radius = np.broadcast_to(np.asarray(radius, dtype=float), self.lat.shape).copy()
        self.resolution = resolution
        self.check_lmax = check_lmax
        self.tabulate = tabulate
        self._cache = {}  # one entry: (Love-number object, its configuration) -> kernel objects
        self._tree = None

    @classmethod
    def from_grid(cls, lat, lon, resolution, mask=None, kernel='point_mass', **kwargs):
        """
        Cells of a regular cell-centred grid (1d lat, lon), optionally restricted to mask (nlat, nlon) > 0.
        Point mass: exact cell area a^2 cos(lat) dlat dlon;  disk: equal-area radius (RegularGrid.grid2radius);
        cluster_disk: pass subgrid=... (e.g. 0.1).
        """
        from pygwater.mesh.regular_grid import regular
        lo2, la2 = np.meshgrid(np.asarray(lon, dtype=float), np.asarray(lat, dtype=float))
        sel = np.ones(la2.shape, bool) if mask is None else np.asarray(mask) > 0
        la, lo = la2[sel], lo2[sel]
        if kernel == 'point_mass':
            area = np.cos(np.deg2rad(la)) * np.deg2rad(resolution) ** 2 * EarthConstant.radiusm ** 2
            return cls(la, lo, kernel, area=area, resolution=resolution, **kwargs)
        if kernel == 'cluster_disk':
            return cls(la, lo, kernel, resolution=resolution, **kwargs)
        radius = np.round(regular.grid2radius(lat_center=la, grid_size=resolution)[0], 12)
        return cls(la, lo, kernel, radius=radius, resolution=resolution, **kwargs)

    @property
    def n_param(self):
        return self.lat.size

    # ------------------------------------------------------------------------------------------------------------
    def _sub_offsets(self):
        """offsets of the sub-cell centres from the cell centre [degree], (ns,)"""
        return (np.arange(self._ns) + 0.5) * self.subgrid - self.resolution / 2

    def sub_cells(self, j):
        """sub-cell centres (lat, lon) [degree] and equal-area disk radii [degree] of cell j (cluster_disk)"""
        from pygwater.mesh.regular_grid import regular
        off = self._sub_offsets()
        lo2, la2 = np.meshgrid(self.lon[j] + off, self.lat[j] + off)
        r = np.round(regular.grid2radius(lat_center=la2, grid_size=self.subgrid)[0], 10)
        return la2, lo2, r

    def _kernels(self, lln):
        """kernel objects for this Earth model; one per distinct disk radius, or one table for all sub-disks"""
        # keyed on the object itself (kept alive by the cache, so its id cannot be reused by another object) and on
        # its configuration (config() on the same object changes the Love numbers)
        key = (id(lln), lln.method, lln.lmax, lln.frame)
        cached = self._cache.get(key)
        if cached is not None and cached[0] is lln:
            return cached[1]
        if self.kernel == 'point_mass':
            kernels = LGF_truncation(lln=lln)
        elif self.kernel == 'disk' and self.tabulate:
            kernels = DiskLoadTable(lln, np.unique(self.radius), thickness=1, check_lmax=self.check_lmax)
        elif self.kernel == 'disk':
            kernels = {r: DiskLoad(lln=lln, radius=r, thickness=1, check_lmax=self.check_lmax)
                       for r in np.unique(self.radius)}
        else:
            from pygwater.mesh.regular_grid import regular
            sub_lat = (self.lat[:, None] + self._sub_offsets()[None, :]).ravel()
            radii = regular.grid2radius(lat_center=np.unique(sub_lat), grid_size=self.subgrid)[0]
            kernels = DiskLoadTable(lln, radii, thickness=1, check_lmax=self.check_lmax)
        self._cache = {key: (lln, kernels)}
        return kernels

    @staticmethod
    def _table_response(T, dlat, dlon, radius, lat, lon, variable):
        """response at the points (N,) of disks of one radius centred at (dlat, dlon) (M,), from a table: (M, N)"""
        kind = {Displacement.Vertical: 'vertical', Displacement.Geoheight: 'geoid', Displacement.Horizontal:
                'horizontal', Displacement.North: 'horizontal', Displacement.East: 'horizontal'}.get(variable)
        if kind is None:
            raise NotImplementedError('disk kernels support Vertical, North, East, Horizontal and Geoheight (got %s); '
                                      'surface gravity is not available for disk kernels' % variable)
        dla, dlo = np.atleast_1d(dlat)[:, None], np.atleast_1d(dlon)[:, None]
        psi = GeoMathKit.angular_distance(dla, dlo, lat[None, :], lon[None, :])
        K = T.evaluate(kind, radius, psi)
        if variable in (Displacement.North, Displacement.East):
            un, ue = GeoMathKit.horizontal_to_north_east(K, dla, dlo, lat[None, :], lon[None, :])
            K = un if variable == Displacement.North else ue
        return K

    def _cluster_column(self, T, j, lat, lon, variable):
        """response of cell j (sum of its sub-disks) at the points"""
        if variable == Displacement.Horizontal:
            raise NotImplementedError('the great-circle Horizontal component is not defined for a cluster; '
                                      'use North and East')
        la2, lo2, r2 = self.sub_cells(j)
        col = np.zeros(lat.size)
        for row in range(la2.shape[0]):  # one sub-cell row shares one radius
            col += self._table_response(T, la2[row], lo2[row], r2[row, 0], lat, lon, variable).sum(axis=0)
        return col

    def design_matrix(self, lat, lon, lln, variable=Displacement.Vertical, **kwargs):
        """
        A (N x n_cell): observable at the points per metre EWH in each cell.
        :param variable: Vertical, North, East, Geoheight (all kernels); Horizontal (point_mass, disk);
                         Gravity (point_mass only; kwarg newtonian)
        """
        lat, lon = np.atleast_1d(lat).astype(float), np.atleast_1d(lon).astype(float)
        K = self._kernels(lln)
        A = np.empty((lat.size, self.n_param))
        for j in range(self.n_param):
            if self.kernel == 'cluster_disk':
                A[:, j] = self._cluster_column(K, j, lat, lon, variable)
                continue
            if self.kernel == 'disk' and self.tabulate:
                A[:, j] = self._table_response(K, self.lat[j], self.lon[j], self.radius[j], lat, lon, variable)[0]
                continue
            if self.kernel == 'point_mass':
                obj, scale = K, EarthConstant.rhow * self.area[j]
            else:
                obj, scale = K[self.radius[j]], 1.0
            A[:, j] = scale * GeoMathKit.kernel_at_points(obj, self.lat[j], self.lon[j], lat, lon, variable,
                                                          **kwargs)
        return A

    def _xyz(self, lat, lon):
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)

    def _cell_index(self, lat, lon):
        """index of the cell containing each point (nearest centre), -1 if outside every cell (needs resolution)"""
        from scipy.spatial import cKDTree
        if self._tree is None:
            self._tree = cKDTree(self._xyz(self.lat, self.lon))
        lat, lon = np.atleast_1d(lat).astype(float), np.atleast_1d(lon).astype(float)
        _, idx = self._tree.query(self._xyz(lat, lon))
        if self.resolution is not None:
            h = self.resolution / 2 * (1 + 1e-9)
            dlon = (lon - self.lon[idx] + 180) % 360 - 180
            idx = np.where((np.abs(lat - self.lat[idx]) <= h) & (np.abs(dlon) <= h), idx, -1)
        return idx

    def load_synthesis(self, x, lat, lon):
        """EWH [m] at the points: value of the cell containing the point (piecewise constant), 0 outside."""
        x = np.asarray(x, dtype=float)
        idx = self._cell_index(lat, lon)
        out = np.take(x, np.maximum(idx, 0), axis=-1)
        return np.where(idx >= 0, out, 0.0)

    def load_parameters(self, ewh, lat, lon):
        """x = EWH of the gridded load at the cell centres (nearest grid node), (n_cell,) or (ntime, n_cell)."""
        ewh = np.asarray(ewh, dtype=float)
        lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
        i = np.abs(lat[:, None] - self.lat[None, :]).argmin(axis=0)
        dlon = (lon[:, None] - self.lon[None, :] + 180) % 360 - 180
        j = np.abs(dlon).argmin(axis=0)
        return ewh[..., i, j]
