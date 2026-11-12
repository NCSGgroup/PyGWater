"""
GRACE/GRACE-FO mascon solutions for comparison with GNSS-inverted loads.

    from pygwater.observations.grace import Mascon, compare
    g = Mascon.read(path)                                   # CSR RL06.x mascons (NetCDF), JPL in the same layout
    g = g.anomalies(reference=('2004-01', '2009-12'))       # relative to the same period as the GNSS data
    s = g.regional_series(lat, lon, inside)                 # area-weighted mean EWH over a region [m]
    stats = compare(months_gnss, s_gnss, g.months, s)       # correlation, RMS, NSE, amplitude ratio

EWH is in metres throughout (the CSR files store centimetres in `lwe_thickness`).
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Mascon:
    months: list                    # 'YYYY-MM'
    lat: np.ndarray                 # (nlat,) increasing
    lon: np.ndarray                 # (nlon,) -180..180, increasing
    ewh: np.ndarray                 # (n_month, nlat, nlon) [m]

    @classmethod
    def read(cls, path, variable='lwe_thickness', units_to_m=None):
        """
        Read a mascon NetCDF (time in days since an epoch, lat, lon, lwe_thickness). Several solutions in one
        calendar month are averaged.
        :param units_to_m: scale to metres; default from the variable's units attribute (cm -> 0.01)
        """
        import netCDF4 as nc
        with nc.Dataset(path) as f:
            t = f['time']
            dates = nc.num2date(t[:], t.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            lat = np.asarray(f['lat'][:], dtype=float)
            lon = np.asarray(f['lon'][:], dtype=float)
            v = f[variable]
            if units_to_m is None:
                u = getattr(v, 'units', 'cm').lower()
                units_to_m = 0.01 if u.startswith('c') else (0.001 if u.startswith('mm') else 1.0)
            data = np.asarray(v[:], dtype=float) * units_to_m
        lon = (lon + 180) % 360 - 180
        order = np.argsort(lon)
        lon, data = lon[order], data[:, :, order]
        if lat[0] > lat[-1]:
            lat, data = lat[::-1], data[:, ::-1]
        keys = ['%04d-%02d' % (d.year, d.month) for d in dates]
        months = sorted(set(keys))
        ewh = np.array([data[[i for i, k in enumerate(keys) if k == m]].mean(axis=0) for m in months])
        return cls(months, lat, lon, ewh)

    def anomalies(self, reference):
        """relative to the mean over the months of the reference period ('YYYY-MM', 'YYYY-MM'), inclusive"""
        sel = [i for i, m in enumerate(self.months) if reference[0] <= m <= reference[1]]
        if not sel:
            raise ValueError('no GRACE month in the reference period %s' % (reference,))
        return Mascon(self.months, self.lat, self.lon, self.ewh - self.ewh[sel].mean(axis=0))

    def sample(self, lat, lon):
        """EWH at points (nearest grid node), (n_month, N)"""
        i = np.clip(np.searchsorted(self.lat, lat), 1, self.lat.size - 1)
        i -= (np.abs(self.lat[i - 1] - lat) < np.abs(self.lat[i] - lat))
        lon = (np.asarray(lon) + 180) % 360 - 180
        j = np.clip(np.searchsorted(self.lon, lon), 1, self.lon.size - 1)
        j -= (np.abs(self.lon[j - 1] - lon) < np.abs(self.lon[j] - lon))
        return self.ewh[:, i, j]

    def regional_series(self, lat, lon, inside=None):
        """area-weighted mean EWH [m] over a region on the grid (lat, lon), as InversionResult.regional_series"""
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        w = np.cos(np.deg2rad(la2[m]))
        return self.sample(la2[m], lo2[m]) @ w / w.sum()


def compare(months_a, a, months_b, b):
    """
    Agreement of two monthly series on their common months (a = estimate, b = reference, e.g. GNSS vs GRACE).
    :return: dict n, r (correlation), rms (of a - b), nse (Nash-Sutcliffe efficiency of a against b),
             amplitude_ratio (std a / std b), bias (mean a - b), months
    """
    ib = {m: i for i, m in enumerate(months_b)}
    common = [m for m in months_a if m in ib]
    ia = {m: i for i, m in enumerate(months_a)}
    x = np.array([a[ia[m]] for m in common], dtype=float)
    y = np.array([b[ib[m]] for m in common], dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return dict(n=int(x.size), r=np.nan, rms=np.nan, nse=np.nan, amplitude_ratio=np.nan, bias=np.nan, months=[])
    return dict(n=int(x.size), r=float(np.corrcoef(x, y)[0, 1]), rms=float(np.sqrt(np.mean((x - y) ** 2))),
                nse=float(1 - np.sum((x - y) ** 2) / np.sum((y - y.mean()) ** 2)),
                amplitude_ratio=float(np.std(x) / np.std(y)), bias=float(np.mean(x - y)),
                months=[m for m, k in zip(common, ok) if k])
