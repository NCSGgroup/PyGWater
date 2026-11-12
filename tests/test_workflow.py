"""
End to end: synthetic daily GNSS files from a known seasonal load -> preprocessing -> monthly inversion -> comparison
with a synthetic GRACE mascon file holding the same load.
"""
import numpy as np
import pytest

import pygwater as pgw
from pygwater.preprocessing import GNSSNetwork
from pygwater.inversion.workflow import MonthlyInversion
from pygwater.observations.grace import Mascon, compare

RES = 0.5
LAT1 = np.arange(43 + RES / 2, 49, RES)
LON1 = np.arange(3 + RES / 2, 10, RES)
MON = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def pattern(lat, lon):
    """spatial pattern of the load [m EWH per unit of the time function]"""
    return 0.15 * np.exp(-((lat - 46) ** 2 + (lon - 6.5) ** 2) / 2.0)


def time_function(t):
    return np.sin(2 * np.pi * t) + 0.3 * np.sin(2 * np.pi * t / 4)       # annual + a 4-year wave


def write_network(folder, cells, lln, n=40, seed=0):
    rng = np.random.default_rng(seed)
    slat, slon = rng.uniform(44, 48, n), rng.uniform(4, 9, n)
    gain = cells.design_matrix(slat, slon, lln, pgw.Displacement.Vertical) @ pattern(cells.lat, cells.lon)
    days = np.arange(np.datetime64('2010-01-01'), np.datetime64('2018-01-01'))
    t = np.array([int(str(d)[:4]) + ((d - np.datetime64(str(d)[:4] + '-01-01')).astype(int) + 0.5) /
                  (366 if int(str(d)[:4]) % 4 == 0 else 365) for d in days])
    head = 'site YYMMMDD yyyy.yyyy __MJD week d reflon _e0(m) __east(m) ____n0(m) _north(m) u0(m) ____up(m) ...\n'
    files = []
    for k in range(n):
        up = gain[k] * time_function(t) + rng.normal(0, 0.003) * (t - 2014) + rng.normal(0, 2e-3, t.size)
        if k == 3:
            up += np.where(t > 2013.7, 0.015, 0.0)                         # an unreported offset
        lines = [head] + ['S%03d %02d%s%02d %.4f 0 0 0 0 0 0.0 0 0.0 0 %.5f 0 0.001 0.001 0.003 0 0 0 %.4f %.4f 0\n'
                          % (k, int(str(d)[2:4]), MON[int(str(d)[5:7]) - 1], int(str(d)[8:10]), tt, u,
                             slat[k], slon[k]) for d, tt, u in zip(days, t, up)]
        f = folder / ('S%03d.tenv3' % k)
        f.write_text(''.join(lines))
        files.append(f)
    return files


def write_grace(path, months):
    import netCDF4 as nc
    lat = np.arange(-89.875, 90, 0.25)
    lon = np.arange(0.125, 360, 0.25)
    lo2, la2 = np.meshgrid((lon + 180) % 360 - 180, lat)
    with nc.Dataset(path, 'w') as f:
        f.createDimension('time', len(months)), f.createDimension('lat', lat.size), f.createDimension('lon', lon.size)
        t = f.createVariable('time', 'f8', ('time',)); t.units = 'days since 2002-01-01 00:00:00'
        f.createVariable('lat', 'f4', ('lat',))[:] = lat
        f.createVariable('lon', 'f4', ('lon',))[:] = lon
        v = f.createVariable('lwe_thickness', 'f4', ('time', 'lat', 'lon')); v.units = 'cm'
        P = pattern(la2, lo2)
        for i, m in enumerate(months):
            tm = int(m[:4]) + (int(m[5:7]) - 0.5) / 12
            t[i] = (np.datetime64(m + '-15') - np.datetime64('2002-01-01')).astype(int)
            v[i] = 100 * P * time_function(tm)


@pytest.fixture(scope='module')
def lln():
    return pgw.LoveNumber().config(lmax=2000, method=pgw.LLN_Data.PREM).get_Love_number()


def test_gnss_to_grace(tmp_path, lln):
    lo2, la2 = np.meshgrid(LON1, LAT1)
    cells = pgw.GreenFunctionBasis.from_grid(LAT1, LON1, RES, kernel='point_mass')
    files = write_network(tmp_path, cells, lln)

    net = GNSSNetwork.from_files(files).clean(detect_jumps=True).monthly().spatial_filter(radius_km=300)
    assert net.stations[3].qc['n_detected_jumps'] == 1 and len(net.stations) == 40
    ds = net.to_dataset(components=('U',), reference=('2011-01', '2016-12'))

    R = pgw.regularization.curvature(cells.lat, cells.lon, RES)
    inv = MonthlyInversion(cells, lln, ds, components=('U',), R=R)
    alpha = inv.choose_alpha(method='gcv')
    res = inv.run(alpha)
    assert len(res.months) == 96 and (res.n_stations == 40).all()

    # regional mean over the station area against the truth
    glat, glon = np.arange(44.125, 48, 0.25), np.arange(4.125, 9, 0.25)
    gl2, gb2 = np.meshgrid(glon, glat)
    w = np.cos(np.deg2rad(gb2))
    tm = np.array([int(m[:4]) + (int(m[5:7]) - 0.5) / 12 for m in res.months])
    truth = np.sum(w * pattern(gb2, gl2)) / w.sum() * time_function(tm)
    truth = truth - truth[[i for i, m in enumerate(res.months) if '2011-01' <= m <= '2016-12']].mean()
    est = res.regional_series(glat, glon)
    assert np.corrcoef(est, truth)[0, 1] > 0.97

    # the same comparison through a GRACE mascon file
    write_grace(tmp_path / 'mascon.nc', res.months)
    g = Mascon.read(tmp_path / 'mascon.nc').anomalies(ds.reference)
    stats = compare(res.months, est, g.months, g.regional_series(glat, glon))
    assert stats['n'] == 96 and stats['r'] > 0.97 and 0.6 < stats['amplitude_ratio'] < 1.4
    mass = res.regional_mass(glat, glon)
    assert mass.shape == (96,) and np.isfinite(mass).all()
