"""Tests of pygwater.preprocessing on a synthetic NGL .tenv3 file (no external data needed)."""
import numpy as np

from pygwater.preprocessing.network import GNSSNetwork
from pygwater.preprocessing.trajectory import fit_trajectory
from pygwater.preprocessing.series import (read_tenv3, read_ngl_steps, _decyear_to_month,
                          _yymmmdd_to_decyear)

HEADER = ('site YYMMMDD yyyy.yyyy __MJD week d reflon _e0(m) __east(m) ____n0(m) _north(m) u0(m) ____up(m) '
          '_ant(m) sig_e(m) sig_n(m) sig_u(m) __corr_en __corr_eu __corr_nu _latitude(deg) _longitude(deg) '
          '__height(m)\n')
MON = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def write_tenv3(path, site='TEST', rate=0.003, annual=0.004, step_at=2018.5, step=0.01, seed=0):
    rng = np.random.default_rng(seed)
    days = np.arange(np.datetime64('2016-01-01'), np.datetime64('2021-01-01'))
    lines = [HEADER]
    truth = []
    for k, d in enumerate(days):
        y = int(str(d)[:4]); m = int(str(d)[5:7]); dd = int(str(d)[8:10])
        ndays = 366 if y % 4 == 0 else 365
        t = y + ((d - np.datetime64('%d-01-01' % y)).astype(int) + 0.5) / ndays
        up = rate * (t - 2018.5) + annual * np.sin(2 * np.pi * (t - 2018.5)) + (step if t >= step_at else 0)
        truth.append(annual * np.sin(2 * np.pi * (t - 2018.5)))
        noise = rng.normal(scale=0.002)
        if k == 100:
            noise = 0.2  # an outlier
        lines.append('%s %02d%s%02d %.4f 57000 2000 0 -1.0 123456 0.1 3456789 0.2 10 %.5f 0.0 0.001 0.001 0.003 '
                     '0 0 0 30.5 101.25 1500.0\n' % (site, y % 100, MON[m - 1], dd, t, up + noise))
    path.write_text(''.join(lines))
    return np.array(truth)


def test_read_clean_monthly_epoch(tmp_path):
    f = tmp_path / 'TEST.tenv3'
    write_tenv3(f)
    (tmp_path / 'steps.txt').write_text('TEST  18JUL02  1  Equipment change\nOTHR  19JAN01  2  quake\n')
    s = read_tenv3(f)
    assert s.name == 'TEST' and abs(s.lat - 30.5) < 1e-9 and abs(s.lon - 101.25) < 1e-9 and len(s.t) == 1827
    steps = read_ngl_steps(tmp_path / 'steps.txt')
    assert set(steps) == {'TEST', 'OTHR'} and abs(steps['TEST'][0] - 2018.5) < 0.01
    p, fit, keep = fit_trajectory(s.t, s.xyz['U'], s.sig['U'], steps['TEST'])
    assert abs(p['rate'] - 0.003) < 3e-4 and abs(p['step0'] - 0.01) < 1e-3 and not keep[100] and keep.sum() > 1800
    net = GNSSNetwork([s]).clean(steps=steps).monthly()
    st = net.stations[0]
    assert len(st.monthly_t) == 60 and st.monthly_t[0] == '2016-01'
    # trend and step removed, seasonal kept: the series follows the annual signal
    assert np.std(st.monthly_xyz['U']) > 0.002 and abs(np.mean(st.monthly_xyz['U'])) < 1e-3
    lat, lon, d, sig, comps = net.epoch('2019-04')
    assert comps == ('U', 'N', 'E') and d.shape == (3,) and sig.shape == (3,) and (sig > 0).all()


def test_dates():
    assert _decyear_to_month(2020.0 + 30.5 / 366) == '2020-01'
    assert _decyear_to_month(2020.0 + 31.5 / 366) == '2020-02'
    assert abs(_yymmmdd_to_decyear('20JAN01') - (2020 + 0.5 / 366)) < 1e-9


def test_clean_can_be_rerun(tmp_path):
    f = tmp_path / 'TEST.tenv3'
    write_tenv3(f)
    s = read_tenv3(f)
    net = GNSSNetwork([s])
    once = net.clean().stations[0].xyz['U'].copy()
    twice = net.clean().stations[0].xyz['U']
    assert np.allclose(once, twice)                      # re-running starts from the series as read
    full = net.clean(remove='none').stations[0].xyz['U']
    assert np.std(full) > np.std(once)                   # ... so a different option gives a different result


def test_epoch_common_reference_period(tmp_path):
    import pytest
    write_tenv3(tmp_path / 'A.tenv3', site='AAAA', seed=1)
    write_tenv3(tmp_path / 'B.tenv3', site='BBBB', seed=2)
    a, b = read_tenv3(tmp_path / 'A.tenv3'), read_tenv3(tmp_path / 'B.tenv3')
    net = GNSSNetwork([a, b]).clean().monthly()
    sb = net.stations[1]                                  # B: shift by 5 cm, and lose 2016-2017
    keep = [i for i, k in enumerate(sb.monthly_t) if k >= '2018-01']
    sb.monthly_t = [sb.monthly_t[i] for i in keep]
    for c in sb.monthly_xyz:
        sb.monthly_xyz[c] = sb.monthly_xyz[c][keep] + 0.05
        sb.monthly_sig[c] = sb.monthly_sig[c][keep]
    # default reference = months shared by all stations (2018-2020): both stations refer to the same period
    lat, lon, d, sig, comps = net.epoch('2019-04', components=('U',))
    assert d.shape == (2,) and abs(d[0] - d[1]) < 0.01
    # an explicit reference period works the same way
    _, _, d2, _, _ = net.epoch('2019-04', components=('U',), reference=('2018-01', '2020-12'))
    assert np.allclose(d, d2)
    # no shared month -> clear error instead of silently mixing reference periods
    sb.monthly_t = ['2030-01']
    with pytest.raises(ValueError):
        net.epoch('2019-04', components=('U',))
