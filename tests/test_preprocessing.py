"""Tests of the GNSS preprocessing steps (synthetic NGL files, no external data needed)."""
import numpy as np
import pytest

from pygwater.preprocessing import GNSSNetwork, GNSSDataset, read_tenv3, read_ngl_steps
from pygwater.preprocessing import ngl
from pygwater.preprocessing.trajectory import clean_segments

MON = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def write_series(path, site, lat=46.0, lon=7.0, start='2012-01-01', end='2019-01-01', annual=4e-3, phase=0.0,
                 rate=2e-3, jumps=(), noise=1.5e-3, seed=0, skip=None, extra=None):
    """daily NGL .tenv3 file; jumps = [(decimal year, size m)], skip(t) -> bool drops days, extra(t) -> m added"""
    rng = np.random.default_rng(seed)
    lines = ['site YYMMMDD yyyy.yyyy __MJD week d reflon _e0(m) __east(m) ____n0(m) _north(m) u0(m) ____up(m) _ant(m) '
             'sig_e(m) sig_n(m) sig_u(m) __corr_en __corr_eu __corr_nu _latitude(deg) _longitude(deg) __height(m)\n']
    for d in np.arange(np.datetime64(start), np.datetime64(end)):
        y, m, dd = int(str(d)[:4]), int(str(d)[5:7]), int(str(d)[8:10])
        t = y + ((d - np.datetime64('%d-01-01' % y)).astype(int) + 0.5) / (366 if y % 4 == 0 else 365)
        if skip is not None and skip(t):
            continue
        up = rate * (t - 2015) + annual * np.sin(2 * np.pi * t + phase) + rng.normal(0, noise)
        up += sum(s for tj, s in jumps if t >= tj) + (extra(t) if extra else 0.0)
        lines.append('%s %02d%s%02d %.4f 57000 2000 0 %.1f 1 0.0 2 0.0 3 %.5f 0.0 0.001 0.001 0.003 0 0 0 %.4f %.4f '
                     '100.0\n' % (site, y % 100, MON[m - 1], dd, t, lon, up, lat, lon))
    path.write_text(''.join(lines))
    return path


def test_unreported_jump_is_detected(tmp_path):
    f = write_series(tmp_path / 'JUMP.tenv3', 'JUMP', jumps=[(2016.4, 0.012)])
    net = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=True).monthly()
    s = net.stations[0]
    assert s.qc['n_detected_jumps'] == 1 and abs(s.detected_jumps[0] - 2016.4) < 5 / 365
    # after removing it, what is left is the 4 mm annual signal
    t = np.array([int(k[:4]) + (int(k[5:7]) - 0.5) / 12 for k in s.monthly_t])
    resid = s.monthly_xyz['U'] - 4e-3 * np.sin(2 * np.pi * t)
    assert np.std(resid) < 1.5e-3
    # without detection the 12 mm jump stays in the series
    s2 = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=False, max_residual=1).monthly().stations[0]
    resid2 = s2.monthly_xyz['U'] - 4e-3 * np.sin(2 * np.pi * t)
    assert np.std(resid2) > 2 * np.std(resid)


def test_segments_and_close_offsets():
    t = np.concatenate([2010 + np.arange(30) / 365.25, 2012 + np.arange(1500) / 365.25])
    keep, st = clean_segments(t, steps=[2013.0, 2013.02])
    assert not keep[:30].any() and keep[30:].sum() < 1500       # isolated month dropped, data between offsets too
    assert st == [2013.0]


def test_minimum_data_and_completeness(tmp_path):
    short = write_series(tmp_path / 'SHRT.tenv3', 'SHRT', start='2016-01-01', end='2018-01-01')
    gappy = write_series(tmp_path / 'GAPY.tenv3', 'GAPY', start='2010-01-01', end='2019-01-01',
                         skip=lambda t: (t * 365.25) % 3 > 0.9)     # one day in three
    net = GNSSNetwork([read_tenv3(short), read_tenv3(gappy)]).clean(detect_jumps=False)
    assert net.stations == [] and set(net.rejected) == {'SHRT', 'GAPY'}


def test_loading_correction_is_subtracted(tmp_path):
    atm = lambda t: 3e-3 * np.cos(2 * np.pi * t * 4)                  # a 3-month "atmospheric" signal
    f = write_series(tmp_path / 'CORR.tenv3', 'CORR', extra=atm, noise=0.5e-3)
    t = np.arange(2011, 2020, 0.01)
    corr = {'CORR': (t, {'U': atm(t)})}
    a = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=False, annual=True).stations[0]
    b = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=False, corrections=corr).stations[0]
    assert b.qc['loading_corrected'] and b.qc['rms_U'] < 0.5 * a.qc['rms_U']


def test_monthly_sigma_rms_is_more_conservative(tmp_path):
    f = write_series(tmp_path / 'SIGM.tenv3', 'SIGM')
    rms = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=False).monthly(sigma='rms').stations[0]
    sem = GNSSNetwork([read_tenv3(f)]).clean(detect_jumps=False).monthly(sigma='sem').stations[0]
    # here the standard error is floored by the 3 mm formal sigma / sqrt(n); the RMS is still about 3 times larger
    assert np.median(rms.monthly_sig['U']) > 2.5 * np.median(sem.monthly_sig['U'])


def test_spatial_filter_rejects_station_moving_against_its_neighbours(tmp_path):
    files = []
    for k in range(6):
        files.append(write_series(tmp_path / ('N%03d.tenv3' % k), 'N%03d' % k, lat=46 + 0.3 * k, lon=7 + 0.2 * k,
                                  seed=k))
    # a poroelastic site: rises when the region subsides
    files.append(write_series(tmp_path / 'PORO.tenv3', 'PORO', lat=46.5, lon=7.5, phase=np.pi, seed=9))
    files.append(write_series(tmp_path / 'FAR0.tenv3', 'FAR0', lat=10, lon=100, seed=10))   # no neighbours
    net = GNSSNetwork.from_files(files).clean(detect_jumps=False).monthly().spatial_filter(radius_km=500)
    names = net.names()
    assert 'PORO' not in names and 'against' in net.rejected['PORO']
    assert len([n for n in names if n.startswith('N')]) == 6 and 'FAR0' in names
    far = [s for s in net.stations if s.name == 'FAR0'][0]
    assert far.qc['spatial_filter'] == 'too few neighbours'


def test_dataset_reference_and_roundtrip(tmp_path):
    a = write_series(tmp_path / 'AAAA.tenv3', 'AAAA', seed=1)
    b = write_series(tmp_path / 'BBBB.tenv3', 'BBBB', start='2014-01-01', seed=2)
    net = GNSSNetwork.from_files([a, b]).clean(detect_jumps=False).monthly()
    ds = net.to_dataset(components=('U',))
    assert ds.reference[0] >= '2014-01'                                # months shared by both stations
    lat, lon, d, sig, comps = ds.epoch('2013-06')                      # only A has data then
    assert d.shape == (1,)
    ds.save(tmp_path / 'ds.npz')
    ds2 = GNSSDataset.load(tmp_path / 'ds.npz')
    assert ds2.months == ds.months and np.allclose(np.nan_to_num(ds2.values['U']), np.nan_to_num(ds.values['U']))
    net.save_qc(tmp_path / 'qc.csv')
    assert 'rms_U' in (tmp_path / 'qc.csv').read_text().splitlines()[0]


def test_ngl_catalogue_and_steps(tmp_path):
    (tmp_path / 'DataHoldings.txt').write_text(
        'Sta  Lat(deg)   Long(deg) Hgt(m)  X(m) Y(m) Z(m) Dtbeg Dtend Dtmod NumSol StaOrigName\n'
        'AAAA  46.0  7.0 100.0 1 2 3 2010-01-01 2020-01-01 2025-01-01 3600\n'
        'BBBB  46.5  367.5 100.0 1 2 3 2018-01-01 2020-01-01 2025-01-01 700\n'
        'CCCC  10.0  100.0 100.0 1 2 3 2005-01-01 2020-01-01 2025-01-01 5000\n')
    st = ngl.read_holdings(tmp_path / 'DataHoldings.txt')
    assert [s.name for s in st] == ['AAAA', 'BBBB', 'CCCC'] and abs(st[1].lon - 7.5) < 1e-9
    assert ngl.select_stations(st, box=(45, 47, 6, 8), start='2012-01', end='2019-12', min_years=4) == ['AAAA']
    (tmp_path / 'steps.txt').write_text('AAAA  15MAR01  1  Antenna_Type_Changed\n'
                                        'AAAA  16JUN10  2  0.5  12.3  6.4  us1000abcd\n')
    ev = read_ngl_steps(tmp_path / 'steps.txt', as_events=True)['AAAA']
    assert [e.kind for e in ev] == ['equipment', 'earthquake'] and ev[1].info.split()[1] == '6.4'
    assert len(read_ngl_steps(tmp_path / 'steps.txt', kinds=('equipment',))['AAAA']) == 1
