"""
Terrestrial water storage from GNSS, compared with GRACE, for one region (settings in config.py):

  1. preprocessing  : read the NGL series, offsets (steps file + detection), trajectory model, quality control,
                      monthly means, neighbour test, common reference period          -> qc.csv, dataset.npz
  2. inversion      : monthly regularised inversion with the chosen basis                 -> inversion.npz
  3. comparison     : regional mean EWH and mass, GNSS against GRACE mascons             -> series.csv, summary.txt

    python studies/GNSS_GRACE/download.py     # once, needs internet
    python studies/GNSS_GRACE/run.py
"""
import time

import numpy as np

import pygwater as pgw
from pygwater.config import FileTool
from pygwater.preprocessing import GNSSNetwork, GNSSDataset, read_ngl_steps
from pygwater.inversion.workflow import MonthlyInversion
from pygwater.observations.grace import Mascon, compare

from config import CONFIG as C


def decyear(month):
    return int(month[:4]) + (int(month[5:7]) - 1) / 12


def preprocess(out):
    src = FileTool.get_input_dir('gnss/NGL')
    files = sorted((src / 'tenv3').glob('*.tenv3'))
    if not files:
        raise FileNotFoundError('no .tenv3 files in %s; run download.py first' % (src / 'tenv3'))
    steps = read_ngl_steps(src / 'steps.txt', as_events=True) if (src / 'steps.txt').exists() else None
    net = GNSSNetwork.from_files(files).select(box=C['box'])
    net.clean(steps=steps, period=(decyear(C['start']), decyear(C['end']) + 1 / 12), min_years=C['min_years'])
    net.monthly()
    if C['spatial_filter_km']:
        net.spatial_filter(radius_km=C['spatial_filter_km'])
    ds = net.to_dataset(components=C['components'], reference=C['reference'])
    net.save_qc(out / 'qc.csv')
    ds.save(out / 'dataset.npz')
    print('stations: %d read, %d kept; %d months' % (len(files), len(ds.names), len(ds.months)))
    return ds


def main():
    out = FileTool.get_output_dir('GNSS_GRACE/%s' % C['name'])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ds = GNSSDataset.load(out / 'dataset.npz') if (out / 'dataset.npz').exists() and C.get('reuse', False) \
        else preprocess(out)

    lat_min, lat_max, lon_min, lon_max = C['box']
    r, b = C['resolution'], C['buffer']
    lat1 = np.arange(lat_min - b + r / 2, lat_max + b, r)
    lon1 = np.arange(lon_min - b + r / 2, lon_max + b, r)
    lln = pgw.LoveNumber().config(lmax=C['lmax'], method=pgw.LLN_Data.Wang).get_Love_number()
    lln.convert(pgw.Frame.CF)
    opts = dict(subgrid=r / 4) if C['kernel'] == 'cluster_disk' else {}
    cells = pgw.GreenFunctionBasis.from_grid(lat1, lon1, r, kernel=C['kernel'], **opts)
    R = pgw.regularization.curvature(cells.lat, cells.lon, r)
    inv = MonthlyInversion(cells, lln, ds, components=C['components'], R=R)
    alpha = inv.choose_alpha(method=C['alpha_method'])
    res = inv.run(alpha, lb=0.0 if C['nonnegative'] else None)
    res.save(out / 'inversion.npz')
    print('inversion: %d cells, alpha %.2e, %d months (%.0f s)' % (cells.n_param, alpha, len(res.months),
                                                                 time.time() - t0))

    glat = np.arange(lat_min + 0.125, lat_max, 0.25)
    glon = np.arange(lon_min + 0.125, lon_max, 0.25)
    gnss = res.regional_series(glat, glon)
    mass = res.regional_mass(glat, glon)
    rows = ['month,gnss_ewh_mm,gnss_mass_gt,grace_ewh_mm,n_stations']
    grace_file = FileTool.get_input_dir(C['grace_file'])
    if grace_file.exists():
        g = Mascon.read(grace_file).anomalies(C['reference'])
        gs = dict(zip(g.months, g.regional_series(glat, glon)))
        stats = compare(res.months, gnss, g.months, g.regional_series(glat, glon))
        summary = ('GNSS vs GRACE over %s, %d common months: r = %.2f, RMS = %.1f mm, NSE = %.2f, amplitude ratio '
                   '= %.2f' % (C['box'], stats['n'], stats['r'], 1e3 * stats['rms'], stats['nse'],
                               stats['amplitude_ratio']))
    else:
        gs, summary = {}, 'no GRACE file at %s: comparison skipped' % grace_file
    for m, e, ms, n in zip(res.months, gnss, mass, res.n_stations):
        rows.append('%s,%.2f,%.2f,%s,%d' % (m, 1e3 * e, ms, '%.2f' % (1e3 * gs[m]) if m in gs else '', n))
    (out / 'series.csv').write_text('\n'.join(rows) + '\n')
    (out / 'summary.txt').write_text(summary + '\n')
    print(summary)
    print('results in %s' % out)


if __name__ == '__main__':
    main()
