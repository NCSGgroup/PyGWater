"""
Real-data reproduction of Fu, Argus & Landerer (2015), Fig. 3: multi-year mean seasonal change in equivalent water
thickness (EWT) from 1 October to 1 April in Washington and Oregon, inverted from GPS vertical displacements.

Data: the seasonal uplift (1 April -> 1 October, mm) of the ~200 continuous GPS sites of their Fig. 1, digitised from
the paper (positions from the vector circles, values from their fill colour on the 0-12 mm colour bar):
<data root>/input/stations/fu2015_seasonal_uplift_digitized.txt. A load INCREASE from October to April causes
SUBSIDENCE, so the data vector is d = -uplift.

Inversion as in Fu et al. (2015) and synthetic.py (same domain, 0.25-deg cells, Laplacian
smoothing, bounded x >= 0), for three bases: D-LGF (theirs), CD-LGF 0.125 and point masses. Regularisation
strength alpha: L-curve corner, GCV, discrepancy principle (sigma = 1.5 mm), and the alpha that best reproduces
their Fig. 3 (to show how much smoothing their beta = 3 implied).
Comparison: their Fig. 3, digitised by hue (hill-shading changes lightness, not hue):
<data root>/input/stations/fu2015_fig3_ewt_digitized.txt (cells with weakly saturated colour, i.e. small values,
are not digitised).

Run from any folder:  python studies/Fu2015_PNW/real_data.py
Results: <data root>/output/inversion_PNW_Fu2015/real_results.npz and real_summary.csv
"""
import os
import sys
import time
import warnings
from pathlib import Path


import numpy as np

from pygwater.config import FileTool
from pygwater.utils.enums import Displacement
from pygwater.utils.constants import EarthConstant
from pygwater.inversion.solver import RegularizedLeastSquares
from pygwater.inversion.regularization import curvature
from synthetic import setup, make_bases, RES, WA_OR

SIGMA_U = 1.5e-3  # [m] assumed error of a site's mean seasonal uplift (GPS + digitising); Fu et al. quote none.
                  # It only matters for the discrepancy principle (and sigma0).


def cell_area(lat):
    return np.cos(np.deg2rad(lat)) * np.deg2rad(RES) ** 2 * EarthConstant.radiusm ** 2


def run(out_dir=None):
    out_dir = FileTool.get_output_dir('inversion_PNW_Fu2015') if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter('ignore', UserWarning)
    t0 = time.time()
    lat1, lon1, land, clat, clon, _, _, lln = setup()
    obs = np.loadtxt(FileTool.get_input_dir('stations/fu2015_seasonal_uplift_digitized.txt'))
    slat, slon, uplift = obs[:, 0], obs[:, 1], obs[:, 2] * 1e-3
    d = -uplift
    fig3 = np.loadtxt(FileTool.get_input_dir('stations/fu2015_fig3_ewt_digitized.txt'))
    key = {(round(a, 3), round(o, 3)): k for k, (a, o) in enumerate(zip(clat, clon))}
    idx = np.array([key.get((round(a, 3), round(o, 3)), -1) for a, o in fig3[:, :2]])
    ok = idx >= 0
    fu_idx, fu_ewt = idx[ok], fig3[ok, 2] * 1e-3
    waor = (clat >= WA_OR[0]) & (clat < WA_OR[1]) & (clon < WA_OR[3])
    print('%d stations, uplift %.1f-%.1f mm; %d cells to compare with Fig. 3' % (len(d), 1e3 * uplift.min(),
                                                                                1e3 * uplift.max(), ok.sum()))
    bases = make_bases(lat1, lon1, land, clat, clon, lln)
    R = curvature(clat, clon, RES)
    alphas = np.logspace(-6, 5, 45)
    rows, est = [], {}
    for bname, (basis, ln) in bases.items():
        A = basis.design_matrix(slat, slon, ln, Displacement.Vertical)
        solver = RegularizedLeastSquares(A, SIGMA_U, R)
        choice = {m: solver.choose_alpha(d, alphas=alphas, method=m)[0] for m in ('lcurve', 'gcv', 'discrepancy')}
        # the alpha that reproduces their Fig. 3 best (needs their result; shows the smoothness they chose)
        errs = [np.sqrt(np.mean((solver.solve(d, a).x[fu_idx] - fu_ewt) ** 2)) for a in alphas]
        choice['matched to Fig. 3'] = float(alphas[int(np.argmin(errs))])
        for method, alpha in choice.items():
            x = solver.solve_bounded(d, alpha, lb=0.0).x      # bounded x >= 0, as Fu et al.
            res = A @ x - d
            mass_gt = np.sum(x[waor] * cell_area(clat[waor])) * EarthConstant.rhow / 1e12
            diff = x[fu_idx] - fu_ewt
            row = dict(basis=bname, method=method, alpha=alpha, dof=solver.solve(d, alpha).dof,
                       fit_rms_mm=1e3 * np.sqrt(np.mean(res ** 2)),
                       ewt_max_mm=1e3 * x[waor].max(), ewt_p95_mm=1e3 * np.percentile(x[waor], 95),
                       ewt_mean_waor_mm=1e3 * x[waor].mean(), mass_waor_gt=mass_gt,
                       vs_fig3_rms_mm=1e3 * np.sqrt(np.mean(diff ** 2)), vs_fig3_bias_mm=1e3 * diff.mean(),
                       vs_fig3_r=float(np.corrcoef(x[fu_idx], fu_ewt)[0, 1]))
            rows.append(row); est[(bname, method)] = x
            print('%-14s %-18s alpha %.1e dof %5.1f fit %.2f mm  EWT p95 %4.0f mean %4.0f mm  %5.1f Gt  '
                  'vs Fig.3: rms %4.0f bias %+4.0f r %.2f' % (bname, method, alpha, row['dof'], row['fit_rms_mm'],
                                                               row['ewt_p95_mm'], row['ewt_mean_waor_mm'], mass_gt,
                                                               row['vs_fig3_rms_mm'], row['vs_fig3_bias_mm'],
                                                               row['vs_fig3_r']))
    keys = list(est)
    np.savez(out_dir / 'real_results.npz', cell_lat=clat, cell_lon=clon, station_lat=slat, station_lon=slon,
             uplift=uplift, fig3_lat=fig3[:, 0], fig3_lon=fig3[:, 1], fig3_ewt=fig3[:, 2] * 1e-3,
             est_keys=np.array(['|'.join(k) for k in keys]), est=np.array([est[k] for k in keys]))
    cols = list(rows[0])
    (out_dir / 'real_summary.csv').write_text(','.join(cols) + '\n' + '\n'.join(
        ','.join(('%.4g' % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) for r in rows) + '\n')
    print('total %.0f s' % (time.time() - t0))
    return rows


if __name__ == '__main__':
    run()
