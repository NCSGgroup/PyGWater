"""
Slepian basis in the Fu et al. (2015) tests (synthetic checkerboard and real data), to set beside the D-LGF, CD-LGF
and point-mass results of synthetic.py and real_data.py.

Slepian functions concentrated on the land cells of the Fu et al. domain (39.5-51.5 N, coast to 114.5 W, 0.25 deg),
band limits L = 120 (resolves the 1.5 deg blocks) and L = 240 (0.75 deg blocks), truncated at the Shannon number
(Simons et al., 2006). Built with SlepianBasis.configure(method='svd').

Parameters are the Slepian coefficients; EWH on the cells is x_cell = S x with S = Y(cells) G[:, :Lc]. Two
regularisers:
    damping   : R = I on the coefficients (the usual Slepian choice: truncation does most of the work)
    Laplacian : R = S^T (L^T L) S, the same roughness penalty on the EWH as the cell bases
alpha: L-curve corner, GCV, oracle (synthetic) and discrepancy principle (real data). No x >= 0 bound (it applies
to cells, not to coefficients).

Run from any folder:  python studies/Fu2015_PNW/slepian.py
Results: <data root>/output/inversion_PNW_Fu2015/slepian_results.npz, slepian_summary.csv, slepian_real_summary.csv
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
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.basis.slepian import SlepianBasis
from pygwater.inversion.solver import RegularizedLeastSquares
from pygwater.inversion.regularization import curvature
from pygwater.observations.synthetic import add_noise
from synthetic import setup, checker, RES, WA_OR, BLOCKS, SIGMA_U, H_RATIO, COMPONENTS
from real_data import SIGMA_U as SIGMA_REAL

LMAX = (120, 240)


def slepian_bases(lat1, lon1, land, clat, clon, lln):
    out = {}
    lo2, la2 = np.meshgrid(lon1, lat1)
    for L in LMAX:
        t = time.time()
        sb = SlepianBasis(L).configure(land.astype(float), la2, lo2, RES, method='svd').truncate(shannon=True)
        sb.G = np.ascontiguousarray(sb.G[:, :sb.Lc])         # keep only the retained functions (memory)
        S = np.vstack([sb.Ymatrix(clat[b], clon[b]) @ sb.G for b in sb._blocks(clat.size)])
        print('Slepian L=%d: Shannon %.1f -> %d functions (%.0f s)' % (L, sb.shannon, sb.Lc, time.time() - t))
        out['Slepian L=%d' % L] = (sb, lln[720], S)
    return out


def solvers(A, sig, S, Rc):
    n = S.shape[1]
    return {'damping': RegularizedLeastSquares(A, sig, np.eye(n)),
            'Laplacian': RegularizedLeastSquares(A, sig, S.T @ Rc @ S)}


def synthetic(lat1, lon1, land, clat, clon, slat, slon, lln, sbases, Rc, seed=7):
    comps_all = COMPONENTS['UNE']
    ns = slat.size
    truth_ops = {'D-LGF (as Fu)': GreenFunctionBasis.from_grid(lat1, lon1, RES, mask=land, kernel='disk'),
                 'CD-LGF 0.05': GreenFunctionBasis(clat, clon, 'cluster_disk', resolution=RES, subgrid=0.05)}
    truth_ops = {k: {c: b.design_matrix(slat, slon, lln[13000 if k.startswith('D') else 26000], c)
                     for c in comps_all} for k, b in truth_ops.items()}
    A = {b: {c: sb.design_matrix(slat, slon, ln, c) for c in comps_all} for b, (sb, ln, _) in sbases.items()}
    west = clon < -121.0
    waor = (clat >= WA_OR[0]) & (clat < WA_OR[1]) & (clon < WA_OR[3])
    alphas = np.logspace(-8, 5, 40)
    rows, est = [], {}
    for blk in BLOCKS:
        x_true = checker(clat, clon, blk)
        for tname, G in truth_ops.items():
            for nname in SIGMA_U:
                d_all = {c: G[c] @ x_true for c in comps_all}
                if nname != 'none':
                    d_all = {c: add_noise(d_all[c], SIGMA_U[nname] * (1 if c == Displacement.Vertical else H_RATIO),
                                          seed=seed + k) for k, c in enumerate(comps_all)}
                for bname, (_, _, S) in sbases.items():
                    for cname, comps in COMPONENTS.items():
                        Am = np.vstack([A[bname][c] for c in comps])
                        sig = np.concatenate([np.full(ns, SIGMA_U[nname] * (1 if c == Displacement.Vertical
                                                                              else H_RATIO)) for c in comps])
                        d = np.concatenate([d_all[c] for c in comps])
                        line = []
                        for rname, sv in solvers(Am, sig, S, Rc).items():
                            choice = {m: sv.choose_alpha(d, alphas=alphas, method=m)[0] for m in ('lcurve', 'gcv')}
                            errs = [np.sqrt(np.mean((S @ sv.solve(d, a).x - x_true)[waor] ** 2)) for a in alphas]
                            choice['oracle'] = float(alphas[int(np.argmin(errs))])
                            for method, alpha in choice.items():
                                rec = S @ sv.solve(d, alpha).x
                                e = rec - x_true
                                row = dict(block=blk, truth=tname, noise=nname, basis=bname, comps=cname,
                                           regulariser=rname, method=method, alpha=alpha)
                                for k, m in {'all': waor, 'W': waor & west, 'E': waor & ~west}.items():
                                    row['rms_' + k + '_mm'] = 1e3 * np.sqrt(np.mean(e[m] ** 2))
                                    row['r_' + k] = float(np.corrcoef(rec[m], x_true[m])[0, 1])
                                rows.append(row)
                                est[(blk, tname, nname, bname, cname, rname, method)] = rec
                                line.append('%s/%s %.2f/%.2f' % (rname[:3], method, row['r_W'], row['r_E']))
                        print('%.2f %-13s %-4s %-14s %-3s | %s' % (blk, tname, nname, bname, cname, '  '.join(line)))
    return rows, est


def real(clat, clon, sbases, Rc):
    obs = np.loadtxt(FileTool.get_input_dir('stations/fu2015_seasonal_uplift_digitized.txt'))
    slat, slon, d = obs[:, 0], obs[:, 1], -obs[:, 2] * 1e-3
    fig3 = np.loadtxt(FileTool.get_input_dir('stations/fu2015_fig3_ewt_digitized.txt'))
    key = {(round(a, 3), round(o, 3)): k for k, (a, o) in enumerate(zip(clat, clon))}
    idx = np.array([key.get((round(a, 3), round(o, 3)), -1) for a, o in fig3[:, :2]])
    fu_idx, fu_ewt = idx[idx >= 0], fig3[idx >= 0, 2] * 1e-3
    waor = (clat >= WA_OR[0]) & (clat < WA_OR[1]) & (clon < WA_OR[3])
    area = np.cos(np.deg2rad(clat)) * np.deg2rad(RES) ** 2 * EarthConstant.radiusm ** 2
    alphas = np.logspace(-8, 5, 53)
    rows, est = [], {}
    for bname, (sb, ln, S) in sbases.items():
        A = sb.design_matrix(slat, slon, ln, Displacement.Vertical)
        for rname, sv in solvers(A, np.full(d.size, SIGMA_REAL), S, Rc).items():
            choice = {m: sv.choose_alpha(d, alphas=alphas, method=m)[0] for m in ('lcurve', 'gcv', 'discrepancy')}
            for method, alpha in choice.items():
                sol = sv.solve(d, alpha)
                x = S @ sol.x
                diff = x[fu_idx] - fu_ewt
                row = dict(basis=bname, regulariser=rname, method=method, alpha=alpha, dof=sol.dof,
                           fit_rms_mm=1e3 * np.sqrt(np.mean((A @ sol.x - d) ** 2)),
                           ewt_min_mm=1e3 * x[waor].min(), ewt_mean_waor_mm=1e3 * x[waor].mean(),
                           mass_waor_gt=np.sum(x[waor] * area[waor]) * EarthConstant.rhow / 1e12,
                           vs_fig3_rms_mm=1e3 * np.sqrt(np.mean(diff ** 2)), vs_fig3_bias_mm=1e3 * diff.mean(),
                           vs_fig3_r=float(np.corrcoef(x[fu_idx], fu_ewt)[0, 1]))
                rows.append(row); est[(bname, rname, method)] = x
                print('real %-14s %-9s %-11s alpha %.1e dof %5.1f fit %.2f mm  min %5.0f mm  %6.1f Gt  '
                      'vs Fig.3 rms %4.0f r %.2f' % (bname, rname, method, alpha, row['dof'], row['fit_rms_mm'],
                                                     row['ewt_min_mm'], row['mass_waor_gt'], row['vs_fig3_rms_mm'],
                                                     row['vs_fig3_r']))
    return rows, est


def write_csv(path, rows):
    cols = list(rows[0])
    path.write_text(','.join(cols) + '\n' + '\n'.join(
        ','.join(('%.4g' % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) for r in rows) + '\n')


def run(out_dir=None):
    out_dir = FileTool.get_output_dir('inversion_PNW_Fu2015') if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter('ignore', UserWarning)
    t0 = time.time()
    lat1, lon1, land, clat, clon, slat, slon, lln = setup()
    sbases = slepian_bases(lat1, lon1, land, clat, clon, lln)
    Rc = curvature(clat, clon, RES)
    rows_r, est_r = real(clat, clon, sbases, Rc)
    write_csv(out_dir / 'slepian_real_summary.csv', rows_r)
    rows_s, est_s = synthetic(lat1, lon1, land, clat, clon, slat, slon, lln, sbases, Rc)
    write_csv(out_dir / 'slepian_summary.csv', rows_s)
    ks, kr = list(est_s), list(est_r)
    np.savez(out_dir / 'slepian_results.npz', cell_lat=clat, cell_lon=clon,
             est_keys=np.array(['|'.join(map(str, k)) for k in ks]), est=np.array([est_s[k] for k in ks]),
             real_keys=np.array(['|'.join(k) for k in kr]), real=np.array([est_r[k] for k in kr]),
             n_functions=np.array([sbases[b][2].shape[1] for b in sbases]))
    print('total %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    if len(sys.argv) > 1:                                  # e.g. ... 360 480 -> results in slepian_L360_480/
        LMAX = tuple(int(a) for a in sys.argv[1:])
        run(FileTool.get_output_dir('inversion_PNW_Fu2015/slepian_L' + '_'.join(sys.argv[1:])))
    else:
        run()
