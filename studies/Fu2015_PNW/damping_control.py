"""
Control run for slepian.py: the cell bases of synthetic.py with the
regulariser swapped from the Laplacian (roughness) to area-weighted damping (R = diag(cell area), i.e. minimum
norm), so that Slepian (damping) and cells can be compared with the same kind of prior.
Results: <data root>/output/inversion_PNW_Fu2015/damping_summary.csv
"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np
from pygwater.config import FileTool
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.inversion.solver import RegularizedLeastSquares
from pygwater.inversion.regularization import area_weighted_identity
from synthetic import setup, make_bases, checker, RES, WA_OR, BLOCKS, SIGMA_U, H_RATIO, COMPONENTS
from slepian import write_csv
from pygwater.observations.synthetic import add_noise
from pygwater.utils.enums import Displacement


def run(seed=7):
    warnings.simplefilter('ignore', UserWarning)
    t0 = time.time()
    lat1, lon1, land, clat, clon, slat, slon, lln = setup()
    ns, comps_all = slat.size, COMPONENTS['UNE']
    bases = make_bases(lat1, lon1, land, clat, clon, lln)
    A = {b: {c: basis.design_matrix(slat, slon, ln, c) for c in comps_all} for b, (basis, ln) in bases.items()}
    cd05 = GreenFunctionBasis(clat, clon, 'cluster_disk', resolution=RES, subgrid=0.05)
    G = {c: cd05.design_matrix(slat, slon, lln[26000], c) for c in comps_all}
    R = area_weighted_identity(clat)
    R = R / np.mean(np.diag(R))
    west = clon < -121.0
    waor = (clat >= WA_OR[0]) & (clat < WA_OR[1]) & (clon < WA_OR[3])
    alphas = np.logspace(-8, 5, 40)
    rows, est = [], {}
    for blk in BLOCKS:
        x_true = checker(clat, clon, blk)
        for nname in SIGMA_U:
            d_all = {c: G[c] @ x_true for c in comps_all}
            if nname != 'none':
                d_all = {c: add_noise(d_all[c], SIGMA_U[nname] * (1 if c == Displacement.Vertical else H_RATIO),
                                      seed=seed + k) for k, c in enumerate(comps_all)}
            for bname in bases:
                for cname, comps in COMPONENTS.items():
                    Am = np.vstack([A[bname][c] for c in comps])
                    sig = np.concatenate([np.full(ns, SIGMA_U[nname] * (1 if c == Displacement.Vertical
                                                                          else H_RATIO)) for c in comps])
                    d = np.concatenate([d_all[c] for c in comps])
                    sv = RegularizedLeastSquares(Am, sig, R)
                    choice = {m: sv.choose_alpha(d, alphas=alphas, method=m)[0] for m in ('lcurve', 'gcv')}
                    errs = [np.sqrt(np.mean((sv.solve(d, a).x - x_true)[waor] ** 2)) for a in alphas]
                    choice['oracle'] = float(alphas[int(np.argmin(errs))])
                    line = []
                    for method, alpha in choice.items():
                        rec = sv.solve(d, alpha).x
                        row = dict(block=blk, truth='CD-LGF 0.05', noise=nname, basis=bname, comps=cname,
                                   regulariser='damping', method=method, alpha=alpha)
                        for k, m in {'all': waor, 'W': waor & west, 'E': waor & ~west}.items():
                            row['rms_' + k + '_mm'] = 1e3 * np.sqrt(np.mean((rec - x_true)[m] ** 2))
                            row['r_' + k] = float(np.corrcoef(rec[m], x_true[m])[0, 1])
                        rows.append(row)
                        est['|'.join(map(str, (blk, nname, bname, cname, method)))] = rec
                        line.append('%s %.2f/%.2f' % (method, row['r_W'], row['r_E']))
                    print('%.2f %-4s %-14s %-3s | %s' % (blk, nname, bname, cname, '  '.join(line)))
    out = FileTool.get_output_dir('inversion_PNW_Fu2015')
    write_csv(out / 'damping_summary.csv', rows)
    np.savez(out / 'damping_results.npz', cell_lat=clat, cell_lon=clon, est_keys=np.array(list(est)),
             est=np.array(list(est.values())))
    print('total %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    run()
