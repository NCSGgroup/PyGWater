"""
Synthetic checkerboard test of Fu, Argus & Landerer (2015), JGR Solid Earth 120, 552-566, Sect. 4.2 / Fig. 11a-b,
reproduced with PyGWater, and extended to other bases and to three-component data.

Fu et al. (2015) setup (Sects. 2.2, 4.2), followed here:
    domain      : 39.5-51.5 N, Pacific coast to 114.5 W (Washington and Oregon + 2.5 deg buffer), land only
    parameters  : 0.25 x 0.25 deg cells, each a uniform disk (D-LGF), PREM Green's functions (Wang et al., 2012)
    data        : GPS vertical displacement only (seasonal change, 1 Apr - 1 Oct)
    stations    : the ~200 continuous GPS sites of their Fig. 1, digitised from the vector figure of the paper
                  (positions accurate to about 0.01 deg; <data root>/input/stations/fu2015_stations_digitized.txt)
    inversion   : min ||W (G x - b)||^2 + beta^2 ||L x||^2, L = Laplacian between neighbouring cells, bounded
                  (x >= 0) least squares (Stark & Parker, 1995); beta from the misfit-roughness trade-off curve
    checkerboard: 0.75 deg blocks (3 x 3 cells) of 500 mm and 0 mm over Washington and Oregon (42-49 N, west of
                  116.5 W), read from their Fig. 11a; zero elsewhere

Additions in this test:
    truth       : generated (a) with the D-LGF itself, as in a standard checkerboard test like theirs, and (b) with
                  CD-LGF (0.05 deg sub-disks, N = 26 000), i.e. a uniform load over each square cell, so that the
                  inversion operator is not the one that made the data
    noise       : none, as in the paper, and sigma_U = 1 mm (seasonal displacements from multi-year series)
    bases       : D-LGF (their choice), CD-LGF 0.125 and point masses (finite LGF, n = 720)
    components  : U only (their choice) and U + N + E (sigma_N = sigma_E = sigma_U / 3)
    alpha       : L-curve corner (their trade-off-curve choice) with BVLS (x >= 0, theirs) and unbounded; GCV; and
                  the oracle alpha (smallest true error), to separate the parameter choice from the basis
    blocks      : 0.75 deg (their test) and 1.5 deg

Run from any folder:  python studies/Fu2015_PNW/synthetic.py
Results: <data root>/output/inversion_PNW_Fu2015/results.npz and summary.csv
"""
import os
import sys
import time
import warnings
from pathlib import Path


import numpy as np

from pygwater.config import FileTool
from pygwater.utils.enums import LLN_Data, Displacement
from pygwater.earth.love_numbers import LoveNumber
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.inversion.solver import RegularizedLeastSquares
from pygwater.inversion.regularization import curvature
from pygwater.observations.synthetic import add_noise

RES = 0.25
DOMAIN = (39.5, 51.5, -126.0, -114.5)            # lat_min, lat_max, lon_min, lon_max (west edge: land/coast)
WA_OR = (42.0, 49.0, -126.0, -116.5)            # where the checkerboard is placed
BLOCKS = (0.75, 1.5)
AMPLITUDE = 0.5                                  # 500 mm in the "on" blocks, 0 in the "off" blocks
SIGMA_U = {'none': 1e-3, '1 mm': 1e-3}          # weights; 'none' adds no noise
H_RATIO = 1 / 3
COMPONENTS = {'U': (Displacement.Vertical,), 'UNE': (Displacement.Vertical, Displacement.North, Displacement.East)}


def land_mask(lat, lon, shp=None):
    import geopandas as gpd
    import shapely
    shp = FileTool.get_input_dir('shp/NorthAmerica/NorthAmerica.shp') if shp is None else shp
    return shapely.contains_xy(gpd.read_file(shp).union_all(), lon, lat)


def checker(lat, lon, block):
    """Fu et al. (2015) Fig. 11a: 0 / AMPLITUDE blocks inside Washington-Oregon, zero elsewhere"""
    i = np.floor((lat - WA_OR[0]) / block).astype(int)
    j = np.floor((lon - WA_OR[2]) / block).astype(int)
    on = ((i + j) % 2 == 0) & (lat >= WA_OR[0]) & (lat < WA_OR[1]) & (lon < WA_OR[3])
    return np.where(on, AMPLITUDE, 0.0)


def setup():
    """cells (land, 0.25 deg), stations (digitised from Fu et al., 2015, Fig. 1), Love numbers"""
    lat1 = np.arange(DOMAIN[0] + RES / 2, DOMAIN[1], RES)
    lon1 = np.arange(DOMAIN[2] + RES / 2, DOMAIN[3], RES)
    lo2, la2 = np.meshgrid(lon1, lat1)
    land = land_mask(la2, lo2)
    st = np.loadtxt(FileTool.get_input_dir('stations/fu2015_stations_digitized.txt'))
    lln = {n: LoveNumber().config(lmax=n, method=LLN_Data.PREM).get_Love_number() for n in (720, 13000, 26000)}
    return lat1, lon1, land, la2[land], lo2[land], st[:, 0], st[:, 1], lln


def make_bases(lat1, lon1, land, clat, clon, lln):
    """inversion bases: the D-LGF of Fu et al. (2015) and two alternatives"""
    return {'D-LGF (Fu)': (GreenFunctionBasis.from_grid(lat1, lon1, RES, mask=land, kernel='disk'), lln[13000]),
            'CD-LGF 0.125': (GreenFunctionBasis(clat, clon, 'cluster_disk', resolution=RES, subgrid=0.125),
                             lln[13000]),
            'finite LGF 720': (GreenFunctionBasis.from_grid(lat1, lon1, RES, mask=land, kernel='point_mass'),
                               lln[720])}


def invert(solver, d, alphas, x_true=None, mask=None):
    """L-curve alpha with BVLS (x >= 0, Fu et al.), L-curve and GCV unbounded, and the oracle alpha (synthetic)"""
    out = {}
    a_l, curve = solver.choose_alpha(d, alphas=alphas, method='lcurve')
    a_g, _ = solver.choose_alpha(d, alphas=alphas, method='gcv')
    out['lcurve+bvls'] = (a_l, solver.solve_bounded(d, a_l, lb=0.0).x)
    out['lcurve'] = (a_l, solver.solve(d, a_l).x)
    out['gcv'] = (a_g, solver.solve(d, a_g).x)
    if x_true is not None:
        errs = [np.sqrt(np.mean((solver.solve(d, a).x - x_true)[mask] ** 2)) for a in alphas]
        a_o = float(alphas[int(np.argmin(errs))])
        out['oracle'] = (a_o, solver.solve(d, a_o).x)
    return out, curve


def run(out_dir=None, seed=7):
    out_dir = FileTool.get_output_dir('inversion_PNW_Fu2015') if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter('ignore', UserWarning)
    t0 = time.time()
    lat1, lon1, land, clat, clon, slat, slon, lln = setup()
    ns = slat.size
    print('%d land cells of %.2f deg, %d stations' % (clat.size, RES, ns))
    comps_all = COMPONENTS['UNE']

    t = time.time()
    bases = make_bases(lat1, lon1, land, clat, clon, lln)
    A = {b: {c: basis.design_matrix(slat, slon, ln, c) for c in comps_all} for b, (basis, ln) in bases.items()}
    cd05 = GreenFunctionBasis(clat, clon, 'cluster_disk', resolution=RES, subgrid=0.05)
    truth_ops = {'D-LGF (as Fu)': A['D-LGF (Fu)'],
                 'CD-LGF 0.05': {c: cd05.design_matrix(slat, slon, lln[26000], c) for c in comps_all}}
    R = curvature(clat, clon, RES)
    print('operators %.0f s' % (time.time() - t))

    rows, est = [], {}
    west = clon < -121.0
    waor = (clat >= WA_OR[0]) & (clat < WA_OR[1]) & (clon < WA_OR[3])
    alphas = np.logspace(-6, 5, 34)
    for blk in BLOCKS:
        x_true = checker(clat, clon, blk)
        for tname, G in truth_ops.items():
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
                        res, _ = invert(RegularizedLeastSquares(Am, sig, R), d, alphas, x_true, waor)
                        for method, (alpha, rec) in res.items():
                            e = rec - x_true
                            sel = {'all': waor, 'W': waor & west, 'E': waor & ~west}
                            row = dict(block=blk, truth=tname, noise=nname, basis=bname, comps=cname, method=method,
                                       alpha=alpha)
                            for k, m in sel.items():
                                row['rms_' + k + '_mm'] = 1e3 * np.sqrt(np.mean(e[m] ** 2))
                                row['r_' + k] = float(np.corrcoef(rec[m], x_true[m])[0, 1])
                            rows.append(row)
                            est[(blk, tname, nname, bname, cname, method)] = rec
                        r = rows[-4:]
                        print('%.2f %-13s %-4s %-14s %-3s | ' % (blk, tname, nname, bname, cname) + '  '.join(
                            '%s r %.2f/%.2f' % (q['method'], q['r_W'], q['r_E']) for q in r))
    keys = list(est)
    np.savez(out_dir / 'results.npz', cell_lat=clat, cell_lon=clon, station_lat=slat, station_lon=slon,
             blocks=np.array(BLOCKS), truths=np.array([checker(clat, clon, b) for b in BLOCKS]),
             est_keys=np.array(['|'.join(map(str, k)) for k in keys]), est=np.array([est[k] for k in keys]))
    cols = list(rows[0])
    (out_dir / 'summary.csv').write_text(','.join(cols) + '\n' + '\n'.join(
        ','.join(('%.4g' % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) for r in rows) + '\n')
    print('total %.0f s' % (time.time() - t0))
    return rows


if __name__ == '__main__':
    run()
