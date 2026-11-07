"""
Closed-loop and checkerboard tests of the GNSS inversion for terrestrial water storage over China.

Nature run (the "truth")
    load      : (a) CSR RL06.3 mascon, January 2020, 2 deg, inside the parameter domain only ('CSR')
                (b) the same, but the GLOBAL load generates the data ('CSR global'): load outside the domain leaks in
                (c) checkerboards of +/-10 cm EWH with 4, 6 and 8 deg blocks inside the domain
    forward   : CD-LGF with 0.1-deg sub-disks, N = 13 000 (cells outside the domain: D-LGF)
    network   : 260 synthetic stations inside China, >= 1.5 deg apart (Poisson disk; CMONOC-like density)
    noise     : white, three levels sigma_U = 0.2, 1 and 2 mm, with sigma_N = sigma_E = sigma_U / 3
                (monthly GNSS means; the CSR month gives 3 mm RMS vertical, a 10-cm checkerboard 1.2-1.6 mm)

Inversion
    parameters: 2-deg cells within BUFFER deg of China (mascons), or Slepian functions concentrated on that domain
    bases     : CD-LGF 0.2 (consistent with the truth), finite LGF 90 (point masses), Slepian L = 60 (beta > 0.1)
    data      : vertical only ('U') or three components ('UNE')
    solver    : src_Inversion.Solver.RegularizedLeastSquares; R = first-difference smoothness (mascons) or identity
                (Slepian); alpha by GCV on the noisy data (option: 'wgcv', 'lcurve'); the 'oracle' alpha with
                the smallest true error is reported alongside, to separate parameter choice from basis effects;
                safeguard: only alphas giving a recovered EWH RMS below MAX_EWH_RMS (30 cm; CSR gives 11.5 cm) are eligible
    metrics   : recovered vs true EWH at the 2-deg cells inside China: RMS error [cm], correlation

Run from any folder:  python src_demo/demo_inversion_China.py [gcv|wgcv|lcurve]   (about 2 min on 2 cores)
Results: <data root>/res/inversion_China/<alpha method>/results.npz and summary.csv
"""
import os
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE_ROOT = HERE.parent
sys.path.insert(0, str(CODE_ROOT))
os.chdir(HERE)

import numpy as np

from src_auxiliary.FileTool import FileTool
from src_auxiliary.EnumClasses import LLN_Data, Displacement
from src_auxiliary.GeoMathKit import GeoMathKit
from src_EarthModel.LoveNumber import LoveNumber
from src_Basis.GreenFunction.GFBasis import GreenFunctionBasis
from src_Basis.Slepian.Slepian import SlepianBasis
from src_Inversion.Solver import RegularizedLeastSquares
from src_Inversion.Regularization import smoothness, identity
from src_OBS.Synthetic import synthetic_network, add_noise, checkerboard
from demo_China_basis_comparison import china_mask, load_mascon

RES = 2.0
BUFFER = 4.0                    # parameter domain: cells within this distance of China [deg]
NOISE = {'0.2 mm': 2e-4, '1 mm': 1e-3, '2 mm': 2e-3}   # sigma_U [m]; horizontal = sigma_U / 3
H_RATIO = 1 / 3
OMEGA = 1.2                     # weight of the weighted GCV ('wgcv')
MAX_EWH_RMS = 0.3               # plausibility safeguard: monthly EWH anomaly RMS over the domain below 30 cm
COMPONENT_SETS = {'U': (Displacement.Vertical,),
                  'UNE': (Displacement.Vertical, Displacement.North, Displacement.East)}
BLOCKS = (4.0, 6.0, 8.0)


def stack(A_by_comp, comps):
    return np.concatenate([A_by_comp[c] for c in comps], axis=0)


def sigma_of(c, sigma_u):
    return sigma_u if c == Displacement.Vertical else sigma_u * H_RATIO


def sig(comps, n, sigma_u):
    return np.concatenate([np.full(n, sigma_of(c, sigma_u)) for c in comps])


def run(shp=None, out_dir=None, n_stations=260, seed=1, alpha_method='gcv', safeguard=True):
    out_dir = FileTool.get_res_dir('inversion_China/%s%s' % (alpha_method, '' if safeguard else '_nosafeguard')) \
        if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter('ignore', UserWarning)
    t0 = time.time()

    '''grid, China cells and parameter domain'''
    ewh_global, lat, lon = load_mascon()
    lo2, la2 = np.meshgrid(lon, lat)
    in_china = china_mask(la2, lo2, shp)
    cl, co = la2[in_china], lo2[in_china]
    dist = np.array([GeoMathKit.angular_distance(a, o, cl, co).min() for a, o in zip(la2.ravel(), lo2.ravel())])
    domain = (dist <= BUFFER).reshape(la2.shape)
    dlat, dlon = la2[domain], lo2[domain]
    china_in_domain = in_china[domain]
    print('%d China cells, %d domain cells (buffer %g deg)' % (in_china.sum(), domain.sum(), BUFFER))

    '''synthetic network'''
    slat, slon = synthetic_network(lambda a, o: china_mask(a, o, shp), n_target=n_stations, min_dist=1.5,
                                   bbox=(15, 56, 70, 140), seed=seed)
    ns = slat.size
    print('%d stations' % ns)

    lln = {n: LoveNumber().config(lmax=n, method=LLN_Data.PREM).get_Love_number() for n in (60, 90, 13000)}
    comps_all = COMPONENT_SETS['UNE']

    '''truth operator: CD-LGF 0.1 on the domain, D-LGF outside'''
    t = time.time()
    truth_in = GreenFunctionBasis(dlat, dlon, 'cluster_disk', resolution=RES, subgrid=0.1)
    G_in = {c: truth_in.design_matrix(slat, slon, lln[13000], c) for c in comps_all}
    outside = ~domain
    truth_out = GreenFunctionBasis.from_grid(lat, lon, RES, mask=outside, kernel='disk')
    G_out = {c: truth_out.design_matrix(slat, slon, lln[13000], c) for c in comps_all}
    print('truth operators %.0f s' % (time.time() - t))

    '''true loads (on the domain cells) and data'''
    truths = {'CSR': ewh_global[domain], 'CSR global': ewh_global[domain]}
    for b in BLOCKS:
        truths['checker %g' % b] = checkerboard(dlat, dlon, b, 0.1)
    data = {}
    for name, x_true in truths.items():
        d = {c: G_in[c] @ x_true for c in comps_all}
        if name == 'CSR global':
            for c in comps_all:
                d[c] = d[c] + G_out[c] @ ewh_global[outside]
        for nname, su in NOISE.items():
            data[name, nname] = {c: add_noise(d[c], sigma_of(c, su), seed=seed + 10 + k)
                                 for k, c in enumerate(comps_all)}

    '''inversion bases'''
    t = time.time()
    bases = {}
    cd = GreenFunctionBasis(dlat, dlon, 'cluster_disk', resolution=RES, subgrid=0.2)
    bases['CD-LGF 0.2'] = (cd, lln[13000], smoothness(dlat, dlon, RES))
    pm = GreenFunctionBasis.from_grid(lat, lon, RES, mask=domain, kernel='point_mass')   # same cell order
    bases['finite LGF 90'] = (pm, lln[90], smoothness(dlat, dlon, RES))
    mlat, mlon = GeoMathKit.global_equal_distance_grid(0.5)
    mlo2, mla2 = np.meshgrid(mlon, mlat)
    cell_of = lambda a, o: domain[np.clip(((a + 90) // RES).astype(int), 0, domain.shape[0] - 1),
                                  np.clip(((o + 180) // RES).astype(int), 0, domain.shape[1] - 1)]
    sb = SlepianBasis(60).configure(mask=cell_of(mla2, mlo2), lat=mlat, lon=mlon, resolution=0.5)
    sb.truncate(beta_min=0.1)
    print('Slepian L=60 on the domain: Shannon %.1f, %d functions (beta > 0.1)' % (sb.shannon, sb.Lc))
    bases['Slepian L60'] = (sb, lln[60], identity(sb.Lc))
    A = {name: {c: b.design_matrix(slat, slon, ln, c) for c in comps_all} for name, (b, ln, R) in bases.items()}
    print('inversion operators %.0f s' % (time.time() - t))

    '''solve every case'''
    rows, est = [], {}
    for tname, x_true in truths.items():
        truth_china = x_true[china_in_domain]
        for bname, (b, ln, R) in bases.items():
            for cname, comps in COMPONENT_SETS.items():
                for nname, su in NOISE.items():
                    solver = RegularizedLeastSquares(stack(A[bname], comps), sig(comps, ns, su), R)
                    d = np.concatenate([data[tname, nname][c] for c in comps])
                    alphas = np.logspace(-7, 3, 41)
                    plaus = (lambda x, b=b: np.sqrt(np.mean(b.load_synthesis(x, dlat, dlon) ** 2)) < MAX_EWH_RMS) \
                        if safeguard else None
                    alpha, _ = solver.choose_alpha(d, alphas=alphas, method=alpha_method, omega=OMEGA, plausible=plaus)
                    s = solver.solve(d, alpha)
                    '''oracle: the alpha with the smallest true error (only possible in a synthetic test)'''
                    X = np.array([solver.solve(d, a).x for a in alphas])
                    E = b.load_synthesis(X, dlat, dlon)[:, china_in_domain] - truth_china
                    oracle = np.sqrt(np.mean(E ** 2, axis=1))
                    rec = b.load_synthesis(s.x, dlat, dlon)[china_in_domain]  # EWH at the China cell centres
                    err = rec - truth_china
                    rows.append(dict(truth=tname, basis=bname, comps=cname, noise=nname, alpha=alpha, dof=s.dof,
                                     sigma0=s.sigma0, rms_err_cm=100 * np.sqrt(np.mean(err ** 2)),
                                     oracle_err_cm=100 * oracle.min(), oracle_alpha=float(alphas[oracle.argmin()]),
                                     rms_true_cm=100 * np.sqrt(np.mean(truth_china ** 2)),
                                     corr=float(np.corrcoef(rec, truth_china)[0, 1])))
                    est[(tname, bname, cname, nname)] = rec
                    print('%-11s %-13s %-3s %-6s alpha %.1e dof %5.1f s0 %.2f  err %5.2f cm (oracle %5.2f, true %5.2f)'
                          '  r %.2f' % (tname, bname, cname, nname, alpha, s.dof, s.sigma0, rows[-1]['rms_err_cm'],
                                        rows[-1]['oracle_err_cm'], rows[-1]['rms_true_cm'], rows[-1]['corr']))

    keys = list(est)
    np.savez(out_dir / 'results.npz', station_lat=slat, station_lon=slon, cell_lat=dlat[china_in_domain],
             cell_lon=dlon[china_in_domain], domain_lat=dlat, domain_lon=dlon,
             truth_names=list(truths), truths=np.array([truths[k][china_in_domain] for k in truths]),
             est_keys=np.array(['|'.join(k) for k in keys]), est=np.array([est[k] for k in keys]),
             n_slepian=sb.Lc, shannon=sb.shannon)
    cols = list(rows[0])
    (out_dir / 'summary.csv').write_text(','.join(cols) + '\n' + '\n'.join(
        ','.join(('%.4g' % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) for r in rows) + '\n')
    print('total %.0f s' % (time.time() - t0))
    return rows


if __name__ == '__main__':
    run(alpha_method=sys.argv[1] if len(sys.argv) > 1 else 'gcv')
