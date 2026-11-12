"""
Forward-model comparison of basis functions over China (cf. Yang et al., 2025, applied regionally).

Load   : CSR RL06.3 mascon, January 2020 (time index 180), averaged/sampled to a global 2-deg grid
         (<data root>/input/mascon/mascon.txt: lat, lon, EWH [cm]); the WHOLE globe is loaded, so far-field loads count.
Points : 1-deg grid nodes inside China (polygon <data root>/input/shp/China/China.shp, or any shapefile given).
Model  : PREM load Love numbers (Wang et al., 2012), CE frame.
Output : vertical, north and east displacement at the points for each basis; saved to
         <data root>/output/China_basis_comparison/results.npz (+ a summary printed and written as CSV).

Bases (all through the common Basis interface of pygwater.basis):
    SH90            spherical harmonics, L = 90 (analysis 'auto'), = Nyquist of the 2-deg load
    finite LGF 90   point masses, Love numbers to n = 90  (equivalent to SH90, Sect. 4)
    finite LGF 360  point masses, Love numbers to n = 360 (truncation beyond the load resolution)
    D-LGF           one equal-area disk per 2-deg cell, N = 13 000 (Eq. 24)
    CD-LGF s        each 2-deg cell split into s x s deg sub-disks (s = 0.5, 0.2, 0.1), N = 13 000;
                    cells farther than FAR_DEG from every point use the D-LGF (identical there: < 2e-5 mm per
                    metre of load per cell beyond 20 deg)
    Slepian         China-concentrated Slepian functions, L = 90, truncated at beta > 0.1 and beta > 0.01
Reference: CD-LGF 0.1 (uniform load per mascon, the configuration recommended by Yang et al., 2025).

Run from any folder:  python studies/China/basis_comparison.py          (about 10 min on 2 cores)
"""
import os
import sys
import time
import warnings
from pathlib import Path


import numpy as np

from pygwater.config import FileTool
from pygwater.utils.enums import LLN_Data, Displacement
from pygwater.utils.geomath import GeoMathKit
from pygwater.earth.love_numbers import LoveNumber
from pygwater.basis.spherical_harmonic.sh_basis import SphericalHarmonicBasis
from pygwater.basis.green.gf_basis import GreenFunctionBasis
from pygwater.basis.slepian import SlepianBasis

RES_LOAD = 2.0          # load grid [deg]
RES_PTS = 1.0           # evaluation grid [deg]
FAR_DEG = 20.0          # CD-LGF beyond this distance -> D-LGF
N_DISK = 13000          # Love-number degree for disk kernels (Eq. 24 for 0.1-deg sub-disks)
COMPONENTS = (Displacement.Vertical, Displacement.North, Displacement.East)


def china_mask(lat, lon, shp=None):
    """True where (lat, lon) (2d arrays) lies inside the China polygon(s) of the shapefile"""
    import geopandas as gpd
    import shapely
    shp = FileTool.get_input_dir('shp/China/China.shp') if shp is None else shp
    geom = gpd.read_file(shp).union_all()
    return shapely.contains_xy(geom, lon, lat)


def load_mascon(file=None):
    """global 2-deg EWH [m] on the cell-centred grid of GeoMathKit"""
    file = FileTool.get_mascon_dir('mascon.txt') if file is None else file
    lat, lon = GeoMathKit.global_equal_distance_grid(RES_LOAD)
    ewh = np.loadtxt(file)[:, 2].reshape(len(lat), len(lon)) * 0.01
    return ewh, lat, lon


def run(shp=None, out_dir=None, slepian_L=90):
    out_dir = FileTool.get_output_dir('China_basis_comparison') if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter('ignore', UserWarning)

    ewh, lat, lon = load_mascon()
    plat1, plon1 = GeoMathKit.global_equal_distance_grid(RES_PTS)
    plo2, pla2 = np.meshgrid(plon1, plat1)
    inside = china_mask(pla2, plo2, shp)
    plat, plon = pla2[inside], plo2[inside]
    print('load %s, %d evaluation points in China' % (ewh.shape, plat.size))

    lln = {n: LoveNumber().config(lmax=n, method=LLN_Data.PREM).get_Love_number() for n in (90, 360, N_DISK)}
    out, timing = {}, {}

    def record(name, basis, x, ln):
        t = time.time()
        out[name] = np.array([basis.predict(x, plat, plon, ln, v) for v in COMPONENTS])  # (3, N) [m]
        timing[name] = time.time() - t
        print('%-16s %6.1f s   vertical range %.2f .. %.2f mm' % (name, timing[name], 1e3 * out[name][0].min(),
                                                                   1e3 * out[name][0].max()))

    '''spherical harmonics'''
    sh = SphericalHarmonicBasis(90)
    record('SH90', sh, sh.load_parameters(ewh, lat, lon), lln[90])

    '''finite LGF (point masses), global'''
    pm = GreenFunctionBasis.from_grid(lat, lon, RES_LOAD, kernel='point_mass')
    xg = pm.load_parameters(ewh, lat, lon)
    for n in (90, 360):
        record('finite LGF %d' % n, pm, xg, lln[n])

    '''D-LGF, global (tabulated kernels)'''
    dk = GreenFunctionBasis.from_grid(lat, lon, RES_LOAD, kernel='disk')
    t = time.time()
    A_disk = np.array([dk.design_matrix(plat, plon, lln[N_DISK], v) for v in COMPONENTS])  # (3, N, ncell)
    out['D-LGF'] = A_disk @ xg
    timing['D-LGF'] = time.time() - t
    print('%-16s %6.1f s' % ('D-LGF', timing['D-LGF']))

    '''CD-LGF: near cells clustered, far cells as D-LGF'''
    dist = np.array([GeoMathKit.angular_distance(dk.lat[j], dk.lon[j], plat, plon).min() for j in range(dk.n_param)])
    near = dist <= FAR_DEG
    far_part = A_disk[:, :, ~near] @ xg[~near]
    for s in (0.5, 0.2, 0.1):
        cd = GreenFunctionBasis(dk.lat[near], dk.lon[near], 'cluster_disk', resolution=RES_LOAD, subgrid=s)
        t = time.time()
        out['CD-LGF %.1f' % s] = np.array([cd.predict(xg[near], plat, plon, lln[N_DISK], v) for v in COMPONENTS]) \
                                 + far_part
        timing['CD-LGF %.1f' % s] = time.time() - t
        print('%-16s %6.1f s   (%d near cells)' % ('CD-LGF %.1f' % s, timing['CD-LGF %.1f' % s], near.sum()))

    '''Slepian functions concentrated on China'''
    t = time.time()
    mlat1, mlon1 = GeoMathKit.global_equal_distance_grid(0.5)
    mlo2, mla2 = np.meshgrid(mlon1, mlat1)
    sb = SlepianBasis(slepian_L).configure(mask=china_mask(mla2, mlo2, shp), lat=mlat1, lon=mlon1, resolution=0.5)
    print('Slepian L=%d: Shannon %.1f, beta>0.1: %d, beta>0.01: %d  (%.0f s)'
          % (slepian_L, sb.shannon, (sb.beta > 0.1).sum(), (sb.beta > 0.01).sum(), time.time() - t))
    vec = sb.sh_coefficients(ewh, lat, lon)
    ln = lln[90] if slepian_L == 90 else LoveNumber().config(lmax=slepian_L, method=LLN_Data.PREM).get_Love_number()
    for b in (0.1, 0.01):
        sb.truncate(beta_min=b)
        record('Slepian b>%g (%d)' % (b, sb.Lc), sb, sb.slepian_coefficients(vec), ln)

    np.savez(out_dir / 'results.npz', lat=plat, lon=plon, components=[c.name for c in COMPONENTS],
             names=list(out), values=np.array(list(out.values())), load=ewh, load_lat=lat, load_lon=lon,
             timing=np.array([timing[k] for k in out]), beta=sb.beta, shannon=sb.shannon)
    summary(out_dir / 'results.npz')
    return out


def summary(npz, ref='CD-LGF 0.1'):
    """RMS and max |difference| [mm] of every basis against the reference, per component; also as CSV"""
    r = np.load(npz)
    names, vals = list(r['names']), r['values'] * 1e3
    R = vals[names.index(ref)]
    lines = ['basis,component,rms_mm,max_mm,signal_rms_mm']
    print('\n%-22s %28s' % ('difference to ' + ref, '   '.join(r['components'])))
    for k, name in enumerate(names):
        d = vals[k] - R
        row = []
        for c, comp in enumerate(r['components']):
            rms, mx = np.sqrt(np.mean(d[c] ** 2)), np.abs(d[c]).max()
            row.append('%5.2f/%5.2f' % (rms, mx))
            lines.append('%s,%s,%.4f,%.4f,%.4f' % (name, comp, rms, mx, np.sqrt(np.mean(R[c] ** 2))))
        print('%-22s %s   (rms/max mm)' % (name, '  '.join(row)))
    (Path(npz).parent / 'summary.csv').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    run()
