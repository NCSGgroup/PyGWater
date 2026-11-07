"""
Tests of the common Basis interface (src_Basis/base.py) and the basis-independent ForwardModel.
"""
import numpy as np
import pytest

from src_auxiliary.EnumClasses import LLN_Data, Displacement
from src_auxiliary.GeoMathKit import GeoMathKit
from src_auxiliary.constants import EarthConstant
from src_EarthModel.LoveNumber import LoveNumber
from src_EarthModel import Observables
from src_Basis.base import Basis
from src_Basis.SphericalHarmonic.SHBasis import SphericalHarmonicBasis
from src_Basis.Slepian.Slepian import SlepianBasis
from src_Basis.GreenFunction.GFBasis import GreenFunctionBasis
from src_Forward.ForwardModel import ForwardModel
from src_Forward import PointMass, Disk

RES = 2.0
LMAX = 90
STA_LAT = np.array([44.0, 47.5, 52.0, 30.0])
STA_LON = np.array([3.0, 6.5, 12.0, -20.0])


@pytest.fixture(scope='module')
def lln():
    return LoveNumber().config(lmax=LMAX, method=LLN_Data.PREM).get_Love_number()


@pytest.fixture(scope='module')
def load():
    lat, lon = GeoMathKit.global_equal_distance_grid(RES)
    ewh = np.zeros((len(lat), len(lon)))
    ewh[np.ix_((lat > 40) & (lat < 50), (lon > 0) & (lon < 10))] = 0.1
    return ewh, lat, lon


def test_interface():
    for B in (SphericalHarmonicBasis, SlepianBasis, GreenFunctionBasis):
        assert issubclass(B, Basis)


def test_sh_basis_equals_grid_synthesis(lln, load):
    """point-wise SH basis == SaGEA-style grid synthesis of the same coefficients"""
    ewh, lat, lon = load
    sh = SphericalHarmonicBasis(LMAX)
    x = sh.load_parameters(ewh, lat, lon)
    cq, sq = sh.vec_to_cs2d(x)
    i, j = [10, 66, 70], [5, 92, 95]
    for v in (Displacement.Vertical, Displacement.North, Displacement.East, Displacement.Geoheight):
        grid = Observables.sh_synthesis(cq, sq, lln, v, lat, lon)
        pts = sh.predict(x, lat[i], lon[j], lln, v)
        np.testing.assert_allclose(pts, grid[i, j], atol=1e-12)


def test_slepian_full_equals_sh(lln, load):
    ewh, lat, lon = load
    lo2, la2 = np.meshgrid(lon, lat)
    sb = SlepianBasis(30).configure(mask=(la2 > 30) & (la2 < 60) & (lo2 > -10) & (lo2 < 30), lat=lat, lon=lon)
    sh = SphericalHarmonicBasis(30)
    lln30 = LoveNumber().config(lmax=30, method=LLN_Data.PREM).get_Love_number()
    d_sh = sh.predict(sh.load_parameters(ewh, lat, lon), STA_LAT, STA_LON, lln30)
    d_sl = sb.predict(sb.load_parameters(ewh, lat, lon), STA_LAT, STA_LON, lln30)
    np.testing.assert_allclose(d_sl, d_sh, atol=1e-12)
    assert sb.n_param == sb.ncoef and sb.truncate(beta_min=0.1).n_param < sb.ncoef


def test_point_mass_basis_equals_convolution(lln, load):
    ewh, lat, lon = load
    lo2, la2 = np.meshgrid(lon, lat)
    mask = ewh > 0
    gf = GreenFunctionBasis.from_grid(lat, lon, RES, mask=mask, kernel='point_mass')
    x = gf.load_parameters(ewh, lat, lon)
    np.testing.assert_allclose(x, 0.1)
    grids = {'lat': la2[mask], 'lon': lo2[mask],
             'area': np.cos(np.deg2rad(la2[mask])) * np.deg2rad(RES) ** 2 * EarthConstant.radiusm ** 2,
             'EWH': (ewh[mask] * EarthConstant.rhow)[:, None]}
    for v in (Displacement.Vertical, Displacement.North, Displacement.East):
        ref = PointMass.GFA_displacement(lln).configure(grids).evaluation({'lat': STA_LAT, 'lon': STA_LON}, v)[:, 0]
        np.testing.assert_allclose(gf.predict(x, STA_LAT, STA_LON, lln, v), ref, rtol=1e-12, atol=1e-15)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_disk_basis_equals_convolution(load):
    ewh, lat, lon = load
    lln = LoveNumber().config(lmax=2000, method=LLN_Data.PREM).get_Love_number()
    lo2, la2 = np.meshgrid(lon, lat)
    mask = ewh > 0
    gf = GreenFunctionBasis.from_grid(lat, lon, RES, mask=mask, kernel='disk', check_lmax=False, tabulate=False)
    x = gf.load_parameters(ewh, lat, lon)
    grids = {'lat': la2[mask], 'lon': lo2[mask], 'radius': gf.radius, 'EWH': ewh[mask][:, None]}
    pts = {'lat': STA_LAT, 'lon': STA_LON}
    for v in (Displacement.Vertical, Displacement.North, Displacement.East):
        ref = Disk.GFA_displacement(lln).configure(grids).evaluation(pts, v)[:, 0]
        np.testing.assert_allclose(gf.predict(x, STA_LAT, STA_LON, lln, v), ref, rtol=1e-10)
        tab = GreenFunctionBasis.from_grid(lat, lon, RES, mask=mask, kernel='disk', check_lmax=False)  # tabulated
        assert np.abs(tab.predict(x, STA_LAT, STA_LON, lln, v) - ref).max() < 1e-4 * np.abs(ref).max()


def test_forward_model_bases_agree(lln, load):
    """same load, same Earth model: SH and finite point-mass LGF agree (Fig. 2 of Yang et al., 2025)"""
    ewh, lat, lon = load
    obs = (Displacement.Vertical, Displacement.North, Displacement.East)
    out = {}
    for basis in (SphericalHarmonicBasis(LMAX),
                  GreenFunctionBasis.from_grid(lat, lon, RES, mask=ewh > 0, kernel='point_mass')):
        fm = ForwardModel(basis, lln)
        out[basis.name] = fm.run(fm.parameterize(ewh, lat, lon), STA_LAT, STA_LON, observables=obs)
    for v in obs:
        np.testing.assert_allclose(out['spherical_harmonics'][v], out['green_function'][v], atol=5e-5)  # 0.05 mm
    A = ForwardModel(SphericalHarmonicBasis(LMAX), lln).design_matrix(STA_LAT, STA_LON, observables=obs)
    assert A.shape == (3 * STA_LAT.size, (LMAX + 1) ** 2)


def test_los_and_load_synthesis(lln, load):
    ewh, lat, lon = load
    gf = GreenFunctionBasis.from_grid(lat, lon, RES, mask=ewh > 0, kernel='point_mass')
    x = gf.load_parameters(ewh, lat, lon)
    np.testing.assert_allclose(gf.load_synthesis(x, [45.0, 45.0, -30.0], [5.0, 9.9, 100.0]), [0.1, 0.1, 0.0])
    look = Observables.look_vector(39.0, -12.0)
    los = gf.design_matrix_los(STA_LAT, STA_LON, lln, look) @ x
    e, n, u = (gf.predict(x, STA_LAT, STA_LON, lln, v) for v in
               (Displacement.East, Displacement.North, Displacement.Vertical))
    np.testing.assert_allclose(los, Observables.line_of_sight(e, n, u, look), atol=1e-15)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cluster_disk_basis_equals_templates(tmp_path):
    """CD-LGF basis (tabulated kernels, any points) == direct sum of sub-disks in src_Forward/DiskCluster.py"""
    from src_Forward.DiskCluster import cluster_disk
    from src_Basis.GreenFunction.DiskLoad import DiskLoad, DiskLoadTable
    lln = LoveNumber().config(lmax=3000, method=LLN_Data.PREM).get_Love_number()

    psi = np.r_[np.linspace(0, 0.5, 400)[1:], np.linspace(0.5, 180, 400)]
    T = DiskLoadTable(lln, [0.1])
    d = DiskLoad(lln, 0.1).configure(residence=psi)
    for kind, f in (('vertical', d.getVertical), ('geoid', d.getGeoidheight), ('horizontal', d.getHorizental)):
        ref = f()
        assert np.abs(T.evaluate(kind, 0.1, psi) - ref).max() < 5e-5 * np.abs(ref).max()

    pts = {'lat': np.array([45.0, 45.3, 45.9, 47.0, 50.0, 30.0]), 'lon': np.array([5.0, 5.1, 4.6, 6.0, 10.0, -20.0])}
    gf = GreenFunctionBasis([45.0, 46.0], [5.0, 5.0], 'cluster_disk', resolution=1.0, subgrid=0.25)
    for v in (Displacement.Vertical, Displacement.Geoheight, Displacement.North, Displacement.East):
        ref = []
        for c in (45.0, 46.0):
            cd = cluster_disk(lln)
            cd.configure({'lat': c, 'lon': 5.0, 'delta_lat': 1.0, 'delta_lon': 1.0}, small_disk=0.25)
            ref.append(cd.template(pts, v, out_dir=tmp_path))
        ref = np.array(ref).T
        A = gf.design_matrix(pts['lat'], pts['lon'], lln, v)
        assert np.abs(A - ref).max() < 5e-5 * np.abs(ref).max()
    with pytest.raises(NotImplementedError):
        gf.design_matrix(pts['lat'], pts['lon'], lln, Displacement.Gravity)
