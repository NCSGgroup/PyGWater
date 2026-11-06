"""
Consistency tests for PyGWater (run from the project root:  python -m pytest tests -q).
They need data/LLN/PREM-LLNs-complete.dat.
"""
import warnings

import numpy as np
import pytest

from src_auxiliary.constants import EarthConstant
from src_auxiliary.EnumClasses import LLN_Data, LLN_variable, Frame
from src_auxiliary.GeoMathKit import GeoMathKit
from src_EarthModel.LoveNumber import LoveNumber
from src_Basis.GreenFunction.PointLoad import LGF_truncation, LGF_Kummer
from src_Basis.GreenFunction.DiskLoad import DiskLoad
from src_Basis.SphericalHarmonic.Harmonic import Harmonic


def test_constants_consistent():
    """Green's function (a/M) and spherical harmonic / disk (3/rho_e) normalisations must agree (Eq. 11)."""
    a, M, rho_e = EarthConstant.radiusm, EarthConstant.Mass, EarthConstant.rhoear
    assert np.isclose(4 * np.pi * a ** 3 / M, 3 / rho_e, rtol=1e-12)


def test_lln_lazy_and_frames():
    lln = LoveNumber()
    assert lln._LLN is None  # nothing read in __init__
    lln.config(lmax=100, method=LLN_Data.PREM).get_Love_number()
    h1, l1, k1 = (lln.LLN[v][1] for v in (LLN_variable.h, LLN_variable.l, LLN_variable.k))
    assert lln.frame == Frame.CE and k1 == 0
    lln.convert(Frame.CM)
    assert np.isclose(lln.LLN[LLN_variable.k][1], -1.0)
    assert np.isclose(lln.LLN[LLN_variable.h][1], h1 - 1)
    lln.convert(Frame.CF)  # CF: h1 + 2 l1 = 0 (Blewitt, 2003)
    h, l = lln.LLN[LLN_variable.h][1], lln.LLN[LLN_variable.l][1]
    assert np.isclose(h + 2 * l, 0, atol=1e-12)
    lln.convert(Frame.CE)
    assert np.isclose(lln.LLN[LLN_variable.h][1], h1) and np.isclose(lln.LLN[LLN_variable.l][1], l1)


def test_kernels_shapes_and_no_mutation():
    lln = LoveNumber().config(lmax=200, method=LLN_Data.PREM).get_Love_number()
    theta = np.array([0.0, 0.5, 10.0, 90.0])
    keep = theta.copy()
    for K in (LGF_truncation, LGF_Kummer):
        k = K(lln).configure(residence=theta)
        for f in (k.getVertical, k.getHorizental, k.getGeoidheight):
            assert f().shape == (4,)
        assert np.ndim(K(lln).configure(residence=5.0).getVertical()) == 0
    np.testing.assert_array_equal(theta, keep)  # caller's array untouched


def test_disk_lmax_rule_warning():
    lln = LoveNumber().config(lmax=1000, method=LLN_Data.PREM).get_Love_number()
    assert DiskLoad.required_lmax(0.1) == 10000 and DiskLoad.required_lmax(0.05) == 14400
    with pytest.warns(UserWarning):
        DiskLoad(lln, radius=0.5)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        d = DiskLoad(lln, radius=0.5, check_lmax=False).configure(residence=np.array([0.1, 1.0]))
        assert d.getVertical().shape == (2,)


def test_harmonic_nyquist_and_least_squares():
    lat, lon = GeoMathKit.global_equal_distance_grid(2)
    h = Harmonic(lat, lon, 90, option=1)
    c, s = np.zeros((91, 91)), np.zeros((91, 91))
    s[90, 90] = 1.0
    assert np.isclose(h.analysis(h.synthesis(c, s), method='midpoint')[1][90, 90], 1.0)  # was 2.0 before

    for res, L in ((2, 44), (1, 89)):  # SaGEA / Driscoll-Healy weights: exact for lmax <= nlat/2 - 1
        lat_, lon_ = GeoMathKit.global_equal_distance_grid(res)
        h = Harmonic(lat_, lon_, L, option=1)
        rng = np.random.default_rng(1)
        c, s = np.tril(rng.normal(size=(L + 1, L + 1))), np.tril(rng.normal(size=(L + 1, L + 1)))
        s[:, 0] = 0
        c2, s2 = h.analysis(h.synthesis(c, s))
        assert np.allclose(c2, c, atol=1e-10) and np.allclose(s2, s, atol=1e-10)

    h = Harmonic(lat, lon, 89, option=1)
    rng = np.random.default_rng(0)
    c, s = rng.normal(size=(90, 90)), rng.normal(size=(90, 90))
    c, s = np.tril(c), np.tril(s)
    s[:, 0] = 0
    c2, s2 = h.analysis(h.synthesis(c, s), method='least_squares')
    assert np.allclose(c2, c, atol=1e-10) and np.allclose(s2, s, atol=1e-10)


def test_slepian_cap_and_equivalence():
    """Shannon number of a polar cap, orthonormality, and exact equivalence with SH when Lc = (L+1)^2 (Eq. 30)."""
    from src_Basis.Slepian.Slepian import SlepianBasis
    L, res, cap = 15, 0.5, 30.
    lat, lon = GeoMathKit.global_equal_distance_grid(res)
    lon, lat = np.meshgrid(lon, lat)
    sb = SlepianBasis(L).configure(mask=(lat > 90 - cap), lat=lat, lon=lon, resolution=res)
    assert np.isclose(sb.shannon, (L + 1) ** 2 * (1 - np.cos(np.deg2rad(cap))) / 2, rtol=1e-3)
    assert np.allclose(sb.G.T @ sb.G, np.eye(sb.ncoef), atol=1e-10)

    lln = LoveNumber().config(lmax=L, method=LLN_Data.PREM).get_Love_number()
    rng = np.random.default_rng(0)
    c, s = np.tril(rng.normal(size=(L + 1, L + 1))), np.tril(rng.normal(size=(L + 1, L + 1)))
    s[:, 0] = 0
    plat, plon = np.array([70., 80., 65.]), np.array([0., 120., -60.])
    u_sh = (sb.cs2d_to_vec(c, s) * sb._displacement_factor(lln)) @ sb.Ymatrix(plat, plon).T
    assert np.allclose(sb.forward(c, s, lln, plat, plon), u_sh, rtol=1e-10, atol=1e-14)
    cc, ss = sb.vec_to_cs2d(sb.cs2d_to_vec(c, s))
    assert np.allclose(cc, c) and np.allclose(ss, s)

    sb.truncate(beta_min=0.5)
    A = sb.design_matrix(plat, plon, lln)
    assert A.shape == (3, sb.Lc)


def test_harmonic_auto_method():
    """'auto' = 'dh' within the sampling limit, 'midpoint' beyond it (no warning)."""
    lat, lon = GeoMathKit.global_equal_distance_grid(2)
    rng = np.random.default_rng(3)
    g = rng.normal(size=(90, 180))
    for L, ref in ((44, 'dh'), (90, 'midpoint')):
        h = Harmonic(lat, lon, L, option=1)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            c_auto, s_auto = h.analysis(g)
        c_ref, s_ref = h.analysis(g, method=ref)
        assert np.array_equal(c_auto, c_ref) and np.array_equal(s_auto, s_ref)


def _ring(lat0, lon0, dist_deg, n=8):
    az = np.deg2rad(np.arange(n) * 360 / n)
    d, p0 = np.deg2rad(dist_deg), np.deg2rad(lat0)
    plat = np.arcsin(np.sin(p0) * np.cos(d) + np.cos(p0) * np.sin(d) * np.cos(az))
    plon = lon0 + np.rad2deg(np.arctan2(np.sin(az) * np.sin(d) * np.cos(p0), np.cos(d) - np.sin(p0) * np.sin(plat)))
    return np.rad2deg(plat), plon


def test_north_east_projection_point_and_disk():
    """North/East from LGF and D-LGF: constant magnitude on a ring, pointing towards the load, |u_h| = |v|."""
    from src_auxiliary.EnumClasses import Displacement
    from src_Forward.Disk import GFA_displacement as DiskFM
    from src_Forward.PointMass import GFA_displacement as PointFM
    lln = LoveNumber().config(lmax=2000, method=LLN_Data.PREM).get_Love_number()
    lat0, lon0 = 45., 10.
    plat, plon = _ring(lat0, lon0, 5.)
    pts = {'lat': plat, 'lon': plon}
    to_load = GeoMathKit.azimuth(plat, plon, lat0, lon0)
    for FM, grids in ((DiskFM, {'radius': np.array([0.5])}), (PointFM, {'area': np.array([1e9])})):
        grids.update({'lat': np.array([lat0]), 'lon': np.array([lon0]), 'EWH': np.array([[1.0]])})
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fm = FM(lln)
            fm.configure(grids=grids)
            un = fm.evaluation(points=pts, variable=Displacement.North).ravel()
            ue = fm.evaluation(points=pts, variable=Displacement.East).ravel()
            v = fm.evaluation(points=pts, variable=Displacement.Horizontal).ravel()
        mag = np.hypot(un, ue)
        assert np.allclose(mag, np.abs(v), rtol=1e-10) and np.allclose(mag, mag[0], rtol=1e-8)
        assert np.allclose((un * np.cos(to_load) + ue * np.sin(to_load)) / mag, 1.0)  # towards the load (v < 0)


def test_slepian_all_variables_equal_sh():
    """With Lc = (L+1)^2 Slepian = SH for Vertical, Geoheight, North and East (Harmonic gradient synthesis)."""
    from src_auxiliary.EnumClasses import Displacement
    from src_Basis.Slepian.Slepian import SlepianBasis
    import src_auxiliary.EnumClasses as Enums
    from src_Basis.SphericalHarmonic.DataClass import SHC
    L, res = 20, 5.
    lat, lon = GeoMathKit.global_equal_distance_grid(res)
    lon2, lat2 = np.meshgrid(lon, lat)
    sb = SlepianBasis(L).configure(mask=(lat2 > 30), lat=lat2, lon=lon2, resolution=res)
    lln = LoveNumber().config(lmax=L, method=LLN_Data.PREM).get_Love_number()
    rng = np.random.default_rng(5)
    c, s = np.tril(rng.normal(size=(L + 1, L + 1))) * 0.01, np.tril(rng.normal(size=(L + 1, L + 1))) * 0.01
    s[:, 0] = 0
    shc = SHC(c, s, LLN=lln)
    PD = Enums.PhysicalDimensions
    for var, pd in ((Displacement.Vertical, PD.VerticalDisplacement), (Displacement.Geoheight, PD.Geoid),
                    (Displacement.North, PD.HorizontalDisplacementNorth),
                    (Displacement.East, PD.HorizontalDisplacementEast)):
        ref = np.asarray(shc.synthesis(grid_space=res, from_type=PD.EWH, to_type=pd).value).reshape(lat2.shape)
        out = sb.forward(c, s, lln, lat2.ravel(), lon2.ravel(), variable=var).reshape(lat2.shape)
        assert np.allclose(out, ref, rtol=1e-10, atol=1e-14), var


def test_gravity_kernel_vs_wang2012_table():
    """Kummer surface-gravity Green's function reproduces Wang et al. (2012) PREM-LGFs.dat (g2 column)."""
    from src_auxiliary.FileTool import FileTool
    f = FileTool.get_data_dir('LLN/PREM-LGFs.dat')
    if not f.exists():
        pytest.skip('PREM-LGFs.dat not available')
    T = np.loadtxt(f, skiprows=1)
    th = T[:, 0].copy()
    th[-1] -= 0.01
    lln = LoveNumber().config(lmax=40000, method=LLN_Data.PREM).get_Love_number()
    g = LGF_Kummer(lln).configure(th).getGravity(include_degree0=False)
    ratio = g * EarthConstant.radiusm * np.deg2rad(th) * 1e18 / T[:, 4]
    assert np.all(np.abs(ratio - 1) < 2e-3)


def test_observables_sh_slepian_gravity_and_los():
    from src_auxiliary.EnumClasses import Displacement
    from src_EarthModel import Observables
    from src_Basis.Slepian.Slepian import SlepianBasis
    from src_Basis.SphericalHarmonic.DataClass import SHC
    L, res = 20, 5.
    lat, lon = GeoMathKit.global_equal_distance_grid(res)
    lon2, lat2 = np.meshgrid(lon, lat)
    lln = LoveNumber().config(lmax=L, method=LLN_Data.PREM).get_Love_number()
    rng = np.random.default_rng(7)
    c, s = np.tril(rng.normal(size=(L + 1, L + 1))) * 0.01, np.tril(rng.normal(size=(L + 1, L + 1))) * 0.01
    s[:, 0] = 0
    shc = SHC(c, s, LLN=lln)
    sb = SlepianBasis(L)
    sb.G, sb.Lc = np.eye(sb.ncoef), sb.ncoef  # identity basis = SH, tests the Slepian evaluation path
    for var, kw in ((Displacement.Gravity, {'newtonian': 'below'}), (Displacement.Gravity, {'newtonian': 'above'}),
                    (Displacement.GravityDisturbance, {'height': 450e3})):
        ref = np.asarray(shc.synthesis_observable(res, var, lln=lln, **kw).value).reshape(lat2.shape)
        out = sb.forward(c, s, lln, lat2.ravel(), lon2.ravel(), variable=var, **kw).reshape(lat2.shape)
        assert np.allclose(out, ref, rtol=1e-10, atol=1e-16), (var, kw)

    '''above - below = 4 pi G rho_w * EWH (band-limited)'''
    from src_Basis.SphericalHarmonic.Harmonic import Harmonic
    ewh = Harmonic(lat, lon, L, option=1).synthesis(c, s)
    d = (shc.synthesis_observable(res, Displacement.Gravity, lln=lln, newtonian='above').value
         - shc.synthesis_observable(res, Displacement.Gravity, lln=lln, newtonian='below').value)
    assert np.allclose(np.asarray(d).reshape(ewh.shape), 4 * np.pi * EarthConstant.ggg * EarthConstant.rhow * ewh,
                       atol=1e-15)

    '''Sentinel-1 ascending (heading -12, incidence 39): LOS points up and west'''
    look = Observables.look_vector(39., -12.)
    assert np.isclose(np.linalg.norm(look), 1) and look[0] < -0.6 and look[2] > 0.77
