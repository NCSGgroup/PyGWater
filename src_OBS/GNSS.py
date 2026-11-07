"""
GNSS station time series: reading, trajectory model, outliers, monthly means, and the data vector for inversion.

Typical use (NGL / UNR time series, Blewitt et al., 2018):

    from src_OBS.GNSS import read_tenv3, read_ngl_steps, GNSSNetwork
    steps = read_ngl_steps('steps.txt')                         # optional: equipment / earthquake offsets
    net = GNSSNetwork([read_tenv3(f) for f in files])
    net.clean(steps=steps, remove='trend+steps', outlier=5.0)   # keeps seasonal + residual (the loading signal)
    net.monthly()                                               # monthly means with standard errors
    lat, lon, d, sigma, comps = net.epoch('2020-01')            # stacked U, N, E for ForwardModel / inversion

Conventions
-----------
- Units: metres (displacements, sigmas), decimal years (time), degrees (positions).
- Components are returned in the order of COMPONENTS = ('U', 'N', 'E'), matching
  EnumClasses.Displacement.Vertical / North / East in src_Inversion and src_Forward.
- Reference frame: NGL series are in IGS14/IGS20, whose origin is the centre of mass for the long-term trend but
  close to the centre of figure (CF) for seasonal loading. Compare with load Love numbers converted to the CF frame
  (LoveNumber.convert(Frame.CF)) or estimate a degree-1 / network translation in the inversion.
- Common-mode error: removing a network-wide common mode also removes part of the regional loading signal, so
  common_mode() is provided but not applied by clean().
"""
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

COMPONENTS = ('U', 'N', 'E')
_MONTHS = {m: i + 1 for i, m in enumerate(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT',
                                           'NOV', 'DEC'])}


@dataclass
class StationSeries:
    name: str
    lat: float                      # [deg]
    lon: float                      # [deg], -180..180
    t: np.ndarray                   # decimal year
    xyz: dict                       # {'U': (n,), 'N': (n,), 'E': (n,)} [m]
    sig: dict                       # same keys, formal sigmas [m]
    model: dict = field(default_factory=dict)   # trajectory-model parameters per component, after clean()
    monthly_t: list = None          # 'YYYY-MM' after monthly()
    monthly_xyz: dict = None
    monthly_sig: dict = None


def _yymmmdd_to_decyear(s):
    """'20JAN15' -> decimal year (mid-day)"""
    yy, mon, dd = int(s[:2]), _MONTHS[s[2:5].upper()], int(s[5:7])
    year = 2000 + yy if yy < 50 else 1900 + yy
    doy = (np.datetime64('%04d-%02d-%02d' % (year, mon, dd)) - np.datetime64('%04d-01-01' % year)).astype(int)
    ndays = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    return year + (doy + 0.5) / ndays


def _decyear_to_month(t):
    """decimal year -> 'YYYY-MM' (calendar month of that day)"""
    year = int(np.floor(t))
    ndays = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    day = np.datetime64('%04d-01-01' % year) + np.timedelta64(int((t - year) * ndays), 'D')
    return str(day.astype('datetime64[M]'))


def read_tenv3(path, name=None):
    """
    Read a NGL/UNR .tenv3 file: site YYMMMDD yyyy.yyyy __MJD week d reflon _e0(m) __east(m) ____n0(m) _north(m)
    u0(m) ____up(m) _ant(m) sig_e(m) sig_n(m) sig_u(m) __corr_en __corr_eu __corr_nu _latitude(deg)
    _longitude(deg) __height(m). East/north/up are the integer part (e0, n0, u0) plus the fraction.
    :return: StationSeries (positions relative to the first epoch)
    """
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or p[0].lower() == 'site' or len(p) < 23:
                continue
            rows.append(p)
    if not rows:
        raise ValueError('no data in %s' % path)
    a = np.array([[float(v) for v in r[2:23]] for r in rows])
    t = a[:, 0]
    e, n, u = a[:, 5] + a[:, 6], a[:, 7] + a[:, 8], a[:, 9] + a[:, 10]
    lat, lon = float(np.median(a[:, 18])), float(np.median(a[:, 19]))
    lon = (lon + 180) % 360 - 180
    xyz = {'U': u - u[0], 'N': n - n[0], 'E': e - e[0]}
    sig = {'U': a[:, 14], 'N': a[:, 13], 'E': a[:, 12]}
    return StationSeries(name or rows[0][0], lat, lon, t, xyz, sig)


def read_ngl_steps(path):
    """
    NGL steps file ('site YYMMMDD code ...', code 1 = equipment change, 2 = earthquake).
    :return: {site: [decimal years]}
    """
    steps = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3 and re.fullmatch(r'\d{2}[A-Za-z]{3}\d{2}', p[1]):
                steps.setdefault(p[0], []).append(_yymmmdd_to_decyear(p[1]))
    return {k: sorted(set(v)) for k, v in steps.items()}


def trajectory_design(t, steps=(), annual=True, semiannual=True, t0=None):
    """
    Columns of the standard trajectory model: bias, rate, Heaviside steps, annual and semi-annual sine/cosine.
    :return: G (n, p), list of column names
    """
    t = np.asarray(t, dtype=float)
    t0 = np.mean(t) if t0 is None else t0
    cols, names = [np.ones_like(t), t - t0], ['bias', 'rate']
    for k, ts in enumerate(steps):
        if t.min() < ts < t.max():
            cols.append((t >= ts).astype(float)); names.append('step%d' % k)
    for per, on in ((1.0, annual), (0.5, semiannual)):
        if on:
            w = 2 * np.pi / per
            cols += [np.sin(w * (t - t0)), np.cos(w * (t - t0))]
            names += ['sin%g' % per, 'cos%g' % per]
    return np.column_stack(cols), names


def fit_trajectory(t, y, sig=None, steps=(), annual=True, semiannual=True, outlier=5.0, max_iter=5):
    """
    Weighted least-squares trajectory fit with iterative outlier rejection (|residual| > outlier * 1.4826 MAD).
    :return: params (dict name -> value), fitted model (n,), keep mask (n,)
    """
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    w = np.ones_like(y) if sig is None else 1.0 / np.maximum(np.asarray(sig, dtype=float), 1e-6)
    G, names = trajectory_design(t, steps, annual, semiannual)
    keep = np.isfinite(y)
    for _ in range(max_iter):
        x, *_ = np.linalg.lstsq(G[keep] * w[keep, None], y[keep] * w[keep], rcond=None)
        r = y - G @ x
        if outlier is None:
            break
        mad = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new = np.isfinite(y) & (np.abs(r) <= outlier * max(mad, 1e-12))
        if np.array_equal(new, keep):
            break
        keep = new
    return dict(zip(names, x)), G @ x, keep


class GNSSNetwork:

    def __init__(self, stations):
        self.stations = list(stations)

    def clean(self, steps=None, remove='trend+steps', outlier=5.0, annual=True, semiannual=True, min_span=2.0):
        """
        Fit the trajectory model per station and component, flag outliers and remove the chosen terms.
        :param steps: {site: [decimal years]} (e.g. read_ngl_steps); None = no offsets
        :param remove: 'trend+steps' (keep seasonal + residual, i.e. the loading signal; default),
                       'trend+steps+seasonal' (keep the residual only) or 'none'
        :param min_span: drop stations whose series is shorter than this [years]
        """
        keep_st = []
        for s in self.stations:
            if s.t.max() - s.t.min() < min_span:
                continue
            ok = np.ones_like(s.t, dtype=bool)
            for c in COMPONENTS:
                p, fit, k = fit_trajectory(s.t, s.xyz[c], s.sig[c], (steps or {}).get(s.name, ()), annual,
                                           semiannual, outlier)
                s.model[c] = p
                ok &= k
                drop = ['bias', 'rate'] + [n for n in p if n.startswith('step')]
                if remove == 'trend+steps+seasonal':
                    drop += [n for n in p if n.startswith(('sin', 'cos'))]
                if remove != 'none':
                    G, names = trajectory_design(s.t, (steps or {}).get(s.name, ()), annual, semiannual)
                    s.xyz[c] = s.xyz[c] - G[:, [names.index(n) for n in drop]] @ np.array([p[n] for n in drop])
            s.t = s.t[ok]
            for c in COMPONENTS:
                s.xyz[c], s.sig[c] = s.xyz[c][ok], s.sig[c][ok]
            keep_st.append(s)
        self.stations = keep_st
        return self

    def monthly(self, min_days=10):
        """
        Calendar-month means; sigma = max(std / sqrt(n), mean formal sigma / sqrt(n)). Months with fewer than
        min_days daily solutions are dropped.
        """
        for s in self.stations:
            keys = np.array([_decyear_to_month(x) for x in s.t])
            s.monthly_t, s.monthly_xyz, s.monthly_sig = [], {c: [] for c in COMPONENTS}, {c: [] for c in COMPONENTS}
            for k in sorted(set(keys)):
                sel = keys == k
                n = sel.sum()
                if n < min_days:
                    continue
                s.monthly_t.append(k)
                for c in COMPONENTS:
                    v = s.xyz[c][sel]
                    s.monthly_xyz[c].append(v.mean())
                    s.monthly_sig[c].append(max(v.std(ddof=1) / np.sqrt(n), s.sig[c][sel].mean() / np.sqrt(n)))
            for c in COMPONENTS:
                s.monthly_xyz[c], s.monthly_sig[c] = np.array(s.monthly_xyz[c]), np.array(s.monthly_sig[c])
        return self

    def common_mode(self, month_key=None):
        """
        Network common mode per month and component (median over stations) of the monthly series. Not removed
        automatically: it also contains the regional loading signal.
        :return: {month: {comp: value}}
        """
        months = sorted({k for s in self.stations for k in (s.monthly_t or [])})
        out = {}
        for m in months if month_key is None else [month_key]:
            out[m] = {}
            for c in COMPONENTS:
                v = [s.monthly_xyz[c][s.monthly_t.index(m)] for s in self.stations if m in s.monthly_t]
                out[m][c] = float(np.median(v)) if v else np.nan
        return out

    def epoch(self, month, components=COMPONENTS, demean=True):
        """
        Data vector for one month, stacked component by component (as ForwardModel.design_matrix stacks rows).
        :param demean: subtract each station's mean over all months (the inversion then estimates the load
                       anomaly with respect to the same mean as, e.g., GRACE anomalies)
        :return: lat, lon (stations with data), d (n_comp * N,), sigma (n_comp * N,), components
        """
        st = [s for s in self.stations if s.monthly_t and month in s.monthly_t]
        lat = np.array([s.lat for s in st]); lon = np.array([s.lon for s in st])
        d, sg = [], []
        for c in components:
            for s in st:
                i = s.monthly_t.index(month)
                v = s.monthly_xyz[c][i] - (s.monthly_xyz[c].mean() if demean else 0.0)
                d.append(v); sg.append(s.monthly_sig[c][i])
        return lat, lon, np.array(d), np.array(sg), tuple(components)

    def months(self):
        return sorted({k for s in self.stations for k in (s.monthly_t or [])})
