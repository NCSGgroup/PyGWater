"""
GNSS station time series: the StationSeries container, readers (NGL .tenv3) and offsets (NGL steps file).

Units: metres (positions, sigmas), decimal years (time), degrees (station coordinates). Components are kept in the
order COMPONENTS = ('U', 'N', 'E'), matching enums.Displacement.Vertical / North / East.
"""
import re
from dataclasses import dataclass, field

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
    raw: tuple = None               # (t, xyz, sig) as read; clean() always starts from it, so it can be re-run
    qc: dict = field(default_factory=dict)      # quality-control figures and flags, after clean()
    detected_jumps: list = field(default_factory=list)   # offsets found by clean(detect_jumps=True) [decimal years]


def _is_leap(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def _yymmmdd_to_decyear(s):
    """'20JAN15' -> decimal year (mid-day)"""
    yy, mon, dd = int(s[:2]), _MONTHS[s[2:5].upper()], int(s[5:7])
    year = 2000 + yy if yy < 50 else 1900 + yy
    doy = (np.datetime64('%04d-%02d-%02d' % (year, mon, dd)) - np.datetime64('%04d-01-01' % year)).astype(int)
    return year + (doy + 0.5) / (366 if _is_leap(year) else 365)


def _decyear_to_month(t):
    """decimal year -> 'YYYY-MM' (calendar month of that day)"""
    year = int(np.floor(t))
    day = np.datetime64('%04d-01-01' % year) + np.timedelta64(int((t - year) * (366 if _is_leap(year) else 365)), 'D')
    return str(day.astype('datetime64[M]'))


def month_to_decyear(month):
    """'YYYY-MM' -> decimal year of the middle of that month"""
    y, m = int(month[:4]), int(month[5:7])
    start = np.datetime64('%04d-%02d-01' % (y, m))
    end = start.astype('datetime64[M]') + 1
    mid = (start - np.datetime64('%04d-01-01' % y)).astype(int) + ((end.astype('datetime64[D]') - start).astype(int)) / 2
    return y + mid / (366 if _is_leap(y) else 365)


def read_tenv3(path, name=None):
    """
    Read a NGL/UNR .tenv3 file (Blewitt et al., 2018): site YYMMMDD yyyy.yyyy MJD week d reflon e0 east n0 north
    u0 up ant sig_e sig_n sig_u corr_en corr_eu corr_nu latitude longitude height. East/north/up are the integer
    part (e0, n0, u0) plus the fraction.
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
    order = np.argsort(a[:, 0], kind='stable')
    a = a[order]
    t = a[:, 0]
    e, n, u = a[:, 5] + a[:, 6], a[:, 7] + a[:, 8], a[:, 9] + a[:, 10]
    lat, lon = float(np.median(a[:, 18])), float(np.median(a[:, 19]))
    lon = (lon + 180) % 360 - 180
    xyz = {'U': u - u[0], 'N': n - n[0], 'E': e - e[0]}
    sig = {'U': a[:, 14], 'N': a[:, 13], 'E': a[:, 12]}
    return StationSeries(name or rows[0][0], lat, lon, t, xyz, sig)


@dataclass
class Step:
    t: float                        # decimal year
    kind: str                       # 'equipment' or 'earthquake'
    info: str = ''                  # equipment change type, or 'distance_km magnitude event'


def read_ngl_steps(path, kinds=('equipment', 'earthquake'), as_events=False):
    """
    NGL steps file (https://geodesy.unr.edu/NGLStationPages/steps.txt). Two line types:
        site YYMMMDD 1 change_type                                    equipment change
        site YYMMMDD 2 threshold_km distance_km magnitude event_id    earthquake within the threshold distance
    :param kinds: which offsets to keep
    :param as_events: False -> {site: [decimal years]} (one step per date, for the trajectory model);
                      True  -> {site: [Step, ...]} with kind and details
    """
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 3 or not re.fullmatch(r'\d{2}[A-Za-z]{3}\d{2}', p[1]):
                continue
            kind = 'earthquake' if p[2] == '2' else 'equipment'
            if kind not in kinds:
                continue
            info = ' '.join(p[4:7]) if kind == 'earthquake' else ' '.join(p[3:])
            out.setdefault(p[0], []).append(Step(_yymmmdd_to_decyear(p[1]), kind, info))
    if as_events:
        return {k: sorted(v, key=lambda s: s.t) for k, v in out.items()}
    return {k: sorted({s.t for s in v}) for k, v in out.items()}
