"""
Nevada Geodetic Laboratory (NGL, University of Nevada, Reno; Blewitt et al., 2018): station list, station selection
and download of daily .tenv3 time series and the steps (offsets) file.

    from pygwater.preprocessing import ngl
    ngl.download_catalogue(out_dir)                                   # DataHoldings.txt and steps.txt
    st = ngl.read_holdings(out_dir / 'DataHoldings.txt')
    sites = ngl.select_stations(st, box=(42, 49, -125, -116), start='2006-01', end='2020-12', min_years=5)
    ngl.download_series(sites, out_dir / 'tenv3')                     # SITE.tenv3, skips files already present

Downloads need internet access to geodesy.unr.edu; run them on your own computer. Files go to
<data root>/input/gnss/NGL/ by default (see the study studies/GNSS_GRACE/download.py).

NGL series are processed in the IGS20 frame (IGS14 for the older archive). The URL of a station's .tenv3 file is
tried in the order of TENV3_URLS; the first that answers with a valid series is used. The IGS20 layout
(gps_timeseries/IGS20/tenv3/IGS20/SITE.tenv3) was checked in November 2026; the region folders next to it (NA/, EU/,
...) hold plate-fixed series, which are not used here.
"""
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BASE = 'https://geodesy.unr.edu'
HOLDINGS_URL = BASE + '/NGLStationPages/DataHoldings.txt'
STEPS_URL = BASE + '/NGLStationPages/steps.txt'
TENV3_URLS = {
    'IGS20': (BASE + '/gps_timeseries/IGS20/tenv3/IGS20/{site}.tenv3',
              BASE + '/gps_timeseries/tenv3/IGS20/{site}.tenv3'),
    'IGS14': (BASE + '/gps_timeseries/IGS14/tenv3/IGS14/{site}.tenv3',
              BASE + '/gps_timeseries/tenv3/IGS14/{site}.tenv3'),
}


@dataclass
class StationInfo:
    name: str
    lat: float
    lon: float                      # -180..180
    height: float
    start: str                      # 'YYYY-MM-DD'
    end: str
    n_solutions: int


def read_holdings(path):
    """
    DataHoldings.txt: Sta Lat(deg) Long(deg) Hgt(m) X Y Z Dtbeg Dtend Dtmod NumSol [StaOrigName]
    :return: list of StationInfo
    """
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 11 or p[0] == 'Sta':
                continue
            try:
                lon = (float(p[2]) + 180) % 360 - 180
                out.append(StationInfo(p[0], float(p[1]), lon, float(p[3]), p[7], p[8], int(p[10])))
            except ValueError:
                continue
    return out


def select_stations(stations, box=None, inside=None, start=None, end=None, min_years=3.0, min_solutions=0):
    """
    Stations in a region whose data cover at least `min_years` within [start, end].
    :param box: (lat_min, lat_max, lon_min, lon_max) [degree]
    :param inside: or a function (lat, lon) -> bool array (e.g. a shapefile test)
    :param start, end: 'YYYY-MM' or 'YYYY-MM-DD' (None = open)
    :return: list of station names
    """
    s0 = np.datetime64(start or '1990-01', 'D')
    s1 = np.datetime64(end or '2100-01', 'D')
    out = []
    for s in stations:
        if box is not None:
            lon = (s.lon - box[2]) % 360 + box[2]
            if not (box[0] <= s.lat <= box[1] and box[2] <= lon <= box[3]):
                continue
        if inside is not None and not bool(inside(np.array([s.lat]), np.array([s.lon]))[0]):
            continue
        a, b = max(np.datetime64(s.start, 'D'), s0), min(np.datetime64(s.end, 'D'), s1)
        if (b - a).astype(int) / 365.25 < min_years or s.n_solutions < min_solutions:
            continue
        out.append(s.name)
    return out


def _fetch(url, path, timeout=60, check=None):
    req = urllib.request.Request(url, headers={'User-Agent': 'PyGWater'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    if data[:2] == b'\x1f\x8b':                       # gzip-compressed
        import gzip
        data = gzip.decompress(data)
    if check is not None and not check(data):
        raise ValueError('not a valid file: %s' % url)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_bytes(data)
    return len(data)


def _is_tenv3(data):
    """a .tenv3 file starts with the header line 'site YYMMMDD ...' and has data lines"""
    head = data[:200].decode('ascii', 'replace')
    return head.startswith('site') and data.count(b'\n') > 1


def download_catalogue(out_dir):
    """DataHoldings.txt and steps.txt into out_dir"""
    out_dir = Path(out_dir)
    for url in (HOLDINGS_URL, STEPS_URL):
        _fetch(url, out_dir / url.rsplit('/', 1)[1])
    return out_dir


def download_series(sites, out_dir, frame='IGS20', overwrite=False, pause=0.2):
    """
    Daily .tenv3 files of the stations, saved as out_dir/SITE.tenv3.
    :return: dict {site: 'ok' | 'present' | error message}
    """
    out_dir = Path(out_dir)
    status = {}
    for site in sites:
        path = out_dir / ('%s.tenv3' % site)
        if path.exists() and not overwrite:
            status[site] = 'present'
            continue
        err = None
        for pattern in TENV3_URLS[frame]:
            try:
                _fetch(pattern.format(site=site), path, check=_is_tenv3)
                status[site], err = 'ok', None
                break
            except Exception as e:  # try the next URL pattern
                err = str(e)
        if err is not None:
            status[site] = err
        time.sleep(pause)
    return status
