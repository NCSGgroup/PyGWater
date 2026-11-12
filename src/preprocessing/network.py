"""
A GNSS network from raw daily series to monthly displacement anomalies ready for inversion:

    net = GNSSNetwork.from_files(files)                              # read_tenv3 for each file
    net.clean(steps=read_ngl_steps(path, as_events=True))            # offsets, trajectory, outliers, quality control
    net.monthly()                                                    # calendar-month means and realistic sigmas
    net.spatial_filter(radius_km=500)                                # drop stations dominated by local effects
    ds = net.to_dataset(reference=('2004-01', '2009-12'))            # GNSSDataset: (month x station) arrays
    lat, lon, d, sigma, comps = ds.epoch('2020-01')                  # one month, stacked U, N, E

clean() removes trend, offsets and (optionally) post-seismic terms and keeps the seasonal signal and the residual,
which hold the loading signal. It also computes quality-control figures for every station and rejects stations that
fail the thresholds; net.qc_table() lists them, net.rejected says why the others were dropped.

Reference frame: NGL series are in IGS14/IGS20, whose origin is the centre of mass for the long-term trend but close
to the centre of figure (CF) for seasonal loading. Compare with load Love numbers converted to the CF frame
(LoveNumber.convert(Frame.CF)). Common-mode error: removing a network-wide common mode also removes part of the
regional loading signal, so common_mode() is provided but not applied.
"""
import csv
from dataclasses import dataclass

import numpy as np

from pygwater.preprocessing.series import COMPONENTS, Step, read_tenv3, _decyear_to_month
from pygwater.preprocessing.trajectory import (trajectory_design, fit_trajectory, seasonal_amplitude, clean_segments,
                                               detect_jumps as find_jumps)


def _split_steps(entries, postseismic_magnitude):
    """{site: [float | Step]} entry -> (step times, post-seismic times)"""
    steps, post = [], []
    for s in entries or ():
        if isinstance(s, Step):
            steps.append(s.t)
            if s.kind == 'earthquake' and postseismic_magnitude is not None:
                try:
                    mag = float(s.info.split()[1])
                except (IndexError, ValueError):
                    mag = 0.0
                if mag >= postseismic_magnitude:
                    post.append(s.t)
        else:
            steps.append(float(s))
    return sorted(set(steps)), sorted(set(post))


class GNSSNetwork:

    def __init__(self, stations):
        self.stations = list(stations)
        self.rejected = {}          # name -> reason

    @classmethod
    def from_files(cls, files):
        """one station per .tenv3 file"""
        return cls([read_tenv3(f) for f in files])

    def names(self):
        return [s.name for s in self.stations]

    def select(self, box=None, inside=None):
        """
        Keep the stations in a region.
        :param box: (lat_min, lat_max, lon_min, lon_max) [degree]
        :param inside: or a function (lat, lon) -> bool, e.g. from a shapefile (shapely.contains_xy)
        """
        keep = []
        for s in self.stations:
            ok = True
            if box is not None:
                lon = (s.lon - box[2]) % 360 + box[2]
                ok = box[0] <= s.lat <= box[1] and box[2] <= lon <= box[3]
            if ok and inside is not None:
                ok = bool(inside(np.array([s.lat]), np.array([s.lon]))[0])
            (keep.append(s) if ok else self.rejected.__setitem__(s.name, 'outside the region'))
        self.stations = keep
        return self

    # --------------------------------------------------------------------------------------------------------------
    def clean(self, steps=None, remove='trend+steps', outlier=5.0, annual=True, semiannual=True, min_years=4.0,
              min_completeness=0.5, max_sigma=0.02, max_residual=0.01, postseismic_magnitude=None, period=None,
              corrections=None, detect_jumps=True, jump_options=None, segment_options=None):
        """
        From the daily series as read to the loading signal, station by station:
          1. keep the period of interest; drop days with a large formal vertical sigma;
          2. subtract external loading models (e.g. non-tidal atmosphere and ocean), if given;
          3. drop short isolated pieces of data and the data between offsets less than two weeks apart, and require
             at least `min_years` x 365 daily solutions (Schumacher et al., 2018, Sect. 3.1.2);
          4. fit the trajectory model (trend, offsets, optional post-seismic terms, annual and semi-annual terms) with
             outlier rejection; optionally detect offsets missing from the steps database (Sect. 3.1.4) and refit;
          5. remove the chosen terms and record quality-control figures; reject stations that fail the thresholds.
        Re-running clean() (e.g. with other options) starts again from the series as read.
        :param steps: {site: [decimal years]} or {site: [Step]} (read_ngl_steps(..., as_events=True)); None = none
        :param remove: 'trend+steps' (keep seasonal + residual, the loading signal; default),
                       'trend+steps+seasonal' (keep the residual only) or 'none'
        :param min_years: minimum number of daily solutions, in years (count, not span; Schumacher et al.: 4)
        :param min_completeness: minimum fraction of days with a solution over the span of the series
        :param max_sigma: daily solutions with a formal vertical sigma above this [m] are dropped
        :param max_residual: reject stations whose vertical residual RMS (after the model) exceeds this [m]
        :param postseismic_magnitude: earthquakes of at least this magnitude (needs Step events) also get a
                                      logarithmic post-seismic term; None = offsets only
        :param period: (t_start, t_end) in decimal years to keep, e.g. (2003.0, 2023.0); None = all
        :param corrections: {site: (t, {'U': ..., 'N': ..., 'E': ...})} displacements of loads that are not wanted in
                            the signal (e.g. non-tidal atmosphere and ocean, in the same frame) [m]; interpolated in
                            time and subtracted
        :param detect_jumps: look for offsets missing from the steps database (trajectory.detect_jumps, vertical)
        :param jump_options, segment_options: keyword arguments of trajectory.detect_jumps / clean_segments
        """
        keep_st = []
        for s in self.stations:
            if s.raw is None:
                s.raw = (s.t.copy(), {c: v.copy() for c, v in s.xyz.items()}, {c: v.copy() for c, v in s.sig.items()})
            t, xyz, sig = s.raw
            sel = np.ones(t.size, bool) if period is None else (t >= period[0]) & (t <= period[1])
            sel &= sig['U'] <= max_sigma
            s.t = t[sel].copy()
            s.xyz = {c: v[sel].copy() for c, v in xyz.items()}
            s.sig = {c: v[sel].copy() for c, v in sig.items()}
            s.model, s.qc = {}, {}
            s.monthly_t = s.monthly_xyz = s.monthly_sig = None
            if s.t.size < 30:
                self.rejected[s.name] = 'fewer than 30 daily solutions'
                continue
            corrected = corrections is not None and s.name in corrections
            if corrected:
                tc, dc = corrections[s.name]
                for c in COMPONENTS:
                    if c in dc:
                        s.xyz[c] = s.xyz[c] - np.interp(s.t, tc, dc[c])
            st, post = _split_steps((steps or {}).get(s.name, ()), postseismic_magnitude)
            seg, st = clean_segments(s.t, st, **(segment_options or {}))
            n_seg = int((~seg).sum())
            self._subset(s, seg)
            if s.t.size < min_years * 365:
                self.rejected[s.name] = 'fewer than %g years of daily solutions' % min_years
                continue
            completeness = s.t.size / ((s.t.max() - s.t.min()) * 365.25 + 1)
            if completeness < min_completeness:
                self.rejected[s.name] = 'completeness %.2f below %.2f' % (completeness, min_completeness)
                continue
            found = []
            if detect_jumps:
                p, fit, k = fit_trajectory(s.t, s.xyz['U'], s.sig['U'], st, annual, semiannual, outlier,
                                           postseismic=post)
                found = find_jumps(s.t[k], (s.xyz['U'] - fit)[k], known=st, **(jump_options or {}))
                st = sorted(st + found)
            ok = np.ones_like(s.t, dtype=bool)
            resid = {}
            for c in COMPONENTS:
                p, fit, k = fit_trajectory(s.t, s.xyz[c], s.sig[c], st, annual, semiannual, outlier,
                                           postseismic=post)
                s.model[c] = p
                ok &= k
                resid[c] = s.xyz[c] - fit
                drop = ['bias', 'rate'] + [n for n in p if n.startswith(('step', 'post'))]
                if remove == 'trend+steps+seasonal':
                    drop += [n for n in p if n.startswith(('sin', 'cos'))]
                if remove != 'none':
                    G, names = trajectory_design(s.t, st, annual, semiannual, postseismic=post)
                    s.xyz[c] = s.xyz[c] - G[:, [names.index(n) for n in drop]] @ np.array([p[n] for n in drop])
            rms = {c: float(np.sqrt(np.mean(resid[c][ok] ** 2))) for c in COMPONENTS}
            s.qc = dict(n_days=int(ok.sum()), span_years=float(s.t.max() - s.t.min()),
                        completeness=float(completeness), outlier_fraction=float(1 - ok.mean()),
                        dropped_segment_days=n_seg, n_steps=len(st), n_detected_jumps=len(found),
                        n_postseismic=len(post), loading_corrected=bool(corrected),
                        median_sigma_U=float(np.median(s.sig['U'])), annual_U=seasonal_amplitude(s.model['U']),
                        rms_U=rms['U'], rms_N=rms['N'], rms_E=rms['E'])
            s.detected_jumps = found
            if rms['U'] > max_residual:
                self.rejected[s.name] = 'vertical residual RMS %.1f mm above %.1f mm' % (1e3 * rms['U'],
                                                                                         1e3 * max_residual)
                continue
            self._subset(s, ok)
            keep_st.append(s)
        self.stations = keep_st
        return self

    @staticmethod
    def _subset(s, mask):
        s.t = s.t[mask]
        for c in COMPONENTS:
            s.xyz[c], s.sig[c] = s.xyz[c][mask], s.sig[c][mask]

    def monthly(self, min_days=10, sigma='rms'):
        """
        Calendar-month means. Months with fewer than min_days daily solutions are dropped.
        :param sigma: uncertainty of a monthly value:
                      'rms' : RMS of the daily values about the monthly mean (Khan et al., 2016; Schumacher et al.,
                              2018), realistic for the time-correlated noise of GNSS series (default)
                      'sem' : standard error of the mean, max(std, mean formal sigma) / sqrt(n); assumes white noise
                              and is optimistic by a factor of 3-5
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
                    if sigma == 'rms':
                        sg = max(np.sqrt(np.mean((v - v.mean()) ** 2)), s.sig[c][sel].mean() / np.sqrt(n))
                    else:
                        sg = max(v.std(ddof=1), s.sig[c][sel].mean()) / np.sqrt(n)
                    s.monthly_xyz[c].append(v.mean())
                    s.monthly_sig[c].append(sg)
            for c in COMPONENTS:
                s.monthly_xyz[c], s.monthly_sig[c] = np.array(s.monthly_xyz[c]), np.array(s.monthly_sig[c])
        return self

    def spatial_filter(self, radius_km=500.0, component='U', min_neighbours=3, min_corr=0.0, max_misfit=3.0,
                       min_common=24):
        """
        Remove stations dominated by local effects (e.g. poroelastic aquifer response, pumping, local
        instabilities), adapting the 500 km median filter of Schumacher et al. (2018, Sect. 3.2.3) from trends to
        monthly series: each station is compared with the month-by-month median of its neighbours within
        `radius_km`. A station is rejected if
          - its correlation with the neighbour median is below `min_corr` (default 0: it moves against the region,
            as poroelastic sites do under a load), or
          - its RMS difference from the neighbour median exceeds `max_misfit` times the median RMS difference of the
            network.
        Stations with fewer than `min_neighbours` neighbours (or `min_common` common months) are kept and flagged.
        Run after monthly(). :return: self
        """
        from pygwater.utils.geomath import GeoMathKit
        st = self.stations
        if len(st) < min_neighbours + 1:
            return self
        lat, lon = np.array([s.lat for s in st]), np.array([s.lon for s in st])
        months = self.months()
        idx = {k: i for i, k in enumerate(months)}
        Y = np.full((len(st), len(months)), np.nan)
        for j, s in enumerate(st):
            v = s.monthly_xyz[component]
            Y[j, [idx[k] for k in s.monthly_t]] = v - np.mean(v)
        R = 6371.0
        corr, misfit = np.full(len(st), np.nan), np.full(len(st), np.nan)
        for j, s in enumerate(st):
            dist = R * np.deg2rad(GeoMathKit.angular_distance(lat[j], lon[j], lat, lon))
            nb = np.nonzero((dist <= radius_km) & (np.arange(len(st)) != j))[0]
            s.qc['n_neighbours'] = int(nb.size)
            if nb.size < min_neighbours:
                s.qc['spatial_filter'] = 'too few neighbours'
                continue
            med = np.nanmedian(Y[nb], axis=0)
            cnt = np.sum(np.isfinite(Y[nb]), axis=0)
            use = np.isfinite(Y[j]) & (cnt >= min_neighbours)
            if use.sum() < min_common:
                s.qc['spatial_filter'] = 'too few common months'
                continue
            a, b = Y[j, use] - Y[j, use].mean(), med[use] - med[use].mean()
            corr[j] = float(np.corrcoef(a, b)[0, 1])
            misfit[j] = float(np.sqrt(np.mean((a - b) ** 2)))
            s.qc['neighbour_corr'], s.qc['neighbour_rms'] = corr[j], misfit[j]
        typical = np.nanmedian(misfit) if np.isfinite(misfit).any() else np.inf
        keep = []
        for j, s in enumerate(st):
            if np.isfinite(corr[j]) and corr[j] < min_corr:
                self.rejected[s.name] = 'moves against its neighbours (r = %.2f)' % corr[j]
            elif np.isfinite(misfit[j]) and misfit[j] > max_misfit * typical:
                self.rejected[s.name] = 'differs from its neighbours (RMS %.1f mm, network median %.1f mm)' % (
                    1e3 * misfit[j], 1e3 * typical)
            else:
                s.qc.setdefault('spatial_filter', 'passed')
                keep.append(s)
        self.stations = keep
        return self

    def months(self):
        return sorted({k for s in self.stations for k in (s.monthly_t or [])})

    def common_mode(self, month_key=None):
        """
        Network common mode per month and component (median over stations) of the monthly series. Not removed
        automatically: it also contains the regional loading signal.
        :return: {month: {comp: value}}
        """
        out = {}
        for m in self.months() if month_key is None else [month_key]:
            out[m] = {}
            for c in COMPONENTS:
                v = [s.monthly_xyz[c][s.monthly_t.index(m)] for s in self.stations if m in s.monthly_t]
                out[m][c] = float(np.median(v)) if v else np.nan
        return out

    def epoch(self, month, components=COMPONENTS, demean=True, reference=None, min_months=12):
        """
        Data vector for one month (see GNSSDataset.epoch); kept for convenience, equivalent to
        self.to_dataset(components, reference, min_months).epoch(month) when demean is True.
        """
        if not demean:
            st = [s for s in self.stations if s.monthly_t and month in s.monthly_t]
            d = np.array([s.monthly_xyz[c][s.monthly_t.index(month)] for c in components for s in st])
            sg = np.array([s.monthly_sig[c][s.monthly_t.index(month)] for c in components for s in st])
            return (np.array([s.lat for s in st]), np.array([s.lon for s in st]), d, sg, tuple(components))
        return self.to_dataset(components, reference, min_months).epoch(month)

    # --------------------------------------------------------------------------------------------------------------
    def to_dataset(self, components=COMPONENTS, reference=None, min_months=12):
        """
        Monthly anomalies of all stations as arrays, each station relative to its mean over a COMMON reference
        period, so that every station (and the inversion) refers to the same period, as GRACE anomalies do.
        :param reference: ('YYYY-MM', 'YYYY-MM'), first and last month of the reference period (inclusive), e.g.
                          ('2004-01', '2009-12') for the GRACE baseline. Default: the months in which every station
                          has data
        :param min_months: stations with fewer monthly values in the reference period are left out
        :return: GNSSDataset
        """
        if any(s.monthly_t is None for s in self.stations):
            self.monthly()
        if reference is None:
            ref = set.intersection(*[set(s.monthly_t) for s in self.stations]) if self.stations else set()
            if not ref:
                raise ValueError('the stations share no month; give a reference period, e.g. '
                                 "reference=('2004-01', '2009-12')")
        else:
            ref = {k for k in self.months() if reference[0] <= k <= reference[1]}
        st, means = [], []
        for s in self.stations:
            sel = [i for i, k in enumerate(s.monthly_t) if k in ref]
            if len(sel) >= min(min_months, len(ref)) and sel:
                st.append(s)
                means.append({c: s.monthly_xyz[c][sel].mean() for c in components})
            else:
                self.rejected.setdefault(s.name, 'fewer than %d months in the reference period' % min_months)
        months = sorted({k for s in st for k in s.monthly_t})
        index = {k: i for i, k in enumerate(months)}
        val = {c: np.full((len(months), len(st)), np.nan) for c in components}
        sig = {c: np.full((len(months), len(st)), np.nan) for c in components}
        for j, (s, mu) in enumerate(zip(st, means)):
            rows = [index[k] for k in s.monthly_t]
            for c in components:
                val[c][rows, j] = s.monthly_xyz[c] - mu[c]
                sig[c][rows, j] = s.monthly_sig[c]
        return GNSSDataset(names=[s.name for s in st], lat=np.array([s.lat for s in st]),
                           lon=np.array([s.lon for s in st]), months=months, values=val, sigma=sig,
                           components=tuple(components), reference=(min(ref), max(ref)))

    # --------------------------------------------------------------------------------------------------------------
    def qc_table(self):
        """one row per station (kept and rejected): name, lat, lon, status and the quality-control figures"""
        rows = [dict(name=s.name, lat=s.lat, lon=s.lon, status='kept', **s.qc) for s in self.stations]
        rows += [dict(name=n, status='rejected: ' + r) for n, r in self.rejected.items()]
        return rows

    def save_qc(self, path):
        rows = self.qc_table()
        cols = list(dict.fromkeys(k for r in rows for k in r))
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)


@dataclass
class GNSSDataset:
    """Monthly displacement anomalies of a network: values[c][month, station], NaN where a station has no data."""
    names: list
    lat: np.ndarray
    lon: np.ndarray
    months: list                    # 'YYYY-MM'
    values: dict                    # {component: (n_month, n_station)} [m]
    sigma: dict                     # same shape [m]
    components: tuple
    reference: tuple                # first and last month of the reference period

    def epoch(self, month, components=None):
        """
        One month, stacked component by component (as ForwardModel.design_matrix stacks rows); stations without
        data in that month are left out.
        :return: lat, lon, d (n_comp * N,), sigma (n_comp * N,), components
        """
        comps = tuple(components or self.components)
        i = self.months.index(month)
        ok = np.all([np.isfinite(self.values[c][i]) for c in comps], axis=0)
        d = np.concatenate([self.values[c][i, ok] for c in comps])
        sg = np.concatenate([self.sigma[c][i, ok] for c in comps])
        return self.lat[ok], self.lon[ok], d, sg, comps

    def station_mask(self, month, components=None):
        comps = tuple(components or self.components)
        i = self.months.index(month)
        return np.all([np.isfinite(self.values[c][i]) for c in comps], axis=0)

    def save(self, path):
        np.savez(path, names=np.array(self.names), lat=self.lat, lon=self.lon, months=np.array(self.months),
                 components=np.array(self.components), reference=np.array(self.reference),
                 **{'value_' + c: self.values[c] for c in self.components},
                 **{'sigma_' + c: self.sigma[c] for c in self.components})

    @classmethod
    def load(cls, path):
        f = np.load(path)
        comps = tuple(str(c) for c in f['components'])
        return cls(names=[str(n) for n in f['names']], lat=f['lat'], lon=f['lon'],
                   months=[str(m) for m in f['months']], values={c: f['value_' + c] for c in comps},
                   sigma={c: f['sigma_' + c] for c in comps}, components=comps,
                   reference=tuple(str(r) for r in f['reference']))
