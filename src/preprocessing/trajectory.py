"""
Trajectory model of a GNSS position component (Bevis & Brown, 2014):

    y(t) = bias + rate (t - t0) + sum_k step_k H(t - t_k) + sum_q post_q log(1 + (t - t_q) / tau) H(t - t_q)
           + annual and semi-annual sine/cosine + residual

fitted by weighted least squares with iterative outlier rejection (|residual| > outlier x 1.4826 MAD).
For loading studies the trend, steps and post-seismic terms are removed and the seasonal terms and the residual are
kept (they hold the loading signal).
"""
import numpy as np


def trajectory_design(t, steps=(), annual=True, semiannual=True, t0=None, postseismic=(), tau=10 / 365.25):
    """
    :param steps: decimal years of offsets (equipment changes, earthquakes)
    :param postseismic: decimal years of earthquakes that get a logarithmic post-seismic term
    :param tau: relaxation time of the logarithmic term [years] (default 10 days)
    :return: G (n, p), list of column names
    """
    t = np.asarray(t, dtype=float)
    t0 = np.mean(t) if t0 is None else t0
    cols, names = [np.ones_like(t), t - t0], ['bias', 'rate']
    for k, ts in enumerate(steps):
        if t.min() < ts < t.max():
            cols.append((t >= ts).astype(float)); names.append('step%d' % k)
    for k, tq in enumerate(postseismic):
        if t.min() < tq < t.max():
            cols.append(np.where(t >= tq, np.log1p(np.maximum(t - tq, 0) / tau), 0.0)); names.append('post%d' % k)
    for per, on in ((1.0, annual), (0.5, semiannual)):
        if on:
            w = 2 * np.pi / per
            cols += [np.sin(w * (t - t0)), np.cos(w * (t - t0))]
            names += ['sin%g' % per, 'cos%g' % per]
    return np.column_stack(cols), names


def fit_trajectory(t, y, sig=None, steps=(), annual=True, semiannual=True, outlier=5.0, max_iter=5, postseismic=()):
    """
    Weighted least-squares trajectory fit with iterative outlier rejection.
    :return: params (dict name -> value), fitted model (n,), keep mask (n,)
    """
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    w = np.ones_like(y) if sig is None else 1.0 / np.maximum(np.asarray(sig, dtype=float), 1e-6)
    G, names = trajectory_design(t, steps, annual, semiannual, postseismic=postseismic)
    keep = np.isfinite(y)
    x = np.zeros(G.shape[1])
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


def seasonal_amplitude(params, period=1.0):
    """amplitude of the annual (period=1) or semi-annual (0.5) term of a fitted model [same unit as y]"""
    return float(np.hypot(params.get('sin%g' % period, 0.0), params.get('cos%g' % period, 0.0)))


def detect_jumps(t, r, window=7, n_sigma=5.0, min_size=3.5e-3, known=(), exclusion=7):
    """
    Offsets missing from the steps database, after Schumacher et al. (2018, GJI 214, Sect. 3.1.4): an offset changes
    the mean of the residuals between the `window` days before and after it.
        D(k) = mean(r over the next `window` days) - mean(r over the previous `window` days)
    Days where |D| exceeds max(n_sigma * sigma_D, min_size) form groups; within a group the offset is placed at the
    largest |D|. sigma_D is a robust (MAD) standard deviation of D itself, i.e. of the noise of a 7-day mean
    difference (Schumacher et al. set 3 sigma of the series with a 3.5 mm floor; the daily scatter is too large a
    yardstick for 7-day means, so the threshold is taken from D, with a higher default factor).
    :param t: decimal years (sorted); r: residuals of a trajectory fit that already includes the known offsets [m]
    :param known: decimal years of offsets already modelled; detections within `exclusion` days of them are ignored
    :return: list of decimal years of the detected offsets
    """
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    if t.size < 4 * window:
        return []
    day = np.round((t - t[0]) * 365.25).astype(int)
    n = day.max() + 1
    s, c = np.zeros(n), np.zeros(n)
    np.add.at(s, day, r)
    np.add.at(c, day, 1)
    cs, cc = np.concatenate([[0], np.cumsum(s)]), np.concatenate([[0], np.cumsum(c)])
    k = np.arange(window, n - window + 1)
    nb, na = cc[k] - cc[k - window], cc[k + window] - cc[k]
    ok = (nb >= window / 2) & (na >= window / 2)
    D = np.full(k.size, 0.0)
    D[ok] = (cs[k + window] - cs[k])[ok] / na[ok] - (cs[k] - cs[k - window])[ok] / nb[ok]
    Dv = D[ok]
    sigma = 1.4826 * np.median(np.abs(Dv - np.median(Dv))) if Dv.size else np.inf
    thr = max(n_sigma * sigma, min_size)
    hit = np.abs(D) > thr
    jumps, i = [], 0
    known_days = [(kt - t[0]) * 365.25 for kt in known]
    while i < k.size:
        if not hit[i]:
            i += 1
            continue
        j = i
        while j + 1 < k.size and hit[j + 1]:
            j += 1
        best = i + int(np.argmax(np.abs(D[i:j + 1])))
        jd = k[best]                                       # the offset happens between day jd-1 and day jd
        if all(abs(jd - kd) > exclusion for kd in known_days):
            jumps.append(float(t[0] + (jd - 0.5) / 365.25))
        i = j + 1
    return jumps


def clean_segments(t, steps=(), min_segment=0.5, gap=90 / 365.25, min_between=14 / 365.25):
    """
    Data-selection rules of Schumacher et al. (2018, Sect. 3.1.2):
      - short isolated pieces of data (shorter than `min_segment` years, separated from the rest by gaps longer than
        `gap`) are dropped, unless they hold most of the data;
      - data between two offsets less than `min_between` apart are dropped and only one offset is kept.
    :return: keep mask (n,), list of the offsets kept
    """
    t = np.asarray(t, dtype=float)
    keep = np.ones(t.size, bool)
    cut = np.nonzero(np.diff(t) > gap)[0]
    edges = np.concatenate([[0], cut + 1, [t.size]])
    segs = [(a, b) for a, b in zip(edges[:-1], edges[1:])]
    longest = max(segs, key=lambda ab: t[ab[1] - 1] - t[ab[0]])
    for a, b in segs:
        if (a, b) != longest and t[b - 1] - t[a] < min_segment:
            keep[a:b] = False
    st = sorted(steps)
    kept = []
    for s in st:
        if kept and s - kept[-1] < min_between:
            keep &= ~((t >= kept[-1]) & (t < s))
            continue
        kept.append(s)
    return keep, kept
