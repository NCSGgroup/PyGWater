"""
Download the NGL station list, the steps file and the daily .tenv3 series of the stations in a region.
Run on a computer with internet access:

    python studies/GNSS_GRACE/download.py

Files go to <data root>/input/gnss/NGL/ (DataHoldings.txt, steps.txt, tenv3/SITE.tenv3); files already present are
skipped, so the script can be re-run to complete a download. The region and period are those of config.py.
"""
from pygwater.config import FileTool
from pygwater.preprocessing import ngl

from config import CONFIG

if __name__ == '__main__':
    out = FileTool.get_input_dir('gnss/NGL')
    if not (out / 'DataHoldings.txt').exists() or CONFIG.get('refresh_catalogue', False):
        ngl.download_catalogue(out)
    stations = ngl.read_holdings(out / 'DataHoldings.txt')
    sites = ngl.select_stations(stations, box=CONFIG['box'], start=CONFIG['start'], end=CONFIG['end'],
                                min_years=CONFIG['min_years'])
    print('%d stations in %s with at least %g years between %s and %s' % (len(sites), CONFIG['box'],
                                                                         CONFIG['min_years'], CONFIG['start'],
                                                                         CONFIG['end']))
    status = ngl.download_series(sites, out / 'tenv3', frame=CONFIG['frame'])
    failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
    print('downloaded %d, already present %d, failed %d' % (sum(v == 'ok' for v in status.values()),
                                                           sum(v == 'present' for v in status.values()), len(failed)))
    for k, v in list(failed.items())[:10]:
        print('  %s: %s' % (k, v))
