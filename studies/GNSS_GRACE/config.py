"""Settings of the GNSS-GRACE study; edit here, then run download.py (once, with internet) and run.py."""
CONFIG = dict(
    name='PNW',                                 # output folder: <data root>/output/GNSS_GRACE/<name>/
    box=(41.0, 50.0, -125.0, -116.0),           # stations: lat_min, lat_max, lon_min, lon_max [deg]
    start='2006-01', end='2022-12',             # period
    min_years=4.0,                              # at least 4 years of daily data (Schumacher et al., 2018)
    frame='IGS20',                              # NGL archive
    reference=('2006-01', '2015-12'),           # common reference period of GNSS and GRACE anomalies
    components=('U',),                          # ('U',) or ('U', 'N', 'E')
    resolution=0.5,                             # cell size of the inversion [deg]
    buffer=2.0,                                 # cells extend this far beyond the station box [deg]
    kernel='disk',                              # 'point_mass', 'disk' or 'cluster_disk'
    lmax=13000,                                 # Love numbers (Wang et al., 2012), CF frame
    alpha_method='gcv',                         # 'gcv', 'discrepancy', 'lcurve'
    nonnegative=False,                          # anomalies can be negative
    spatial_filter_km=500,                      # neighbour test radius (Schumacher et al., 2018)
    grace_file='mascon/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',   # under <data root>/input/
)
