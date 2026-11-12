"""
GNSS preprocessing: from raw daily station series to monthly displacement anomalies ready for inversion.

    ngl          station list, station selection and download from the Nevada Geodetic Laboratory
    series       StationSeries, read_tenv3, read_ngl_steps (equipment changes and earthquakes)
    trajectory   trajectory model (trend, offsets, post-seismic, seasonal) with outlier rejection
    network      GNSSNetwork (clean, quality control, monthly means, common reference period) and GNSSDataset
"""
from pygwater.preprocessing.series import COMPONENTS, StationSeries, Step, read_tenv3, read_ngl_steps
from pygwater.preprocessing.trajectory import trajectory_design, fit_trajectory
from pygwater.preprocessing.network import GNSSNetwork, GNSSDataset
