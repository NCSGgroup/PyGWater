"""
@Company: CGE-HUST, Wuhan, China
@Version: 1.0
@Author: Yang Fan
@Contact : yfan_cge@hust.edu.cn
@Modify Time: 2020/6/4 12:02
@Description:
"""
import calendar
import datetime
import gzip
import math

import numpy as np
# import warnings
#
# warnings.filterwarnings("error")


class GeoMathKit:

    def __init__(self):
        pass

    @staticmethod
    def angular_distance(point1_lat, point1_lon, point2_lat, point2_lon):
        """
        Calculate the angular distance (degree) between two points
        :param point1_lat: [degree]
        :param point1_lon: [degree]
        :param point2_lat: [degree]
        :param point2_lon: [degree]
        :return: [degree]
        """
        beta1 = np.deg2rad(point1_lat)
        alpha1 = np.deg2rad(point1_lon)
        beta2 = np.deg2rad(point2_lat)
        alpha2 = np.deg2rad(point2_lon)

        m = np.cos(beta1) * np.cos(beta2) * np.cos(alpha1 - alpha2) + np.sin(beta1) * np.sin(beta2)
        m = np.clip(m, -1, 1)  # works for scalars and arrays
        theta = np.arccos(m)
        return np.rad2deg(theta)

    @staticmethod
    def azimuth(point1_lat, point1_lon, point2_lat, point2_lon):
        """
        Azimuth at point 1 of the great circle towards point 2, clockwise from north.
        :param point1_lat, point1_lon, point2_lat, point2_lon: [degree], scalars or broadcastable arrays
        :return: [rad]
        """
        p1, p2 = np.deg2rad(point1_lat), np.deg2rad(point2_lat)
        dl = np.deg2rad(np.asarray(point2_lon) - np.asarray(point1_lon))
        return np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))

    @staticmethod
    def horizontal_to_north_east(v, load_lat, load_lon, points_lat, points_lon):
        """
        Project the great-circle horizontal displacement v (positive AWAY from the load, Farrell 1972), evaluated at
        the points, onto local north and east.
            u_N = -v cos(az),  u_E = -v sin(az),   az = azimuth at the point towards the load
        :return: (u_N, u_E)
        """
        az = GeoMathKit.azimuth(points_lat, points_lon, load_lat, load_lon)
        return -v * np.cos(az), -v * np.sin(az)

    @staticmethod
    def kernel_at_points(obj, load_lat, load_lon, points_lat, points_lon, variable, **kwargs):
        """
        Evaluate a Green's-function object (LGF_truncation, LGF_Kummer, DiskLoad, ...) for a load centred at
        (load_lat, load_lon) at the points, for any variable of enums.Displacement:
        Vertical, Geoheight, Horizontal (great-circle, away from load), North, East, Gravity.
        :param kwargs: passed to getGravity (e.g. newtonian='above')
        """
        from pygwater.utils.enums import Displacement
        theta = GeoMathKit.angular_distance(load_lat, load_lon, points_lat, points_lon)
        obj.configure(residence=theta)
        if variable == Displacement.Vertical:
            return obj.getVertical()
        elif variable == Displacement.Geoheight:
            return obj.getGeoidheight()
        elif variable == Displacement.Horizontal:
            return obj.getHorizental()
        elif variable == Displacement.Gravity:
            return obj.getGravity(**kwargs)
        elif variable in (Displacement.North, Displacement.East):
            un, ue = GeoMathKit.horizontal_to_north_east(obj.getHorizental(), load_lat, load_lon,
                                                         points_lat, points_lon)
            return un if variable == Displacement.North else ue
        else:
            raise NotImplementedError(variable)

    @staticmethod
    def haversine(lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        :param lon1:  point 1
        :param lat1:
        :param lon2: point 2
        :param lat2:
        :return:
        """

        '''Degree to radian'''
        lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

        ''' haversine'''
        dlon = lon2 - lon1
        dlat = lat2 - lat1

        a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
        c = 2 * np.arcsin(np.sqrt(a))
        r = 6371  # average radius of the earth, unit [km]

        return c * r * 1000

    @staticmethod
    def getCoLatLoninRad(lat, lon):
        '''
        :param lat: geophysical coordinate in degree
        :param lon: geophysical coordinate in degree
        :return: Co-latitude and longitude in rad
        '''

        theta = (90. - lat) / 180. * np.pi
        phi = lon / 180. * np.pi

        return theta, phi

    @staticmethod
    def CS_1dTo2d(CS: np.ndarray):
        """

        C00 C10 C20 =>     C00
                           C10 C20

        :param CS: one-dimension array
        :return: two dimension array
        """

        def index(N):
            n = (np.round(np.sqrt(2 * N))).astype(np.int) - 1
            m = N - (n * (n + 1) / 2).astype(np.int) - 1
            return n, m

        CS_index = np.arange(len(CS)) + 1
        n, m = index(CS_index)

        dim = index(len(CS))[0] + 1
        CS2d = np.zeros((dim, dim))
        CS2d[n, m] = CS

        return CS2d

    @staticmethod
    def CS_2dTo1d(CS: np.ndarray):
        """
        Transform the CS in 2-dimensional matrix to 1-dimemsion vectors
        example:
        00
        10 11
        20 21 22
        30 31 32 33           =>           00 10 11 20 21 22 30 31 32 33 ....

        :param CS:
        :return:
        """
        shape = np.shape(CS)
        assert len(shape) == 2
        index = np.nonzero(np.tril(np.ones(shape)))

        return CS[index]

    @staticmethod
    def dayListByMonth(begin, end):
        """
        get the date of every day between the given 'begin' month and 'end' month

        :param begin: year, month
        :param end: year,month
        :return:
        """

        daylist = []
        begin_date = datetime.date(begin[0], begin[1], 1)
        end_date = datetime.date(end[0], end[1], calendar.monthrange(end[0], end[1])[1])

        while begin_date <= end_date:
            date_str = begin_date
            daylist.append(date_str)
            begin_date += datetime.timedelta(days=1)

        return daylist

    @staticmethod
    def monthListByMonth(begin, end):
        """
        get the date of every day between the given 'begin' month and 'end' month

        :param begin: '2008-01'
        :param end: '2009-07'
        :return:
        """

        begin = [int(x) for x in begin.split('-')]
        end = [int(x) for x in end.split('-')]
        daylist = []
        begin_date = datetime.date(begin[0], begin[1], 1)
        end_date = datetime.date(end[0], end[1], calendar.monthrange(end[0], end[1])[1])

        while begin_date <= end_date:
            date_str = begin_date
            daylist.append(date_str)
            begin_date += datetime.timedelta(days=1)

        monthlist = []
        for day in daylist:
            if day.day == 1:
                monthlist.append(day)

        return monthlist

    @staticmethod
    def dayListByDay(begin, end):
        """
        get the date of every day between the given 'begin' day and 'end' day

        :param begin: year, month, day. '2009-01-01'
        :param end: year,month,day. '2010-01-01'
        :return:
        """

        daylist = []
        begin_date = datetime.datetime.strptime(begin, "%Y-%m-%d")
        end_date = datetime.datetime.strptime(end, "%Y-%m-%d")

        while begin_date <= end_date:
            date_str = begin_date
            daylist.append(date_str)
            begin_date += datetime.timedelta(days=1)

        return daylist

    @staticmethod
    def un_gz(file_name):

        # aquire the filename and remove the postfix
        f_name = file_name.replace(".gz", "")
        # start uncompress
        g_file = gzip.GzipFile(file_name)
        # read uncompressed files and write down a copy without postfix
        open(f_name, "wb+").write(g_file.read())
        g_file.close()

    @staticmethod
    def year_fraction(date):
        start = datetime.date(date.year, 1, 1).toordinal()
        year_length = datetime.date(date.year + 1, 1, 1).toordinal() - start
        return date.year + float(date.toordinal() - start) / year_length

    @staticmethod
    def global_equal_distance_grid(grid_size=1):
        """
        :param grid_size: unit [degree]
        :return: [degree]
        """
        lon = np.arange(-180 + grid_size / 2, 180, step=grid_size)
        lat = np.arange(-90 + grid_size / 2, 90, step=grid_size)

        return lat, lon


if __name__ == '__main__':
    # CS = np.array([23, 56, 78, 90, 32, 34])
    # GeoMathKit.CS_1dTo2d(CS)

    GeoMathKit.global_equal_distance_grid(1)
