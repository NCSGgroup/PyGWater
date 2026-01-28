import numpy as np
from src.Auxiliary.constants import EarthConstant
from src.Auxiliary.GeoMathKit import GeoMathKit


class regular:

    @staticmethod
    def grid2radius(lat_center, grid_size):
        """
        calculate the radius of disk load that equals to the area of pixel
        :param lat_center: center of the pixel, latitude, [degree]
        :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
        :return: theta_radius [degree], length_radius [m]
        """
        from src.Auxiliary.constants import EarthConstant
        example = 30  # longitude, but indeed the choice could be arbitrary.

        a = GeoMathKit.angular_distance(point1_lat=lat_center + grid_size / 2, point1_lon=example,
                                        point2_lat=lat_center - grid_size / 2, point2_lon=example)

        b = GeoMathKit.angular_distance(point1_lat=lat_center, point1_lon=example, point2_lat=lat_center,
                                        point2_lon=example + grid_size)

        # print(np.deg2rad(a) * np.deg2rad(b))
        r = np.sqrt(np.deg2rad(a) * np.deg2rad(b) / np.pi)

        return np.rad2deg(r), EarthConstant.radiusm * r

    @staticmethod
    def grid2radius_type2(lat_center, grid_size):
        """
        calculate the radius of disk load that equals to the area of pixel
        :param lat_center: center of the pixel, latitude, [degree]
        :param grid_size: equal-distance grid, grid interval, e.g., 1 degree
        :return: theta_radius [degree], length_radius [m]
        """
        from src.Auxiliary.constants import EarthConstant

        area = np.cos(np.deg2rad(lat_center)) * np.deg2rad(grid_size) ** 2

        # print(area)
        x = 1 - area / (2 * np.pi)

        r = np.arccos(x)

        return np.rad2deg(r), EarthConstant.radiusm * r


    @staticmethod
    def global_grid_for_pointmass(resolution=0.5, Earth_radius=EarthConstant.radiusm):
        """

        :param Earth_radius: [meter]
        :param resolution: [degree]
        :return: grid, dict
        """

        lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=resolution)
        lon, lat = np.meshgrid(lon, lat)

        area = np.cos(np.deg2rad(lat)) * np.deg2rad(resolution) ** 2 * Earth_radius ** 2

        grids = {
            'lat': lat.flatten(),
            'lon': lon.flatten(),
            'area': area.flatten()
        }

        return grids