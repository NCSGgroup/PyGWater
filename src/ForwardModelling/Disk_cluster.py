import numpy as np
from src.GreenFunction.LLN import LoveNumber, LLN_Data, LLN_variable, Frame
from src.GreenFunction.DiskLoad import DiskLoad
from src.Auxiliary.EnumClasses import Displacement
from src.Auxiliary.GeoMathKit import GeoMathKit
from src.LoadMesh.RegularGrid import regular
from pathlib import Path


class cluster_disk:

    def __init__(self, lln: LoveNumber):
        self._small_disk = None
        self._small_grid = None
        self._DL = None
        self._lln = lln

    def configure(self, grid: dict, small_disk=0.1):
        """

        :param grid: lat, lon, delta_lat, delta_lon, unit is degree; this is assumed as the center of a single grid
        :param small_disk: radius of small disks to be partitioned from the specified grid.
        :return:
        """
        self._small_disk = small_disk

        lat = [grid['lat'] - grid['delta_lat'] / 2, grid['lat'] + grid['delta_lat'] / 2]
        lon = [grid['lon'] - grid['delta_lon'] / 2, grid['lon'] + grid['delta_lon'] / 2]

        new_lat = np.arange(lat[0] + small_disk / 2, lat[1], small_disk)
        new_lon = np.arange(lon[0] + small_disk / 2, lon[1], small_disk)

        new_lat, new_lon = np.meshgrid(new_lat, new_lon)
        new_lat = new_lat.flatten()
        new_lon = new_lon.flatten()

        RR = regular.grid2radius(lat_center=new_lat, grid_size=small_disk)[0]

        Unit_Disk_loads = {}
        for rr in RR:
            if rr not in Unit_Disk_loads.keys():
                Unit_Disk_loads[rr] = DiskLoad(lln=self._lln, radius=rr, thickness=1)

        self._DL = Unit_Disk_loads
        self._small_grid = [new_lat, new_lon]

        pass

    def template(self, points: dict, variable=Displacement.Vertical, postfix='1'):
        """
        Generate a template, uniform disk with thickness of 1 meter
        :param points: points (many points) to be evaluated. ['lat':, 'lon':], unit is degree
        :return:
        """
        grid_lat, grid_lon = self._small_grid
        dl = self._DL
        small_disk = self._small_disk

        RR = regular.grid2radius(lat_center=grid_lat, grid_size=small_disk)[0]
        dis = np.zeros_like(points['lat'])

        for id, rr in enumerate(list(RR)):
            print(id)
            theta = GeoMathKit.angular_distance(grid_lat[id], grid_lon[id], points['lat'],
                                                points['lon'])

            '''To avoid singularity'''
            theta[theta == 0] += 1e-6
            theta[theta == 180] += 1e-6

            '''calculate the displacement'''
            dis += self._getFunc(DL=dl[rr].configure(residence=theta), variable=variable)

            '''release memory'''
            dl[rr].release()

            pass

        np.save('../../res/LGD_template/LGD_%s.npy' % postfix, dis)

        pass

    def _getFunc(self, variable: Displacement, DL: DiskLoad):
        if variable == Displacement.Vertical:
            return DL.getVertical()
        elif variable == Displacement.Horizontal:
            return DL.getHorizental()
        elif variable == Displacement.Geoheight:
            return DL.getGeoidheight()


class cluster_disk_fast:

    def __init__(self, lln: LoveNumber):
        self._small_disk = None
        self._small_grid = None
        self._DL = None
        self._lln = lln

    def configure(self, grid: dict, small_disk=0.1):
        """

        :param grid: lat, lon, delta_lat, delta_lon, unit is degree; this is assumed as the center of a single grid
        :param small_disk: radius of small disks to be partitioned from the specified grid.
        :return:
        """
        self._small_disk = small_disk

        lat = [grid['lat'] - grid['delta_lat'] / 2, grid['lat'] + grid['delta_lat'] / 2]
        lon = [grid['lon'] - grid['delta_lon'] / 2, grid['lon'] + grid['delta_lon'] / 2]

        new_lat = np.arange(lat[0] + small_disk / 2, lat[1], small_disk)
        new_lon = np.arange(lon[0] + small_disk / 2, lon[1], small_disk)

        new_lat, new_lon = np.meshgrid(new_lat, new_lon)
        new_lat = new_lat.flatten()
        new_lon = new_lon.flatten()

        RR = regular.grid2radius(lat_center=new_lat, grid_size=small_disk)[0]

        Unit_Disk_loads = {}
        for rr in RR:
            if rr not in Unit_Disk_loads.keys():
                Unit_Disk_loads[rr] = DiskLoad(lln=self._lln, radius=rr, thickness=1)

        self._DL = Unit_Disk_loads
        self._small_grid = [new_lat, new_lon]

        pass

    def template(self, points: dict, variable=Displacement.Vertical, postfix='1'):
        """
        Generate a template, uniform disk with thickness of 1 meter
        :param points: points (many points) to be evaluated. ['lat':, 'lon':], unit is degree
        :return:
        """
        grid_lat, grid_lon = self._small_grid
        dl = self._DL
        small_disk = self._small_disk

        RR = regular.grid2radius(lat_center=grid_lat, grid_size=small_disk)[0]
        dis = np.zeros_like(points['lat'])

        '''trial to test how far these points are'''
        id = 100
        theta0 = GeoMathKit.angular_distance(grid_lat[id], grid_lon[id], points['lat'],
                                             points['lon'])
        dist0 = self._getFunc(DL=dl[RR[id]].configure(residence=theta0), variable=variable)
        threshold = 8
        ff1 = theta0 < threshold
        ff2 = ~ff1

        for id, rr in enumerate(list(RR)):
            print(id)
            theta = GeoMathKit.angular_distance(grid_lat[id], grid_lon[id], points['lat'],
                                                points['lon'])

            '''To avoid singularity'''
            theta[theta == 0] += 1e-6
            theta[theta == 180] += 1e-6

            '''calculate the displacement'''
            dis[ff1] += self._getFunc(DL=dl[rr].configure(residence=theta[ff1]), variable=variable)
            dis[ff2] += dist0[ff2]

            '''release memory'''
            dl[rr].release()

            pass

        np.save('../../res/LGD_template/LGD_%s.npy' % postfix, dis)

        pass

    def _getFunc(self, variable: Displacement, DL: DiskLoad):
        if variable == Displacement.Vertical:
            return DL.getVertical()
        elif variable == Displacement.Horizontal:
            return DL.getHorizental()
        elif variable == Displacement.Geoheight:
            return DL.getGeoidheight()


class GFA_displacement_from_template:
    """
    Green Function Approach (GFA): based on disk load green function;
    """

    def __init__(self, template_dir: str):
        self._path = Path(template_dir)

        pass

    def evaluation_general(self, grids: dict):
        """
        Evaluation of desired variable at specified points
        :param grids: 'EWH': ..; 'Template': ...
        :return:
        """
        dis = np.zeros_like(grids['EWH'])
        for id, rr in enumerate(list(grids['EWH'])):
            print(id)
            tl = np.load(self._path / ('LGD_%s.npy' % grids['Template'][id]))
            dis += tl * rr

        return dis

    def evaluation_regular_grid(self, grids: dict, dim=90):
        """
        Evaluation of desired variable at specified points
        :param grids: 'EWH': ..; 'Template': ...
        :param dim: number of points along the longitude direction
        :return:
        """
        dis = np.zeros_like(grids['EWH'])
        flag = -1000
        for id, rr in enumerate(list(grids['EWH'])):
            print(id)
            a = grids['Template'][id]
            if a != flag:
                tl = np.load(self._path / ('LGD_%s.npy' % a))
                flag = a
            else:
                tl = np.roll(tl.reshape((dim, -1)), shift=1, axis=1)
                tl = tl.flatten()

            dis += tl[:, None] * rr

        return dis


def demo1():
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    cd = cluster_disk(lln=lln)

    res = 2
    small_radius = 0.1
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    for id in np.arange(int(180 / res)):
        cd.configure(grid={'lat': lat[id, 0], 'lon': lon[id, 0], 'delta_lat': res, 'delta_lon': res},
                     small_disk=small_radius)
        points = {
            'lat': lat.flatten(),
            'lon': lon.flatten(),
        }
        cd.template(points=points, postfix='%s' % id)
    pass


def demo2():
    a = np.load('../../res/LGD_template/LGD_0.npy')
    b = np.load('../../res/LGD_template/LGD_1.npy')
    pass

def demo3():
    lln = LoveNumber().config(lmax=10000, method=LLN_Data.PREM).get_Love_number()
    cd = cluster_disk(lln=lln)

    res = 2
    small_radius = 0.05
    lat, lon = GeoMathKit.global_equal_distance_grid(grid_size=res)
    lon, lat = np.meshgrid(lon, lat)

    id = 45
    cd.configure(grid={'lat': lat[id, 0], 'lon': lon[id, 0], 'delta_lat': res, 'delta_lon': res},
                 small_disk=small_radius)
    points = {
        'lat': lat[id, 0]*np.ones(1000),
        'lon': lon[id, 0]+np.linspace(start=0.0001, stop=res*5, num=1000),
    }
    cd.template(points=points, postfix='fig14_res_0.05')


if __name__ == '__main__':
    demo3()
