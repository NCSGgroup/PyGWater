import numpy as np
from src_EarthModel.LoveNumber import LoveNumber, LLN_Data, LLN_variable, Frame
from src_Basis.GreenFunction.DiskLoad import DiskLoad
from src_auxiliary.EnumClasses import Displacement
from src_auxiliary.GeoMathKit import GeoMathKit
from src_LoadMesh.RegularGrid import regular
from pathlib import Path
from src_auxiliary.FileTool import FileTool


def template_name(postfix, variable=Displacement.Vertical):
    """file name of a CD-LGF template; vertical keeps the original name LGD_<postfix>.npy"""
    if variable == Displacement.Vertical:
        return 'LGD_%s.npy' % postfix
    return 'LGD_%s_%s.npy' % (variable.name, postfix)


class cluster_disk:

    def __init__(self, lln: LoveNumber, verbose=False):
        self._small_disk = None
        self._small_grid = None
        self._DL = None
        self._lln = lln
        self._verbose = verbose

    def configure(self, grid: dict, small_disk=0.1):
        """

        :param grid: lat, lon, delta_lat, delta_lon, unit is degree; this is assumed as the center of a single grid
        :param small_disk: SUBGRID SIZE [degree] used to partition the cell (e.g. 0.1, as recommended in
                           Yang et al., 2025). Each sub-cell is a disk of equal area, whose angular RADIUS is
                           ~ small_disk/sqrt(pi) at the equator (0.056 deg for 0.1) and smaller towards the poles.
                           The Love-number expansion should satisfy Eq. 24 for that radius, see
                           DiskLoad.required_lmax().
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

    def template(self, points: dict, variable=Displacement.Vertical, postfix='1', out_dir=None):
        """
        Generate a template, uniform disk with thickness of 1 meter
        :param points: points (many points) to be evaluated. ['lat':, 'lon':], unit is degree
        :param out_dir: folder for the template; default <project>/res/LGD_template
        :param variable: Displacement.Vertical, Geoheight, North or East (one template file per variable)
        :return: the template (saved as LGD_<postfix>.npy for Vertical, LGD_<Variable>_<postfix>.npy otherwise)
        """
        grid_lat, grid_lon = self._small_grid
        dl = self._DL
        small_disk = self._small_disk

        RR = regular.grid2radius(lat_center=grid_lat, grid_size=small_disk)[0]
        dis = np.zeros_like(points['lat'])

        for id, rr in enumerate(list(RR)):
            if self._verbose:
                print(id)
            '''calculate the displacement: Vertical, Geoheight, Horizontal, North or East'''
            dis += GeoMathKit.kernel_at_points(dl[rr], grid_lat[id], grid_lon[id], points['lat'], points['lon'],
                                               variable)

            '''release memory'''
            dl[rr].release()

            pass

        out_dir = FileTool.get_project_dir('res/LGD_template') if out_dir is None else Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        np.save(out_dir / template_name(postfix, variable), dis)

        return dis



class cluster_disk_fast:
    """
    Approximate version of cluster_disk: beyond 8 deg (measured from sub-disk no. 100) the response of EVERY sub-disk
    is replaced by that of sub-disk no. 100, so the cell must contain more than 100 sub-cells. Kept for
    reproducibility; for new work use cluster_disk or GreenFunctionBasis(kernel='cluster_disk'), which is exact
    and faster (tabulated kernels, src_Basis/GreenFunction/DiskLoad.DiskLoadTable).
    """

    def __init__(self, lln: LoveNumber, verbose=False):
        self._small_disk = None
        self._small_grid = None
        self._DL = None
        self._lln = lln
        self._verbose = verbose

    def configure(self, grid: dict, small_disk=0.1):
        """

        :param grid: lat, lon, delta_lat, delta_lon, unit is degree; this is assumed as the center of a single grid
        :param small_disk: SUBGRID SIZE [degree] used to partition the cell (e.g. 0.1, as recommended in
                           Yang et al., 2025). Each sub-cell is a disk of equal area, whose angular RADIUS is
                           ~ small_disk/sqrt(pi) at the equator (0.056 deg for 0.1) and smaller towards the poles.
                           The Love-number expansion should satisfy Eq. 24 for that radius, see
                           DiskLoad.required_lmax().
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

    def template(self, points: dict, variable=Displacement.Vertical, postfix='1', out_dir=None):
        """
        Generate a template, uniform disk with thickness of 1 meter
        :param points: points (many points) to be evaluated. ['lat':, 'lon':], unit is degree
        :param out_dir: folder for the template; default <project>/res/LGD_template
        :param variable: Displacement.Vertical, Geoheight, North or East (one template file per variable)
        :return: the template (saved as LGD_<postfix>.npy for Vertical, LGD_<Variable>_<postfix>.npy otherwise)
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
        dist0 = GeoMathKit.kernel_at_points(dl[RR[id]], grid_lat[id], grid_lon[id], points['lat'], points['lon'],
                                            variable)
        threshold = 8
        ff1 = theta0 < threshold
        ff2 = ~ff1

        for id, rr in enumerate(list(RR)):
            if self._verbose:
                print(id)
            '''calculate the displacement'''
            dis[ff1] += GeoMathKit.kernel_at_points(dl[rr], grid_lat[id], grid_lon[id], points['lat'][ff1],
                                                    points['lon'][ff1], variable)
            dis[ff2] += dist0[ff2]

            '''release memory'''
            dl[rr].release()

            pass

        out_dir = FileTool.get_project_dir('res/LGD_template') if out_dir is None else Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        np.save(out_dir / template_name(postfix, variable), dis)

        return dis



class GFA_displacement_from_template:
    """
    Green Function Approach (GFA): based on disk load green function;
    """

    def __init__(self, template_dir: str, verbose=False):
        self._path = Path(template_dir)
        self._verbose = verbose

        pass

    def evaluation_general(self, grids: dict, variable=Displacement.Vertical):
        """
        Evaluation of desired variable at specified points
        :param grids: 'EWH': ..; 'Template': ...
        :param variable: which templates to use (Vertical, Geoheight, North, East)
        :return:
        """
        dis = np.zeros_like(grids['EWH'])
        for id, rr in enumerate(list(grids['EWH'])):
            if self._verbose:
                print(id)
            tl = np.load(self._path / template_name(grids['Template'][id], variable))
            dis += tl * rr

        return dis

    def evaluation_regular_grid(self, grids: dict, dim=90, variable=Displacement.Vertical):
        """
        Evaluation of desired variable at specified points
        :param grids: 'EWH': ..; 'Template': ...
        :param dim: number of LATITUDE rows of the evaluation grid (the template is reshaped to (dim, nlon) and
                    rolled along longitude; valid for North/East as well)
        :param variable: which templates to use (Vertical, Geoheight, North, East)
        :return:
        """
        dis = np.zeros_like(grids['EWH'])
        flag = -1000
        for id, rr in enumerate(list(grids['EWH'])):
            if self._verbose:
                print(id)
            a = grids['Template'][id]
            if a != flag:
                tl = np.load(self._path / template_name(a, variable))
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
    a = np.load(FileTool.get_project_dir('res/LGD_template/LGD_0.npy'))
    b = np.load(FileTool.get_project_dir('res/LGD_template/LGD_1.npy'))
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
