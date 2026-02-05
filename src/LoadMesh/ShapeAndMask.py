import numpy as np
import geopandas as gpd
from shapely import box
from pathlib import Path
import shapely
import shapely.vectorized
import pandas as pd
import h5py


def box2shp(box_area=(70.1, 33.9, -11.1, 45.1)):
    """
    Generate a shp file from a box boundary.

    For example: box = [76.1, 33.9, -11.1, 45.1] ==> [up (lat), down (lat), left (lon), right (lon)]
    """

    poly = box(xmin=box_area[2], ymin=box_area[0], xmax=box_area[3], ymax=box_area[1])
    d = {'ID': [0], 'geometry': [poly]}
    gdf = gpd.GeoDataFrame(d, crs='epsg:4326')

    # '''visualization'''
    # import pygmt
    # import geopandas as gpd
    # fig = pygmt.Figure()
    # pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.2)
    # pygmt.config(FONT_ANNOT='10p', COLOR_NAN='white')
    # pygmt.makecpt(cmap='wysiwyg', series=[0, 56], background='o')
    # region = 'g'
    # pj = "Q30/-20/12c"
    # fig.basemap(region=region, projection=pj,
    #             frame=['WSne', 'xa10f5+lLongitude (\\260 E)', 'ya5f5+lLatitude (\\260 N)'])
    #
    # fig.coast(shorelines="1/0.2p", region=region, projection=pj, water="skyblue")
    # fig.plot(data=gdf.boundary, pen="0.5p,red", projection=pj)
    # fig.show()

    return gdf


class basin2grid_shp:
    """
    Given a basin, to subdivide the basin into grid by grid.
    """

    def __init__(self, grid=(3, 3)):
        self.new_shp = None

        sub_basin = grid

        N = 180 // sub_basin[0]
        M = 360 // sub_basin[1]

        num_sub_basins = N * M

        lat_lu = np.array([90 - i * sub_basin[0] for i in range(N)])
        lon_lu = np.array([-180 + i * sub_basin[1] for i in range(M)])

        lat_ru = lat_lu.copy()
        lon_ru = lon_lu + sub_basin[1]

        lat_ld = lat_lu - sub_basin[0]
        lon_ld = lon_lu.copy()

        # lat_rd = lat_lu - sub_basin[0]
        # lon_rd = lon_lu + sub_basin[1]

        # lon_lu, lat_lu = np.meshgrid(lon_lu, lat_lu)
        lon_ru, lat_ru = np.meshgrid(lon_ru, lat_ru)
        lon_ld, lat_ld = np.meshgrid(lon_ld, lat_ld)
        # lon_rd, lat_rd = np.meshgrid(lon_rd, lat_rd)

        poly = box(xmin=lon_ld, ymin=lat_ld, xmax=lon_ru, ymax=lat_ru)
        ID = [i for i in range(1, num_sub_basins + 1)]
        d = {'ID': ID, 'geometry': list(poly.flatten())}
        self.gdf = gpd.GeoDataFrame(d, crs='epsg:4326')
        pass

    def load_original_shp(self, original_shp):
        if isinstance(original_shp, str):
            self.__old_shp = gpd.read_file(original_shp).to_crs(crs='epsg:4326')
        else:
            self.__old_shp = original_shp

        return self

    def create_shp(self):
        basin_shp = self.__old_shp

        if basin_shp.unary_union.geom_type == 'MultiPolygon':
            bb = max(basin_shp.unary_union.geoms, key=lambda a: a.area)
        else:
            bb = basin_shp.unary_union

        basin_all = gpd.GeoDataFrame({'geometry': gpd.GeoSeries(bb)}, crs='epsg:4326')

        grid = self.gdf

        index = shapely.intersects(bb, grid.geometry).values

        new = []

        for i in np.arange(len(index)):
            if not index[i]:
                continue

            mm = grid[grid.ID == i + 1]
            mm = mm.drop(columns=['ID'])
            new.append(mm.overlay(basin_all, how='intersection', keep_geom_type=True))

        gdf = gpd.GeoDataFrame(pd.concat(new))
        ID = [i for i in range(1, len(gdf) + 1)]
        d = {'ID': ID, 'geometry': gdf.geometry}
        new_shp = gpd.GeoDataFrame(d, crs='epsg:4326')

        self.new_shp = new_shp

        return self

    def save2file(self, new_basin_name='new', out_dir='../temp'):
        dp= Path(out_dir)/new_basin_name
        if not dp.exists():
            dp.mkdir()
        self.new_shp.to_file(dp / ('%s.shp' % new_basin_name))
        pass

    def delete_tiny_grids(self):
        """
        The aim of this function is to delete tiny grids
        """
        gdf = self.new_shp

        '''requirement of minimal area'''
        propotion = 8  # This is totally empirical
        gdf = gdf[gdf.area >= (gdf.area.values.max()) / propotion]

        '''must not be a multiPolygon'''
        # gdf = gdf[gdf.type == 'Polygon']

        '''requirement of the centroid, e.g., must be inside a box'''
        # lon_ld = -130
        # lat_ld = 20
        # lon_ru = -65
        # lat_ru = 50
        # poly = box(xmin=lon_ld, ymin=lat_ld, xmax=lon_ru, ymax=lat_ru)
        # gdf = gdf[poly.contains(gdf.centroid)]

        '''reorganization'''
        gdf['ID'] = np.arange(start=1, stop=gdf.shape[0] + 1)

        self.new_shp = gdf

        return self

    def delete_ocean_grid(self):
        """
        Only keep the land grid. This calculation is based on the high-resolution (0.1-degree) land mask data (coming from
        the W3RA model). #todo: to inspect and likely improve this land mask in future
        :return:
        """
        model_land_mask = '../../data/land_mask/land_mask_res0.1.h5'
        mask = h5py.File(model_land_mask, 'r')['mask'][:-1, :]
        gdf = self.new_shp

        invalids = []

        for id in range(1, gdf.ID.size + 1):
            tt = gdf[gdf.ID == id]
            minx = int(float(tt.bounds.minx) / 0.1) + 1800
            maxx = int(float(tt.bounds.maxx) / 0.1) + 1800
            maxy = 900 - int(float(tt.bounds.miny) / 0.1)
            miny = 900 - int(float(tt.bounds.maxy) / 0.1)
            vv = np.sum(mask[miny:maxy + 1, minx:maxx + 1])
            tg = mask[miny:maxy + 1, minx:maxx + 1].size

            flag = True

            '''in case the land area is too small'''
            if vv / tg < 0.16:
                flag = False

            invalids.append(flag)

            pass

        n = gdf[invalids]
        n.loc[:, 'ID'] = np.arange(np.sum(np.array(invalids))) + 1
        self.new_shp = n
        return self


class basinMask:
    """
    Get the global mask, given the basin. The mask resolution must be defined.
    """

    def __init__(self, mask_res, new_mask_name='MDB', save_dir='../data/basin/mask'):
        self._res = mask_res
        self._basin_name = new_mask_name
        self._save_dir = Path(save_dir)

        self.box_mask = None
        self.__shp = None
        pass

    def configureBox(self, box: dir):
        """
        box = {
        "lat": [
            -9.9,
            -43.8
        ],
        "lon": [
            112.4,
            154.3
        ]}
        """
        if box is None:
            return self

        res = self._res
        err = res / 10
        lat = np.arange(90 - res / 2, -90 + res / 2 - err, -res)
        lon = np.arange(-180 + res / 2, 180 - res / 2 + err, res)

        lati = [np.argmin(np.fabs(lat - max(box['lat']))),
                np.argmin(np.fabs(lat - min(box['lat'])))]
        loni = [np.argmin(np.fabs(lon - min(box['lon']))),
                np.argmin(np.fabs(lon - max(box['lon'])))]
        id = [lati[0], lati[1] + 1, loni[0], loni[1] + 1]

        lon, lat = np.meshgrid(lon, lat)
        mask = np.zeros(np.shape(lon))

        mask[id[0]:id[1], id[2]:id[3]] = 1

        self.box_mask = mask

        return self

    def load_shp(self, shp):
        if isinstance(shp, str):
            self.__shp = gpd.read_file(shp).to_crs(crs='epsg:4326')
        else:
            self.__shp = shp

        return self

    def shp_to_mask(self):
        res = self._res
        basin_name = self._basin_name

        h5fn = str(self._save_dir / ('%s_res_%s.h5' % (basin_name, res)))
        hf = h5py.File(h5fn, 'w')

        mask_basin = {}
        gdf = self.__shp

        err = res / 10
        lat = np.arange(90 - res / 2, -90 + res / 2 - err, -res)
        lon = np.arange(-180 + res / 2, 180 - res / 2 + err, res)

        lon, lat = np.meshgrid(lon, lat)

        mask_all = np.zeros(np.shape(lat))
        for id in np.arange(gdf.ID.size) + 1:
            # if id !=3: continue
            bd1 = gdf[gdf.ID == id]
            mask1 = shapely.vectorized.touches(bd1.geometry.item(), lon, lat)
            mask2 = shapely.vectorized.contains(bd1.geometry.item(), lon, lat)
            mask3 = mask1 + mask2
            # mask3=mask2
            # print('Mask for %s' % bd1.ID.values[0])
            if self.box_mask is not None:
                mask3 = (mask3 * self.box_mask).astype(bool)

            mask_basin[bd1.ID.values[0]] = mask3
            mask_all += mask3

            hf.create_dataset('sub_basin_%s' % id, data=mask3.astype(int))
            pass

        hf.create_dataset('basin', data=mask_all.astype(int))
        hf.create_dataset('lat', data=lat)
        hf.create_dataset('lon', data=lon)
        hf.create_dataset('resolution', data=res)
        hf.close()

        mask_basin[0] = mask_all.astype(bool)

        return mask_basin, lat, lon


class getMask_flow:

    def __init__(self, basin_name='exp', coarse_cell=3, small_cell=0.1, is_ocean_removed=False):
        self.__new_basin_name = basin_name
        self.coarse_cell = coarse_cell
        self.small_cell = small_cell
        self.__is_ocean_removed = is_ocean_removed

        self.shp = None
        self.mask, self.lat, self.lon = None, None, None

        assert small_cell <= coarse_cell
        pass

    def getShp(self, box_area=(71.6, 36.1, -11.1, 42.1), shape_file=None):
        if shape_file is not None:
            self.shp = gpd.read_file(shape_file).to_crs(crs='epsg:4326')
            return self

        '''get the frame/exterior of the study area'''
        fr = box2shp(box_area=box_area)

        '''divide the study area into smaller cells at a resolution of coarse cell'''
        bg = basin2grid_shp(grid=(self.coarse_cell, self.coarse_cell))

        bg.load_original_shp(original_shp=fr).create_shp().delete_tiny_grids()

        if self.__is_ocean_removed:
            bg.delete_ocean_grid()

        self.shp = bg.new_shp
        bg.save2file(new_basin_name=self.__new_basin_name, out_dir='../../res/shp')
        return self

    def shp2mask(self):
        bm = basinMask(mask_res=self.small_cell, new_mask_name=self.__new_basin_name, save_dir='../../res/mask')

        self.mask, self.lat, self.lon = bm.load_shp(shp=self.shp).shp_to_mask()

        pass



def demo1():
    gmf = getMask_flow(basin_name='exp', coarse_cell=3, small_cell=0.5, is_ocean_removed=True)

    gmf.getShp(box_area=(71.6, 36.1, -11.1, 42.1)).shp2mask()
    pass


def demo_visualization():
    box_area = (71.6, 36.1, -11.1, 42.1)
    gdf =gpd.read_file(filename='/home/user/codes/PyGWater/res/shp/exp/exp.shp')

    import pygmt

    fig = pygmt.Figure()
    pygmt.config(MAP_HEADING_OFFSET=0, MAP_TITLE_OFFSET=-0.2)
    pygmt.config(FONT_ANNOT='10p', COLOR_NAN='white')
    pygmt.makecpt(cmap='wysiwyg', series=[0, 56], background='o')

    region = [box_area[2] - 5, box_area[3] + 5, box_area[1] - 5, box_area[0] + 5]
    pj = "Q30/-20/12c"
    fig.basemap(region=region, projection=pj,
                frame=['WSne', 'xa10f5+lLongitude (\\260 E)', 'ya5f5+lLatitude (\\260 N)'])


    # fig.show()

    '''mask visualization'''
    h5 = h5py.File('/home/user/codes/PyGWater/res/mask/exp_res_0.5.h5','r')
    lat = h5['lat'][:]
    lon = h5['lon'][:]
    mask =h5['basin'][:]
    res = h5['resolution'][()].item()



    # mask = h5['sub_basin_23'][:]
    pygmt.makecpt(cmap='polar', series=[-1, 1], background='o')
    data_region = [min(lon[0]), max(lon[0]), min(lat[:,0]), max(lat[:,0])]

    ee = pygmt.xyz2grd(y=lat.flatten(), x=lon.flatten(), z=mask.flatten(),
                       spacing=(res, res), region=data_region)

    # pj = 'Q10c'
    fig.grdimage(
        grid=ee,
        cmap=True,
        # frame=['xaf', 'yaf'] + ['+t %s v.s. %s (%.1f mm)' % (dt1.name, miss2.name+'_'+dt2.name, statistic)],
        frame=['xaf', 'yaf'],
        # frame=['xaf', 'yaf'] + ['+t %s v.s. %s' % (miss1.name + '_' + dt1.name, miss2.name + '_' + dt2.name)],
        # dpi=100,
        projection=pj,
        interpolation='n',
        # nan_transparent=True,
        region=region
    )

    fig.coast(shorelines="1/0.2p", region=region, projection=pj)

    for i in range(1, gdf.shape[0] + 1):
        xy = gdf[gdf.ID == i].centroid
        fig.text(x=xy.x, y=xy.y, text="%s" % i, font='7p,black')

    fig.plot(data=gdf.boundary, pen="0.5p,blue", projection=pj)

    fig.show()
    pass


if __name__ == '__main__':
    # demo1()
    demo_visualization()
