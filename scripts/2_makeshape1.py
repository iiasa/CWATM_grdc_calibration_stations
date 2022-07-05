"""
# -------------------------------------------------------------------------
# Name:        Use MERIT coordinates of upstream area to create shapefiles
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to create shapefiles  from rivernetwork MERIT data
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022

# input:  grdc_MERIT_1111.txt  station with new location fitted to merit UPA
# output: grdc_MERIT_2222.txt: station with new location fitted to merit UPA and shapefile
#       shapefiles of basins in : ashapex_merit  e.g. grdc_basin_merit_1104800.shp

No: Number from 1 ...
GRDC_No: GRDC number
lat: original latitude from GRDC metafile
lon: original longituted from GRDC metafile
area; provided basin area from GRDC metafile
newlat: corrected latitude based on MERIT UPA dataset
newlon: corrected longitute based on MERIT UPA dataset
newarea: basin area based on MERIT UPA dataset
shapearea: area of the shape calculated with geopandas directly from shape
Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar

# uses
# grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
# grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
# Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
"""

# ----------------------------------------------------------------------


# import pyflwdir, some dependencies
import geopandas as gp
import numpy as np
import rasterio
from rasterio import features
import pyflwdir

import sys
import os
import os.path
import warnings

# local convenience methods (see utils.py script in notebooks folder)
#from utils import vectorize  # convenience method to vectorize rasters
#from utils import quickplot, colors, cm  # data specific quick plot method


# convenience method for vectorizing a raster
def vectorize(data, nodata, transform, name="value"):
    feats_gen = features.shapes(
        data,
        mask=data != nodata,
        transform=transform,
        connectivity=8,
    )
    feats = [
        {"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)
    ]

    # parse to geopandas for plotting / writing to file
    gdf = gp.GeoDataFrame.from_features(feats, crs=crs)
    gdf[name] = gdf[name].astype(data.dtype)
    return gdf

#---------------------------------------------------


#-----------------------------------------------

def getBasinInd(grdc_shape,newshape):

    p_inter = gp.overlay(grdc_shape, newshape, how='intersection')
    p_union = gp.overlay(grdc_shape, newshape, how='union')

    pint_area = p_inter.area.sum()
    puni_area = p_union.area.sum()
    indicator = pint_area / puni_area

    # s1 = grdc_shape.symmetric_difference(p2, align=False)
    # sym_min = s1.area

    return indicator

#----------------------------------------------------

#https://gis.stackexchange.com/questions/413349/calculating-area-of-lat-lon-polygons-without-transformation-using-geopandas
def gpd_geographic_area(geodf):
    if not geodf.crs and geodf.crs.is_geographic:
        raise TypeError('geodataframe should have geographic coordinate system')

    geod = geodf.crs.get_geod()

    def area_calc(geom):
        if geom.geom_type not in ['MultiPolygon', 'Polygon']:
            return np.nan

        # For MultiPolygon do each separately
        if geom.geom_type == 'MultiPolygon':
            return np.sum([area_calc(p) for p in geom.geoms])

        # orient to ensure a counter-clockwise traversal.
        # See https://pyproj4.github.io/pyproj/stable/api/geod.html
        # geometry_area_perimeter returns (area, perimeter)
        return geod.geometry_area_perimeter(orient(geom, 1))[0]

    return geodf.geometry.apply(area_calc)





def line_integral_polygon_area(geom, radius=6378137):
    """
    Computes area of spherical polygon, assuming spherical Earth.
    Returns result in ratio of the sphere's area if the radius is specified.
    Otherwise, in the units of provided radius.
    lats and lons are in degrees.

    from https://stackoverflow.com/a/61184491/6615512
    """
    if geom.geom_type not in ['MultiPolygon', 'Polygon']:
        return np.nan

    # For MultiPolygon do each separately
    if geom.geom_type == 'MultiPolygon':
        return np.sum([line_integral_polygon_area(p) for p in geom.geoms])

    # parse out interior rings when present. These are "holes" in polygons.
    if len(geom.interiors) > 0:
        interior_area = np.sum([line_integral_polygon_area(Polygon(g)) for g in geom.interiors])
        geom = Polygon(geom.exterior)
    else:
        interior_area = 0

    # Convert shapely polygon to a 2 column numpy array of lat/lon coordinates.
    geom = np.array(geom.boundary.coords)

    lats = np.deg2rad(geom[:, 1])
    lons = np.deg2rad(geom[:, 0])

    # Line integral based on Green's Theorem, assumes spherical Earth

    # close polygon
    if lats[0] != lats[-1]:
        lats = np.append(lats, lats[0])
        lons = np.append(lons, lons[0])

    # colatitudes relative to (0,0)
    a = np.sin(lats / 2) ** 2 + np.cos(lats) * np.sin(lons / 2) ** 2
    colat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # azimuths relative to (0,0)
    az = np.arctan2(np.cos(lats) * np.sin(lons), np.sin(lats)) % (2 * np.pi)

    # Calculate diffs
    # daz = np.diff(az) % (2*pi)
    daz = np.diff(az)
    daz = (daz + np.pi) % (2 * np.pi) - np.pi

    deltas = np.diff(colat) / 2
    colat = colat[0:-1] + deltas

    # Perform integral
    integrands = (1 - np.cos(colat)) * daz

    # Integrate
    area = abs(sum(integrands)) / (4 * np.pi)

    area = min(area, 1 - area)
    if radius is not None:  # return in units of radius
        return (area * 4 * np.pi * radius ** 2) - interior_area
    else:  # return in ratio of sphere total area
        return area - interior_area


# a wrapper to apply the method to a geo data.frame
def gpd_geographic_area_line_integral(geodf):
    return geodf.geometry.apply(line_integral_polygon_area)

# --------------------------------------------

#ldd = "E:/Data/GRDC_new/testEurope/flow_danube.tif"
#ldd = "E:/Data/GRDC_new/testEurope/flow_irEurope5.tif"
# read and parse data

# c = np.hstack((a, b))
#c = np.vstack((a, b))
#l1 = "E:/Data/GRDC_new/testEurope/ftw1.npy"
#with open(l1, "rb") as handle:
#    k1 = pickle.load(handle)
#flw = pyflwdir.load('E:/Data/GRDC_new/testEurope/flw.pkl')

warnings.filterwarnings("ignore")

root =  "P:/watproject/Datasets/MERIT_yamazaki/dir_"
rootshape = "P:/watproject/Datasets/MERIT_yamazaki/ashapex_merit/grdc_basin_merit_"
rootgrdcshape = "P:/watproject/Datasets/MERIT_yamazaki/ashape_merit/grdc_basin_merit_"
rootgrdcshapeold = "P:/watproject/Datasets/MERIT_yamazaki/grdc_station_upstream/grdc_basins_smoothed_md_no_"

sizenoarea = 5
#----------------------------------------------------------

#grdc_stations = "makeshape_all.txt"
grdc_stations = "results/grdc_Merit_1111.txt"
f = open(grdc_stations, "r") 
grdc = f.readlines()[1:]
f.close()

grdcallfile = "grdc_2022_10701.txt"
grdc_d = {}
with open(grdcallfile) as f:
    for line in f:
       l1 = line.split("\t")
       grdc_d[l1[0]] = l1
headgrdc = grdc_d["grdc_no"]
headgrdc[-1] = headgrdc[-1][0:-1]

# translation file new GRDC ID - old GRDC ID
lehnerfile = "grdc2022_lehner2012.txt"
lehner = {}
with open(lehnerfile) as f:
    for line in f:
       l1 = line.split("\t")
       lehner[l1[1]] = l1[2][0:-1]
f.close()


grdc_Merit = "results/grdc_Merit_2222.txt"
#header = "GRDC_No\tlat\tlon\tarea\tnewlat\tnewlon\tnewarea\tUPS_Indicatort\tshapearea\tLehner_shape_avail\tlehner_shapearea\tnewshape_lehner_comparison\n"
header = "No\tGRDC_No\tlat\tlon\tnewlat\tnewlon\tGRDCarea\tarea\tshapearea\tLehner_shape_avail\tlehner_shapearea\tnewshape_lehner_comparison\n"
f = open(grdc_Merit, "w")
f.write(header)
f.close()



for stationNo in range(len(grdc)):
    # ................................
    #1284
    noerror = True
    
    station = grdc[stationNo].split("\t")

    upsreal = float(station[6])
    # lon,lat  -> changed because of flw library
    coord = [ float(station[5]), float(station[4])]
    grdc_no =  station[1]

    print ("-----------------------")
    print (stationNo,grdc_no)

    try:
        no = lehner[grdc_no]
    except:
        no = grdc_no

    shapename = rootgrdcshapeold + no + ".shp"
    shapetrue = False
    if os.path.isfile(shapename):
        shapetrue = True
        s1 = gp.GeoDataFrame.from_file(shapename)
        try:
            areashape = np.sum(gpd_geographic_area_line_integral(s1) / (1000 * 1000))
        except:
            areashape = 1000 * 1482.5 * np.sum(s1.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[4]))

        # get map extent from shape file
        s1ext = s1.bounds
        xmin = int(s1ext.minx//5) * 5
        xmax = int(s1ext.maxx//5) * 5
        xanz = (xmax - xmin)//5 + 1

        ymin = int(s1ext.miny//5) * 5
        ymax = int(s1ext.maxy//5) * 5
        yanz = (ymax - ymin)//5
    else:
        # no shapefile
        sizenoarea = 5
        if upsreal > 200000:
            sizenoarea = 10

        print ("no areashape")
        areashape=0.0
        xmin = int(coord[0]//5) * 5 - sizenoarea
        xmax = xmin + sizenoarea*2 + sizenoarea
        xanz = (xmax - xmin)//5 + 1
        ymin = int(coord[1]//5) * 5 - sizenoarea
        ymax = ymin + sizenoarea*3
        yanz = (ymax - ymin)//5 + 1



    flwdiry =[]
    transform=[]

    for j in range(yanz,-1,-1):
        flwdirx = []
        for i in range(xanz):
            y = ymin + j * 5
            x = xmin + i * 5

            #upa_n30e000 make dir name
            ydir = y // 30 * 30
            xdir = x // 30 * 30

            ewdir ="e"
            nsdir ="n"
            if xdir < 0:
                xdir = abs(xdir)
                ewdir ="w"
            if ydir < 0:
                ydir = abs(ydir)
                nsdir ="s"
            # make file name

            ew ="e"
            ns ="n"
            if x < 0:
                x = abs(x)
                ew ="w"
            if y < 0:
                y = abs(y)
                ns ="s"

            #ldd = root + ns+f"{y:02d}"+ew+f"{x:03d}"+"_dir.tif"
            ldd = root + nsdir +f"{ydir:02d}" + ewdir + f"{xdir:03d}" + "/" + ns+f"{y:02d}"+ew+f"{x:03d}"+"_dir.tif"
            if os.path.isfile(ldd):

                src = rasterio.open(ldd, "r")

                flwdirx.append(src.read(1))
                transform.append(src.transform)
                latlon = src.crs.to_epsg()
                crs = src.crs
                src.close()
            else:
                empty = np.zeros((6000,6000))
                empty = empty.astype(np.uint8)
                flwdirx.append(empty)
                trans = rasterio.Affine(0.0008333333333333334, 0.0,(xmin + i * 5)-0.0004166666666666667,0.0,-0.0008333333333333334, (ymin + j * 5)+5-0.0004166666666666667)
                transform.append(trans)
                ik=1
                #[Affine(0.0008333333333333334, 0.0, -0.0004166666666666667,  0.0, -0.0008333333333333334, 39.999583333333334)]
                # n35e000_

        flwdiry.append(np.hstack((flwdirx)))

    flwdir = np.vstack((flwdiry))

    print ("loaded")
    try:
        flw = pyflwdir.from_array(flwdir, ftype="d8", transform=transform[0],check_ftype = False, latlon=latlon)
        #flw = pyflwdir.from_array(flwdir, ftype="d8", transform=src.transform,  latlon=crs.is_geographic, cache=True,)
    except Exception as e:
        noerror = False
        print ("Conversion to direction not working")
        s = str(stationNo + 1) + "\t" + grdc_no + "\t" + station[2] + "\t" + station[3] + "\t" + station[4]
        s = s + "\t" + station[5] + "\t" + station[6] + "\t" + station[7] + "\n" + str(e)
        print(s)
        s = str(stationNo) + "\t" + grdc_no + "\tConversion to direction not working\n"
        f = open(grdc_Merit, "a")
        f.write(s)
        f.close()


    #flw = flw2.load('E:/Data/GRDC_new/testEurope/flw.pkl')
    #flw.dump('E:/Data/GRDC_new/testEurope/flw.pkl')

    #coord = [12.015, 48.9475]
    #x, y = np.array([12.015]), np.array([48.9475])


    # define output locations
    #x, y = np.array([4.67916667, 7.60416667]), np.array([51.72083333, 50.3625
    #gdf_out = gp.GeoSeries(gp.points_from_xy(x, y, crs=4326))
    # delineate subbasins

    if noerror:
        print ("basin")
        try:
            subbasins = flw.basins(xy=(coord))
        except Exception as e:
            noerror = False
            print("Conversion to basin not working")
            s = str(stationNo + 1) + "\t" + grdc_no + "\t" + station[2] + "\t" + station[3] + "\t" + station[4]
            s = s + "\t" + station[5] + "\t" + station[6] + "\t" + station[7] + "\n" + str(e)
            print(s)
            s = str(stationNo) + "\t" + grdc_no + "\tConversion to basin not working\n"
            f = open(grdc_Merit, "a")
            f.write(s)
            f.close()


        if noerror:
            print ("vectorize")
            try:
                # vectorize subbasins using the vectorize convenience method from utils.py
                gdf_bas = vectorize(subbasins.astype(np.int32), 0, flw.transform, name="basin")

                a = grdc_d[grdc_no]
                a[-1] = a[-1][:-1]
                j = 0
                for h in headgrdc:
                    gdf_bas[h] = a[j]
                    j += 1
                gdf_bas['lat_merit']= coord[1]
                gdf_bas['lon_merit']= coord[0]
                gdf_bas['area_merit'] = station[7]



                #gdf_bas.head()
            except Exception as e:
                noerror = False
                print("Vectorize not working")
                s = str(stationNo+1) + "\t" + grdc_no + "\t" + station[2] + "\t" + station[3] + "\t" + station[4]
                s = s + "\t" + station[5] + "\t" + station[6] + "\t" + station[7] + "\n"+str(e)
                print(s)
                s = str(stationNo) + "\t" + grdc_no + "\tvectorize not working\n"
                f = open(grdc_Merit, "a")
                f.write(s)
                f.close()

            if noerror:
                print ("load shape")
                shapefile1 = rootshape + grdc_no + ".shp"
                gdf_bas.to_file(shapefile1)
                s2 = gp.GeoDataFrame.from_file(shapefile1)
                try:
                    areashape = np.sum(gpd_geographic_area_line_integral(s2) / (1000 * 1000))
                except:
                    areashape = 1000 * 1482.5 * np.sum(s2.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[4]))

                if shapetrue:
                    shapet = 1
                    s1 = gp.GeoDataFrame.from_file(shapename)
                    try:
                        areashape2 = np.sum(gpd_geographic_area_line_integral(s1) / (1000 * 1000))
                    except:
                        areashape2 = 1000 * 1482.5 * np.sum(s1.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[4]))

                    indicator = getBasinInd(s1,s2)
                else:
                    indicator = 0
                    areashape2 = 0
                    shapet = 0

                #header = "GRDC_No\tlat\tlon\tarea\tnewlat\tnewlon\tnewarea\tshape\tshapearea\tindicator\n"
                s = str(stationNo+1) + "\t" + grdc_no + "\t" + station[2] + "\t" + station[3] + "\t" + station[4] + "\t"+ station[5]
                s = s + "\t" + station[6] + "\t" + station[7] + "\t" + str(areashape)
                s = s + "\t" + str(shapet) + "\t" + str(areashape2) + "\t" + str(indicator) + "\n"


                print(s)
                f = open(grdc_Merit, "a")
                f.write(s)
                f.close()
                #print ("--------- -------------------------")

    #g1 = gdf_bas.simplify(0.001)
    #g1.to_file("out_simpl21.shp")


print ("Done")