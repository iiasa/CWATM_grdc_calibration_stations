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


GRDC_No: GRDC number
lat: original latitude from GRDC metafile
lon: original longituted from GRDC metafile
area; provided basin area from GRDC metafile
newlat: corrected latitude based on MERIT UPA dataset
newlon: corrected longitute based on MERIT UPA dataset
newarea: basin area based on MERIT UPA dataset
UPS_Indicator:  min error in % from MERIT UPA to provided basin area
shapearea: area of the shape calculated with geopandas directly from shape
Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar

# uses
# grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
# grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
# Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
"""

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

from shapely.geometry.polygon import orient
from shapely.geometry import Polygon


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

warnings.filterwarnings("ignore")

root =  "P:/watproject/Datasets/MERIT_yamazaki/upa_"
root2 =  "P:\\watproject\\Datasets\\MERIT_yamazaki\\upa_"

rootshape = "P:/watproject/Datasets/MERIT_yamazaki/ashape_merit/grdc_basin_merit_"
rootshape2 = "P:/watproject/Datasets/MERIT_yamazaki/ashape2x_merit/grdc_basin_merit_"
rootgrdcshape = "P:/watproject/Datasets/MERIT_yamazaki/grdc_station_upstream/"

grdc_stations = "makeshape_all10349.txt"
#----------------------------------------------------------


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
f.close()

# results
grdc_Merit = "grdc_shape_allend2.txt"
header = "No\t GRDC_No\tlat\tlon\tnewlat\tnewlon\tGRDCarea\tarea\tshapearea\toldsha\toldshapearea\tindicator\n"
f = open(grdc_Merit, "w")
f.write(header)
f.close()

for stationNo in range(len(grdc)):

    noerror = True
    station = grdc[stationNo].split("\t")
    grdc_no =  station[1]
    ii =1

    upsreal = float(station[4])
    # lon, lat
    coord = [ float(station[6]), float(station[5])]
    x = int(coord[0]//5) * 5
    y = int(coord[1]//5) * 5

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

    upsname = root2 + nsdir + f"{ydir:02d}" + ewdir + f"{xdir:03d}" + "\\" + ns + f"{y:02d}" + ew + f"{x:03d}" + "_upa.tif"

    src = rasterio.open(upsname, "r")
    ups=src.read(1)
    transform=src.transform
    latlon = src.crs.to_epsg()
    crs = src.crs
    src.close()
    #print ("loaded")

    top = transform[5]
    left = transform[2]
    invcell = 1200.0

    col = int((coord[0] - left) * invcell)
    row = int((top -coord[1]) * invcell)
    ups = ups[row,col]

    indicator = 0
    areashape1 = 0
    areashape2 = 0

    shapefile1 = rootshape + grdc_no + ".shp"
    if os.path.isfile(shapefile1):
        gdf_bas = gp.GeoDataFrame.from_file(shapefile1)
        try:
            areashape = np.sum(gpd_geographic_area_line_integral(gdf_bas) / (1000 * 1000))
        except:
            areashape = 1000 * 1482.5 * np.sum(gdf_bas.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[6]))

    shap = station[11]
    if shap[0:1] == "1":
        shapename = rootgrdcshape + "grdc_basins_smoothed_md_no_"+ grdc_no +".shp"
        s1 = gp.GeoDataFrame.from_file(shapename)
        try:
            areashape1 = np.sum(gpd_geographic_area_line_integral(s1) / (1000 * 1000))
        except:
            areashape1 = 1000 * 1482.5 * np.sum(s1.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[6]))
        indicator = getBasinInd(s1, gdf_bas)

    quot = -999
    quot2 = 0
    save1 = False
    if upsreal < 0:
        save1 = True
    else:
        if upsreal>ups:
            quot = ups/upsreal
        else:
            quot = upsreal/ups
        #if quot >=0.75:
        #    save = True
        #if quot < 0.75:
        if quot>=0.4:
            save1 = True

    save2 = False
    # shapefile
    #shapefile1 = rootshape + grdc_no + ".shp"
    #gdf_bas = gp.GeoDataFrame.from_file(shapefile1)
    #try:
    #    areashape = np.sum(gpd_geographic_area_line_integral(gdf_bas) / (1000 * 1000))
    #except:
    #    areashape = 1000 * 1482.5 * np.sum(gdf_bas.area) * 0.092593 * 0.092593 * np.cos(np.pi / 180.0 * float(station[6]))


    if ups > 3:
        if areashape >3:
            if ups > areashape:
                quot2 = areashape/ups
            else:
                quot2 = ups/areashape
            if quot2>=0.8:
                save2 = True



    #s = str(stationNo) + "\t" + grdc_no + "\t" + f"{float(station[2]):8.4f}" + "\t" + f"{float(station[3]):8.4f}"
    #s = s + "\t" + f"{coord[1]:9.6f}" + "\t" + f"{coord[0]:9.6f}"
    #s = s + "\t" + f"{upsreal:10.0f}" + "\t" + f"{ups:10.0f}"
    save = False
    if save1:
        if save2:
            save = True


    if save:
        # shapefile
        #shapefile1 = rootshape + grdc_no + ".shp"
        #gdf_bas = gp.GeoDataFrame.from_file(shapefile1)
        g1 = gdf_bas.simplify(0.001)

        shapefile1 = rootshape2 +"smooth_"+ grdc_no + ".shp"
        g1.to_file(shapefile1)

        gdf_bas = gp.GeoDataFrame.from_file(shapefile1)
        try:
            a = grdc_d[grdc_no]
            a[-1] = a[-1][:-1]
            j = 0
            for h in headgrdc:
                gdf_bas[h] = a[j]
                j += 1
        except:
            iii =1
        gdf_bas['lat_merit'] = coord[1]
        gdf_bas['lon_merit'] = coord[0]
        gdf_bas['area_merit'] = ups

        shapefile1 = rootshape2 +"smooth_"+ grdc_no + ".shp"
        gdf_bas.to_file(shapefile1)
        #print(s)


    else:

        # delete smooth
        ii =1
        #s = s + "\tNot saved"
    #print(s)

    if shap[0:1] =="1":
        s1 = "TRUE"
    else:
        s1 = "FALS"
    if save:
        s = str(stationNo) +"\t"+ grdc_no +"\t"+f"{float(station[2]):8.4f}" + "\t" + f"{float(station[3]):8.4f}"
        s = s + "\t" +  f"{coord[1]:9.6f}" + "\t" +  f"{coord[0]:9.6f}"
        s = s + "\t" + f"{upsreal:10.0f}" + "\t" +  f"{ups:10.0f}" + "\t" + f"{areashape:10.0f}"
        s = s + "\t" + s1  +"\t" + f"{areashape1:10.0f}" +"\t" + f"{indicator:.3f}"



        print(s)
        s = s + "\t" + upsname + "\n"
        f = open(grdc_Merit, "a")
        f.write(s)
        f.close()
    else:
        s = str(stationNo) +"\t"+ grdc_no +"\t"+f"{float(station[2]):8.4f}" + "\t" + f"{float(station[3]):8.4f}"
        s = s + "\t" +  f"{coord[1]:9.6f}" + "\t" +  f"{coord[0]:9.6f}"
        s = s + "\t" + f"{upsreal:10.0f}" + "\t" +  f"{ups:10.0f}" + "\t" + f"{areashape:10.0f}"
        s = s + "\t" + s1 + "\t" + f"{areashape1:10.0f}" + "\t" + f"{indicator:.3f}"
        s = s + "\t" + f"{quot:9.3f}" + "\t" + f"{quot2:9.3f}"

        print(s)

        shapefile1 = rootshape2 + "smooth_" + grdc_no + ".shp"
        try:
            os.remove(shapefile1)
        except:
            ii = 1



print ("Done")