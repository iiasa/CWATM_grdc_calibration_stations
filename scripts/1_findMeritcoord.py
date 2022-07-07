"""
# -------------------------------------------------------------------------
# Name:        Find MERIT coordinates
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to check and correct station location
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022

input:  grdc_2022_10577.txt   10577 station datasets >= 10km2 upstream area or no area provided
output: grdc_MERIT_1.txt: station with new location fitted to merit UPA

No: Number from 1 ...
GRDC_No: GRDC number
lat: original latitude from GRDC metafile
lon: original longituted from GRDC metafile
newlat: corrected latitude based on MERIT UPA dataset
newlon: corrected longitute based on MERIT UPA dataset
area; provided basin area from GRDC metafile
newarea: basin area based on MERIT UPA dataset
UPS_Indicator:  min error in % from MERIT UPA to provided basin area
dist_Indicator: distance to original pour point in [unit:100m]
Indicator:  ranking criteria: UPS_Indicator + 2 x dist_indicator

uses files:
grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID

# ----------------------------------------------------------------------
"""
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes

from osgeo import gdal
from osgeo import osr
from osgeo import gdalconst

import sys
import os

# ----------------------------------------------------------------

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


#----------------------------------------------
# INPUT
# MERIT Yamazaki et al 2019 - upstream area in km2
root =  "P:/watproject/Datasets/MERIT_yamazaki/upa_"
root2 =  "P:/watproject/Datasets/MERIT_yamazaki/amodel2/"
# GRDC stations >= 10km2
grdc_stations = root2 + "grdc_2022_10577.txt"
shapefolder = "P:/watproject/Datasets/MERIT_yamazaki/grdc_station_upstream/grdc_basins_smoothed_md_no_"

#OUTPUT
grdc_Merit = "results/grdc_Merit_1.txt"

# --------------------------------------------------------------------------------
# cell size: 3 arcsec
cell = 0.000833333333333333333
# 1/cell
invcell = 1200
# search range in cells: 55 = around 5km
rangexy=55

f = open(grdc_stations, "r")
grdc = f.readlines()[1:]
f.close()

#grdc_Merit = "results/grdc_Merit_1.txt"
header = "No\tGRDC_No\tlat\tlon\tnewlat\tnewlon\tGRDCarea\tarea\tUPS_Indicator\tdist_Indicator\tIndicator\n"
f = open(grdc_Merit, "w")
f.write(header)
f.close()


# translation file new GRDC ID - old GRDC ID
lehnerfile = root2 +"grdc2022_lehner2012.txt"
lehner = {}
with open(lehnerfile) as f:
    for line in f:
       l1 = line.split("\t")
       lehner[l1[1]] = l1[2][0:-1]
f.close()



# -----------------------------------
jjj = 0
for stationNo in range(len(grdc)):

   station = grdc[stationNo].split("\t")
   upsreal = float(station[9])
   # lat lon
   coord = [ float(station[7]), float(station[8])]
   grdc_no =  station[1]

   if upsreal<0:
       try:
            no = lehner[grdc_no]
       except:
            no = grdc_no
       shapename = shapefolder + no + ".shp"
       if os.path.isfile(shapename):
            s1 = gpd.GeoDataFrame.from_file(shapename)
            try:
                upsreal = np.sum(gpd_geographic_area_line_integral(s1) / (1000 * 1000))
            except:
                upsreal = 1000 * 1482.5 * np.sum(s1.area) * 0.092593 * 0.092593 * np.cos(
                    np.pi / 180.0 * float(station[7]))


   jjj = jjj + 1
   if jjj > -1:

    # calculate which ups is used - in case the coordinate are at the edge other ups are loaded
    xmin = int(coord[1] // 5) * 5
    ymin = int(coord[0] // 5) * 5

    x1 = coord[1] % 5
    if x1 <  0.15:
        xx = [-1,0]
    elif x1> 4.85:
        xx = [0,1]
    else:
        xx = [0]

    y1 = coord[0] % 5
    if y1 <  0.15:
        yy = [0,-1]
    elif y1> 4.85:
        yy = [1,0]
    else:
        yy = [0]


    upsy =[]
    transform=[]

    for j in yy:
        upsx = []
        for i in xx:
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

            upsname = root + nsdir +f"{ydir:02d}" + ewdir + f"{xdir:03d}" + "/" + ns+f"{y:02d}"+ew+f"{x:03d}"+"_upa.tif"
            src = rasterio.open(upsname, "r")
            upsx.append(src.read(1))
            transform.append(src.transform)
            latlon = src.crs.to_epsg()
            crs = src.crs
            src.close()

        upsy.append(np.hstack((upsx)))

    ups = np.vstack((upsy))
    col = ups.shape[1]
    row = ups.shape[0]
    top = transform[0][5]
    left = transform[0][2]

    col1 = int((coord[1] - left) * invcell)
    row1 = int((top -coord[0]) * invcell)
    ups1 = ups[row1,col1]

    rangexy = 55
    upsups = np.zeros((rangexy*2+1,rangexy*2+1))
    ind = np.zeros((rangexy*2+1,rangexy*2+1))
    upsind = np.zeros((rangexy*2+1,rangexy*2+1))
    diffind = np.zeros((rangexy*2+1,rangexy*2+1))

    colcol = np.arange(col1-rangexy,col1+rangexy+1)
    rowrow = np.arange(row1-rangexy,row1+rangexy+1)

    j =0
    for y in rowrow:
        i = 0
        for x in colcol:
            upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
            upsups[j,i]= ups[y,x]
            diff = np.sqrt((rangexy-i)**2+(rangexy-j)**2)*0.9
            diffind[j,i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
            # if upsind> 50 diff gets a penalty
            if upsind[j, i]>50:
                diffind[j, i] = diffind[j, i] + 500
            ind[j,i] = upsind[j, i] + 2 * diffind[j,i]

            i = i +1
        j = j + 1

    minxy = np.where(ind==np.min(ind))
    y=minxy[0][0]
    x=minxy[1][0]
    j = rowrow[y]
    i = colcol[x]

    yy=coord[0]+(rangexy-y)*cell
    xx=coord[1]-(rangexy-x)*cell
    ups2 = ups[j,i]


    #------------------------------------------------------
    # if still big error increae range
    if ind[y,x] > 50:
        print ("increase range")
        rangexy = 101
        upsups = np.zeros((rangexy*2+1,rangexy*2+1))
        ind = np.zeros((rangexy*2+1,rangexy*2+1))
        upsind = np.zeros((rangexy*2+1,rangexy*2+1))
        diffind = np.zeros((rangexy*2+1,rangexy*2+1))

        colcol = np.arange(col1-rangexy,col1+rangexy+1)
        rowrow = np.arange(row1-rangexy,row1+rangexy+1)

        j =0
        for y in rowrow:
            i = 0
            for x in colcol:
                upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
                upsups[j,i]= ups[y,x]
                diff = np.sqrt((rangexy-i)**2+(rangexy-j)**2)*0.9
                diffind[j,i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
                # if upsind> 50 diff gets a penalty
                if upsind[j, i] > 50:
                    diffind[j, i] = diffind[j, i] + 500
                ind[j,i] = upsind[j, i] + 0.5 * diffind[j,i]

                i = i +1
            j = j + 1

        minxy = np.where(ind==np.min(ind))
        y=minxy[0][0]
        x=minxy[1][0]
        j = rowrow[y]
        i = colcol[x]

        yy=coord[0]+(rangexy-y)*cell
        xx=coord[1]-(rangexy-x)*cell
        ups2 = ups[j,i]
    #-------------------------------------------------

    # ------------------------------------------------------
    # if still big error increae range
        if ind[y, x] > 80:
            print("increase range2")
            rangexy = 151
            upsups = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            ind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            upsind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            diffind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))

            colcol = np.arange(col1 - rangexy, col1 + rangexy + 1)
            rowrow = np.arange(row1 - rangexy, row1 + rangexy + 1)

            j = 0
            for y in rowrow:
                i = 0
                for x in colcol:
                    upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
                    upsups[j, i] = ups[y, x]
                    diff = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.9
                    diffind[j, i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
                    # if upsind> 50 diff gets a penalty
                    if upsind[j, i] > 50:
                        diffind[j, i] = diffind[j, i] + 1000
                    ind[j, i] = upsind[j, i] + 0.25 * diffind[j, i]

                    i = i + 1
                j = j + 1

            minxy = np.where(ind == np.min(ind))
            y = minxy[0][0]
            x = minxy[1][0]
            j = rowrow[y]
            i = colcol[x]

            yy = coord[0] + (rangexy - y) * cell
            xx = coord[1] - (rangexy - x) * cell
            ups2 = ups[j, i]
    # -------------------------------------------------

    #-------------------------------------------------

    s= str(jjj)+"\t"
    s = s + grdc_no + "\t"
    #print(i,j,x,y,ups2,ind[y,x])
    #print(upsind[y,x],diffind[y,x],ind[y,x])
    s = s + f"{yy:.3f}" + "\t" + f"{xx:.3f}" + "\t" + f"{ups2:.0f}"+ "\t"
    #,f"{ups:.1f}",f"{ind[y,x]:.2f}")
    s = s + f"{coord[0]:.3f}"+ "\t" +f"{coord[1]:.3f}" + "\t"+f"{upsreal:.0f}"
    print (s)
    #print("-------------------")

    # 2418 1262 12.015 48.9475 26429.025 0.07174307839912242
    # upsreal = 26448 coord = [12.015, 48.9475]
    #header = "GRDC_No\tlat\tlon\tarea\tnewlat\tnewlon\tnewarea\tUPS_Indicator\tdist_Indicator\tIndicator"
    s = str(stationNo+1) +"\t"+ grdc_no +"\t"+str(coord[0]) + "\t" + str(coord[1])
    s = s + "\t" + str(yy)+"\t"+str(xx) +"\t"+ str(upsreal) + "\t"+ str(ups2)
    s = s + "\t" + str(upsind[y,x]) +"\t"+str(diffind[y,x])+"\t"+str(ind[y,x]) + "\n"
    f = open(grdc_Merit, "a")
    f.write(s)
    f.close()




    ii =1
