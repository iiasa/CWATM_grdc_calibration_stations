# -------------------------------------------------------------------------
# Name:        shapefile in low-res
# Purpose:     creates shapefiles and a list in lower resolution (5 or 30 arcmin)
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022
# ----------------------------------------------------------------------


import geopandas as gp
import numpy as np
import rasterio
from rasterio import features
import pyflwdir

import sys
import os
import os.path
import warnings
#-----------------------------------------------



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
#----------------------------------------------

warnings.filterwarnings("ignore")

root =  "P:/watproject/Datasets/MERIT_yamazaki/"

# INPUT DATA
rootshape = root + "/ashape2_merit/grdc_basin_merit_"
maketxt = "results/grdc_shape_allend_1.txt"

# resolution 30arcmin or 5arcmin
reso30 = True

if reso30:
    # input 30 arcmin
    rootshape2 = root + "/ashape30min/grdc_basin_30min_"
    lddname = "P:/watmodel/CWATM/cwatm_input_isimip3/routing/ldd_new.tif"
    upsname = "P:/watmodel/CWATM/cwatm_input_isimip3/routing/ups_new.tif"
    # output 30min
    grdc_Merit = "results/basins_30min.txt"
    grdc_error = "results/basins_30min_error.txt"
else:
    # input 5 arcmin
    rootshape2 = root + "/ashape5min/grdc_basin_5min_"
    lddname = "P:/watmodel/CWATM/CWAT_input_5min/processing/routing_eilander/kinematic/ldd.tif"
    upsname = "P:/watmodel/CWATM/CWAT_input_5min/processing/routing_eilander/kinematic/ups.tif"
    # output
    grdc_Merit = "results/basins_5min.txt"
    grdc_error = "results/basins_5min_error.txt"


#maketxt = "grdc_shape_allend_1.txt"
f = open(maketxt, "r")
makeshp = f.readlines()[1:]
f.close()

#grdc_Merit = "basins_5min.txt"
header = "No\tGRDC_No\tsimilarity\tareaGRDC\tarea\tlat\tlon\tarea30min\tlat30min\tlon30move\tlat30move\tlon30min\tindloc\tupsloc\tind\tindups\tups1\tindshape\tshapeloc\tind\tindups\tups2\tindshape\n"
f = open(grdc_Merit, "w")
f.write(header)
f.close()

#grdc_error = "basins_5min_error.txt"
header = "No\tGRDC_No\tsimilarity\tareaGRDC\tarea\tlat\tlon\tarea30min\tlat30min\tlon30move\tlat30move\tlon30min\tindloc\tupsloc\tind\tindups\tups1\tindshape\tshapeloc\tind\tindups\tups2\tindshape\n"
f = open(grdc_error, "w")
f.write(header)
f.close()

src = rasterio.open(upsname)
ups = src.read(1)

src = rasterio.open(lddname, "r")
flwdir = src.read(1)
transform1 = src.transform
latlon1 = src.crs.to_epsg()
crs = src.crs
src.close()
flw = pyflwdir.from_array(flwdir, ftype="ldd", transform=transform1, check_ftype=True, latlon=latlon1)


#---------------------------------------
if reso30:
    rangexy = np.arange(-2 * 0.5, 3 * 0.5, 0.5)
    # upstream area in km2 , equal or bigger is used  for this resolution
    threshold = 9000
    reso = "30min"
    mult = 2.0
else:
    rangexy = np.arange(-2*0.08333333333,3*0.08333333333,0.08333333333)
    threshold = 1000
    reso ="5min"
    mult = 12.0

# +++++++++++++++++++++++++++++++++
# run through als stations

for stationNo in range(len(makeshp)):

    station = makeshp[stationNo].split("\t")
    upsreal = float(station[6])
    upsmodel = float(station[7])

    # x,y
    coord = [ float(station[4]), float(station[5])]
    grdc_no =  station[1]
    #if grdc_no =="6343900":

    if upsreal< 0:
        upsreal = float(station[7])

    # only station bigger equal than 1000km2 (5min) 9000km2 (30min)
    if (upsreal >= threshold) or (upsmodel >= threshold):


        shapefile1 = rootshape +"smooth_"+ grdc_no + ".shp"
        if os.path.isfile(shapefile1):
            s2 = gp.GeoDataFrame.from_file(shapefile1)

            j = 0
            ind1 = []
            ind2 = []
            ind = []
            ups25 =[]
            for y in rangexy:
                for x in rangexy:
                    #print (j)
                    xx = coord[1] + x
                    yy = coord[0] + y
                    if (xx>-180) and (xx<180):
                        c = [xx,yy]
                        subbasins = flw.basins(xy=(c))
                        basin1km = vectorize(subbasins.astype(np.int32), 0, flw.transform, name="basin")
                        shapefile2 = "P:/watproject/Datasets/MERIT_yamazaki/amodel/test/basin_"+reso+"_" + grdc_no + "_" + str(j) + ".shp"
                        #basin1km.to_file(shapefile2)

                        # calculate union and intersection of shapes
                        p_inter = gp.overlay(s2, basin1km, how='intersection')
                        p_union = gp.overlay(s2, basin1km, how='union')
                        pint_area = p_inter.area.sum()
                        puni_area = p_union.area.sum()
                        indshape = (pint_area / puni_area)
                        ind1.append(indshape)

                        # get upstream area of coarse grid
                        col = int((xx + 180) * mult)
                        row = int((90 - yy) * mult)
                        upsvalue =  ups[row,col]

                        if upsreal == 0 or upsvalue == 0:
                            indups = 0
                        else:
                            if upsreal < upsvalue:
                                indups = upsreal / upsvalue
                            else:
                                indups = upsvalue / upsreal
                        ind2.append(indups)
                        ups25.append(upsvalue)

                        # calculate lenght to origin (0,0)
                        ind.append(np.sqrt((1-indups)**2+(1-indshape)**2))

                        j += 1

            #---------------------------------
            maxups = np.max(ind2)
            maxshape = np.max(ind1)
            shapemax = np.where(ind1 == maxshape)[0][0]
            upsmax = np.where(ind2 == maxups)[0][0]


            minind = np.min(ind)
            indmin = np.where(ind==minind)[0][0]
            # get location in a 5x5 matrix
            y= indmin // 5
            x = indmin % 5
            yy = coord[0] + rangexy[y]
            xx = coord[1] + rangexy[x]

            subbasins = flw.basins(xy=([xx,yy]))
            basin1km = vectorize(subbasins.astype(np.int32), 0, flw.transform, name="basin")
            shapefile2 = rootshape2 + "basin_" + grdc_no + ".shp"
            #basin1km.to_file(shapefile2)

            p_inter = gp.overlay(s2, basin1km, how='intersection')
            p_union = gp.overlay(s2, basin1km, how='union')
            pint_area = p_inter.area.sum()
            puni_area = p_union.area.sum()
            indicator = pint_area / puni_area

            try:
                #shapefile2 = "P:/watproject/Datasets/MERIT_yamazaki/amodel/test/b_" + grdc_no + "_union.shp"
                shapefile2 = "P:/watproject/Datasets/MERIT_yamazaki/amodel/test/b__union.shp"
                p_union.to_file(shapefile2)
                shapefile2 = "P:/watproject/Datasets/MERIT_yamazaki/amodel/test/b__inter.shp"
                p_inter.to_file(shapefile2)
            except:
                f = open(grdc_error, "a")
                s = str(stationNo) + "\t" + grdc_no + "\t" + f"{indicator:6.3f}" + "\t" + f"{upsreal:7.0f}"
                s = s + "\t" + f"{upsmodel:7.0f}" + "\t" + f"{coord[0]:7.4f}" + "\t" + f"{coord[1]:6.4f}" + "\n"
                f.write(s)
                f.close()


            col = int((xx + 180) * mult)
            row = int((90 - yy) * mult)
            upsvalue = ups[row, col]

            if reso30:
                locx = xx // 0.5 * 0.5 + 0.25
                locy = yy // 0.5 * 0.5 + 0.25
            else:
                locx = xx//0.083333333333*0.083333333333+0.083333333333/2
                locy = yy//0.083333333333*0.083333333333+0.083333333333/2

            s = str(stationNo) +  "\t" + grdc_no +  "\t" +f"{indicator:6.3f}"+ "\t" +f"{upsreal:7.0f}"
            s = s + "\t" +f"{upsmodel:7.0f}" + "\t" +f"{coord[0]:7.4f}" + "\t" + f"{coord[1]:6.4f}"
            s = s + "\t" +f"{upsvalue:9.0f}"  + "\t" +f"{yy:8.3f}" + "\t" + f"{xx:6.3f}" + "\t" + f"{locy:7.2f}" + "\t" + f"{locx:5.2f}"
            s = s + "\t" + str(indmin)
            # loc of area best, indicator (min best), area (max best), shape similarity (max best)
            s = s + "\t" +  str(upsmax) + "\t" + f"{ind[upsmax]:6.3f}" + "\t" + f"{ind2[upsmax]:6.3f}"  + "\t" + f"{ups25[upsmax]:7.0f}" + "\t" + f"{ind1[upsmax]:6.3f}"
            s = s + "\t" +  str(shapemax) + "\t" + f"{ind[shapemax]:6.3f}" + "\t" + f"{ind2[shapemax]:6.3f}"  + "\t" + f"{ups25[shapemax]:7.0f}" + "\t" + f"{ind1[shapemax]:6.3f}"
            print (s)
            #print(stationNo,indicator,upsvalue,upsreal,yy,xx,locy,locx)
            s = s + "\n"
            f = open(grdc_Merit, "a")
            f.write(s)
            f.close()

