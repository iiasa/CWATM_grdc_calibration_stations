"""
# -------------------------------------------------------------------------
# Name:        Create list of stations to calibrate or validate
# Purpose:     using the watershed boundary shaepfile on 30min (or 5min) and the station metadata
#              to dismiss stations which are not qualify
#              and to select vetween calibration and validation stations
#              Validation station are to close to calibration stations, but have a lower score
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022
# ----------------------------------------------------------------------

# input:  basins_30min.txt or basins_5min.txt: station with new location fitted to 30 arcmin or 5 arcmin
# output: results/grdc_station_sel30min.txt

- if a station can be used for calibration it gets a C as indicator
- if a station passes the selection but is to close to another station which is used for calibration it gets a V
- if a station does not pass the selection criteria it gets a N

Selection criteria:
- similarityMin: Minimum similarity between 3 arcsec basin and low-res basin: here: 0.7 (1 = identicaly, 0 = nothing in common)
- areaaccorMin: Minimim area error: area provided from GRDC and area from low-res shapefile should have at least (e.g 0.4)
  area = 0.3 x areaGRDC is dismissed , area = 0.5 x areaGRDC is used
- t_yrsMin = 5:  5 minimum years of record
t_endMin = 1985: last year should be 1985 or younger for daily timeseries
d_endMin = 1985  last year for monthly timeseries

Scoring points:
if 2 (or more) station are too close to each other the one is taken with the most scoring points
- simiScore = 2  # every 2% scoring +1
- areaScore = 2   # every 2% scoring +1
- tyrsScore = 5  # every 5 years a score
- tyrsScoreMax = 10 # max score for end years
- tendScore = 3  # every 3 years a point
- tendScoremax = 100 # maximum score from this can only be 100
- dayscore = 5  # points if daily timeseries
- MissingScore = 5  # neg. points for missing values

grdc_station_sel30min.txt
-------------------------
No: Number from 0 ...
GRDC_No: GRDC number
latGRDC: original latitude from GRDC metafile
lonGRDC: original longituted from GRDC metafil
latcorrected: coorected latitude on high-res
loncorrected: corrected longitude on high-res
lat30min: latitude moved to grid cell centre
lon30min: longitude moved to grid cell centre
areaGRDC: area provided by GRDC
area: area from high-res UPA MERIT at pour point
area30min : area on pour point of best fitting 30arcmin grid cell
similarity: similarity of 30min low-res shape with high-res shape from 3arcsec
area_accordance> comparison of area30min with areaGRDC
calibrate_orNot> can be use for calibration (C), validation (V) or is dismissed (N)
scoring: scoring based on the scoring criteria
grdc_scored_better: grdc number with scored better
better_score: better score of the grdc no which is used for calibration

# uses
# grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
# Shapefiles with low resolution from: ashape30min/grdc_basin_30min_basin_

"""

import geopandas as gp
import numpy as np
import sys
import os
import os.path
import warnings



def getBasinInd(grdc_shape,newshape):

    p_inter = gp.overlay(grdc_shape, newshape, how='intersection')
    p_union = gp.overlay(grdc_shape, newshape, how='union')

    pint_area = p_inter.area.sum()
    puni_area = p_union.area.sum()
    indicator = pint_area / puni_area

    # s1 = grdc_shape.symmetric_difference(p2, align=False)
    # sym_min = s1.area

    return indicator
# ------------------------
def getvalue(inv):
    if inv == "":
        value = 0
    else:
        value = float(inv)
    return value
#----------------------------------------------------

warnings.filterwarnings("ignore")

root ="P:/watproject/Datasets/MERIT_yamazaki/"

# INPUT DATA
grdcallfile = "grdc_2022_10701.txt"

# resolution 30arcmin or 5arcmin
reso30 = True

if reso30:
    # input 30 arcmin
    shapedir = "P:/watproject/Datasets/MERIT_yamazaki/ashape30min/grdc_basin_30min_basin_"
    grdc_stations = "results/basins_30min.txt"
    # output
    grdc_Station = "results/grdc_station_sel30min.txt"
    header = "No\tGRDC_No\tlatGRDC\tlonGRDC\tlatcorrected\tloncorrected\tlat30min\tlon30min\tareaGRDC\tarea\tarea30min\tsimilarity\tarea_accordance\tcalibrate_orNot\tscoring\tgrdc_scored_better\tbetter_score\n"


else:
    # input 5 arcmin
    shapedir = "P:/watproject/Datasets/MERIT_yamazaki/ashape5min/grdc_basin_5min_basin_"
    grdc_stations = "results/basins_5min.txt"
    # output
    grdc_Station = "results/grdc_station_sel5min.txt"
    header = "No\tGRDC_No\tlatGRDC\tlonGRDC\tlatcorrected\tloncorrected\tlat5min\tlon5min\tareaGRDC\tarea\tarea5min\tsimilarity\tarea_accordance\tcalibrate_orNot\tscoring\tgrdc_scored_better\tbetter_score\n"

#----------------------------------------------------------

f = open(grdc_stations, "r")
grdc = f.readlines()[2:]
f.close()

#grdcallfile = "grdc_2022_10701.txt"
grdc_d = {}
with open(grdcallfile) as f:
    for line in f:
       l1 = line.split("\t")
       grdc_d[l1[0]] = l1
headgrdc = grdc_d["grdc_no"]
headgrdc[-1] = headgrdc[-1][0:-1]
f.close()


f = open(grdc_Station, "w")
f.write(header)
f.close()


# ------------------------------------
# ranking
# deselect:
if reso30:
    # 30min
    similarityMin = 0.7
else:
    # 5min
    similarityMin = 0.8


areaaccorMin = 0.4
t_yrsMin = 5    # 5 minimum years of record
t_endMin = 1985  # last year should be 1985 or younger
d_endMin = 1985

# scoring points
simiScore = 2  # every 2% scoring +1
areaScore = 2
tyrsScore = 5  # every 5 years a score
tyrsScoreMax = 10 # max score for end years
tendScore = 3  # every 3 years a point
tendScoremax = 100
dayscore = 5  # points if daily timeseries

MissingScore = 5  # pneg. points for missing values

# ------------------------------------
grdc_area ={}
region = {}
grdcselect = {}
scoreall = {}

lengrdc = len(grdc)
#lengrdc = 100
for stationNo in range(lengrdc):
    station = grdc[stationNo].split("\t")
    grdc_no =  station[1]
    grdcOri = grdc_d[grdc_no]

    upsGRDC = float(station[3])
    upsreal = float(station[7])
    coord = [ float(station[8]), float(station[9])]       # lat, lon
    similarity = float(station[2])
    if upsGRDC > upsreal:
        areaAccor = upsreal / upsGRDC
    else:
        areaAccor = upsGRDC / upsreal
    if areaAccor < 0: areaAccor = 1.0

    t_yrs = float(grdcOri[20])
    t_end = getvalue(grdcOri[19])
    d_end = getvalue(grdcOri[11])
    dmiss = (100-getvalue(grdcOri[13]))/100 * getvalue(grdcOri[12])
    mmiss = (100-getvalue(grdcOri[17]))/100 * getvalue(grdcOri[16])
    miss = max(dmiss,mmiss)

    select = True
    if similarity < similarityMin: select = False
    if areaAccor  < areaaccorMin: select = False
    if t_yrs  < t_yrsMin: select = False
    if t_end < t_endMin: select = False


    # score
    score = []
    score.append(int(100 * (similarity - similarityMin) / simiScore))
    score.append(int(100 * (areaAccor - areaaccorMin) / simiScore))
    score.append(min(int((t_yrs - t_yrsMin) / tyrsScore), tyrsScoreMax))
    score.append(int((t_end - t_endMin) / tendScore))
    if d_end < d_endMin:
        score.append(0)
    else:
        score.append(dayscore)

    sumscore = sum(score)
    scoreall[grdc_no] = [sumscore, score]

    if select:
        entry = []
        entry.append([grdc_no,upsGRDC,float(grdcOri[6]),float(grdcOri[7]),upsreal,coord[0],coord[1],similarity,areaAccor])
        entry.append([t_yrs,t_end ,d_end,miss,sumscore,score])
        grdcselect[grdc_no] = entry

        # area to sort afterwards
        grdc_area[grdc_no] = upsreal

        # region to preselect the basins
        grdc3 = grdc_no[0:3]
        if grdc3 in region:
            region[grdc3].append(grdc_no)
        else:
            region[grdc3] = [grdc_no]

        ii=1

    ii =1
    grdc_sort = sorted(grdc_area.items(), key=lambda x: x[1], reverse=False)

print("preselect done")

validate = {}

for stationNo in range(len(grdc_sort)):
    grdc_no = grdc_sort[stationNo][0]
    if grdc_no in grdc_area:
        shapefile = shapedir + grdc_no + ".shp"
        s1 = gp.GeoDataFrame.from_file(shapefile)
        grdc3 = grdc_no[0:3]
        upsreal = grdcselect[grdc_no][0][4]


        for g1 in region[grdc3]:
            test = True
            if g1 == grdc_no:
                test = False
            upstest = grdcselect[g1][0][4]
            if upsreal < (upstest*0.75):
                test = False
            if upsreal > (upstest*1.25):
                test = False

            if g1 in validate:
                test = False
            if grdc_no in validate:
                test = False

            if test:
                shapefile = shapedir + g1 + ".shp"
                s2 = gp.GeoDataFrame.from_file(shapefile)
                indicator = getBasinInd(s1, s2)
                #print (grdc_no,g1,indicator)
                if indicator >= 0.95:
                    print(grdc_no, g1, indicator,upsreal,upstest)
                    print ("--",grdc_no,grdcselect[grdc_no][1][4],grdcselect[grdc_no][1][5] )
                    print("--", g1, grdcselect[g1][1][4], grdcselect[g1][1][5])
                    if grdcselect[g1][1][4] > grdcselect[grdc_no][1][4]:
                        validate[grdc_no]= [g1,indicator,upsreal,upstest,grdcselect[grdc_no][1][4],grdcselect[g1][1][4]]
                    else:
                        validate[g1]= [grdc_no,indicator,upstest,upsreal,grdcselect[g1][1][4],grdcselect[grdc_no][1][4]]
                    ii =1
            ii=2

print ("calib - valid done")
nocal = 0

for stationNo in range(lengrdc):

    station = grdc[stationNo].split("\t")
    grdc_no =  station[1]
    grdcOri = grdc_d[grdc_no]

    upsGRDC = float(station[3])  # original area provided
    upsreal = float(station[7])  # low-res area
    coord = [ float(station[8]), float(station[9])]       # lat, lon
    similarity = float(station[2])
    if upsGRDC > upsreal:
        areaAccor = upsreal / upsGRDC
    else:
        areaAccor = upsGRDC / upsreal
    if areaAccor < 0: areaAccor = 1.0

    s = str(stationNo) + "\t" + grdc_no + "\t" + f"{float(grdcOri[6]):8.4f}" + "\t" + f"{float(grdcOri[7]):8.4f}"
    s = s + "\t" + f"{float(station[5]):8.4f}" + "\t" + f"{float(station[6]):8.4f}" + "\t" + f"{float(station[10]):8.4f}"+ "\t" + f"{float(station[11]):8.4f}"
    s = s + "\t" + f"{upsGRDC:10.0f}" + "\t" + f"{float(station[4]):10.0f}" + "\t" + f"{float(station[7]):9.6f}"
    s = s + "\t" + f"{similarity:6.3f}" + "\t" + f"{areaAccor:6.3f}"



    #s = str(stationNo) + "\t" + grdc_no + "\t" + f"{upsGRDC:10.0f}" + "\t" + f"{float(grdcOri[6]):8.4f}" + "\t" + f"{float(grdcOri[7]):8.4f}"
    #s = s + "\t" + f"{upsreal:10.0f}" + "\t" + f"{coord[1]:9.6f}" + "\t" + f"{coord[0]:9.6f}"
    #s = s + "\t" + f"{similarity:6.3f}" + "\t" + f"{areaAccor:6.3f}"

    calib = "N"
    if grdc_no in grdcselect:
        calib = "C"
        nocal = nocal +1
    if grdc_no in validate:
        # station which are similar to other stations
        calib = "V"
    s = s + "\t" + calib
    if grdc_no in grdcselect:
        # + sum scoring - higher is better
        s = s + "\t" + str(grdcselect[grdc_no][1][4])
    if grdc_no in validate:
        s = s + "\t" + validate[grdc_no][0] + "\t" + str( validate[grdc_no][5])

    print (s)
    s1 = s + "\n"
    f1 = open(grdc_Station, "a")
    f1.write(s1)
    f1.close()

print ("Done",nocal)