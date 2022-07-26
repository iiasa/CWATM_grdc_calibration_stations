# The use of GRDC gauging stations for calibrating large-scale hydrological models


## Authors

Peter Burek, Mikhail Smilovic
International Institute for Applied Systems Analysis, Laxenburg, Austria

We address a small but annoying problem every grid-based hydrological model must solve to compare simulated with observed discharge. First, observed station location does not always fit with high-res river network. We updated the database with stations based on a new high-res network. Second, the corrected location does not fit a coarser grid-based river network. We use a new concept of shape-file similarity to fit station locations on a coarser grid, reducing the error in assigning stations to the wrong basin.


## Abstract

The Global Runoff Data Centre provides time series of observed discharges that are very valuable for calibrating and validating the results of hydrological models. We address a common issue in large-scale hydrology which, though  investigated several times, has not been satisfactorily solved. Grid-based hydrological models need to fit the reported station location to the river network depending on the resolution, to compare simulated discharge with observed discharge. We introduce an Intersection over Union ratio approach to selected station locations on a coarser grid scale, reducing the errors in assigning stations to the wrong basin. We update the 10-year-old database of watershed boundaries with additional stations based on a high-resolution (3 arcseconds) river network, and we provide source codes and high- and low-resolution watershed boundaries.

## Summary:

This paper describes the procedure used to generate a dataset of stations of observed discharge to be used at different resolutions for the calibration of large-scale hydrological models. It is based on the metadata of GRDC stations and on MERIT Hydro. The Python source code and dataset produced are freely available for download through a GitHub and a Zenodo repository.

The first step toward generating a high-resolution collection of watershed shapefiles was to update the work of Lehner (2012) to include more basins (10,241 stations vs. 7,163), based on a higher resolution river network database (3’’ MERIT Hydro from Yamazaki et al. (2019) vs. 15’’ the HydroSHEDS from Lehner et al. (2008), including the changed GRDC IDs from September 2021.

The second step, of generating a low-resolution collection of watershed shapefiles based on Intersection over Union ratio, was inspired by the ideas of Rezatofighi et al. (2019) and Munier and Decharme (2021). It is a better approach than selecting a station location on low-resolution river network systems based only on the upstream area and distance to the original location. Here, we provide the low-resolution watershed boundaries on 30’and 5’ and the source code to produce results for different resolutions and projection systems. 

The third step, of selecting suitable stations for calibration and validation, was also based on the Intersection over Union ratio. This selection of stations can now be used in a more effective way to calibrate grid-based hydrological models at different resolutions.


## Dataset produced

### Excel file: evaluation.xlsx

10241 station with corrected lat/lon and high-Res shapefile based on MERIT (Yamazaki et al. 2019)

953  stations for 30 arcmin calibration with corrected lat/lon and adjusted lat/lon for 30 arcmin

3917 stations for 5 arcmin calibration with corrected lat/lon and adjusted lat/lon for 5 arcmin


### CSV files of evaluation.xlsx in evaluation_csv

Description of CSV files in folder evaluation_csv and below


### Basin watershed boundaries shapefiles

10241 shapefiles for stations with corrected lat/lon and high-Res shapefile based on MERIT (Yamazaki et al. 2019)

2725 shapefiles for stations for 30 arcmin calibration with corrected lat/lon on 30 arcmin

6416 shapefiles for stations for 5 arcmin calibration with corrected lat/lon on 5 arcmin



## Data used

### GRDC data

Global Runoff Data Centre (GRDC) database
 
“A global hydrological database is essential for research and application-oriented hydrological and climatological studies at global, regional, and basin scales. The Global Runoff Database is a unique collection of river discharge data on a global scale. It contains time series of daily and monthly river discharge data of well over 10,000 stations worldwide. This adds up to around 470,000 station-years with an average record length of 45 years” (GRDC, 2020).

GRDC: https://www.bafg.de/GRDC/EN/Home/homepage_node.html, last access: 15/05.

Shapefile of basins from Lehner (2012)

Used for comparison. All 7163 shapefiles stored in one folder
https://www.bafg.de/GRDC/EN/02_srvcs/22_gslrs/222_WSB/methodology_Lehner.html
https://www.bafg.de/GRDC/EN/02_srvcs/22_gslrs/222_WSB/watershedBoundaries_node.html

Lehner, B.: Derivation of watershed boundaries for GRDC gauging stations based on the HydroSHEDS drainage network - Technical Report prepared for the GRDC, Global Runoff Data Centre, Koblenz, Germany, 18, 2012.

MERIT database (Yamazaki et al. 2019)

http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/

MERIT Hydro is a global hydrography datasets, developed based on the MERIT DEM and multiple inland water maps. It contains flow direction, flow accumulation, hydrologically adjusted elevations, and river channel width.

High-resolution raster hydrography maps are a fundamental data source for many geoscience applications. Here we introduce MERIT Hydro, a new global flow direction map at 3 arc-second resolution (~90 m at the equator) derived from the latest elevation data (MERIT DEM) and water body datasets (G1WBM, GSWO, and OpenStreetMap). We developed a new algorithm to extract river networks near-automatically by separating actual inland basins from dummy depressions caused by the errors in input elevation data. After a minimum amount of hand-editing, the constructed hydrography map shows good agreement with existing quality-controlled river network datasets in terms of flow accumulation area and river basin shape. The location of river streamlines was realistically aligned with existing satellite-based global river channel data. Relative error in the drainage area was smaller than 0.05 for 90% of GRDC gauges, confirming the accuracy of the delineated global river networks. Discrepancies in flow accumulation area were found mostly in arid river basins containing depressions that are occasionally connected at high water levels and thus resulting in uncertain watershed boundaries. MERIT Hydro improves on existing global hydrography datasets in terms of spatial coverage (between N90 and S60) and representation of small streams, mainly due to increased availability of high-quality baseline geospatial datasets. The new flow direction and flow accumulation maps, along with accompanying supplementary layers on hydrologically adjusted elevation and channel width, will advance geoscience studies related to river hydrology at both global and local scales.

- river network LDD in 5degx5deg chunks stored in 30degx30deg folders e.g. dir_n30e030
- upa_  (upstream area in km2 in 5x5deg) stored in 30degx30deg folders e.g. upa_n30e030

No further preprocessing is needed, just downloading and storing in the same structure as on the webpage.

Yamazaki D., D. Ikeshima, J. Sosa, P.D. Bates, G.H. Allen, T.M. Pavelsky, MERIT Hydro: A high-resolution global hydrography map based on latest topography datasets
Water Resources Research, vol.55, pp.5053-5073, 2019, doi: 10.1029/2019WR024873

## Python programmes

- 1_findMeritcoord.py
- 2_makeshape.py
- 3_makeshape.py
- 4_basincomp.py
- 5_stations4calib.py

- Edit the pathes in each script
- Run each script after the other
- Add manual corrected station after step 1


### Python libraries

Python libraries like geopandas, numpy, rasterio, osgeo

Python library flwdir from Eilander et al. (2021) 
https://pypi.org/project/pyflwdir/0.4.3/
used to create upstream maps, watershed maps from MERIT 3 arcsec river network maps

Eilander, D., van Verseveld, W., Yamazaki, D., Weerts, A., Winsemius, H. C., and Ward, P. J.: A hydrography upscaling method for scale-invariant parametrization of distributed hydrological models, Hydrol. Earth Syst. Sci., 25, 5287-5313, 10.5194/hess-25-5287-2021, 2021.

## Output

Excel file: evaluation.xlsx

-	summary from originally 10701 stations to 10241 selected station with corrected location and basin shapefile
-	all 10701 station with 6 categories (1 smaller than 10km2, 2 no basin found, 3 area error > 50%,4 area >= 9000km2,5 area between 1000-9000km2,6 smaller than 1000km2)
-	station quality (area error - comparisaon between upstream area from MERIT UPA and provided upstream area, comparison between provided location and new location based on MERIT UPA)
-	GRDC station selected for hydrological modelling on 30 arcmin
-	GRDC station selected for hydrological modelling on 5 arcmin


### CSV files of evaluation.xlsx in evaluation_csv

-> Final product: 

10241 station with correected lat/lon and high-Res shapefile based on MERIT (Yamazaki et al. 2019)

953  stations for 30 arcmin calibration with corrected lat/lon

3917 stations for 5 arcmin calibration with corrected lat/lon

Smoothen High-Res shapefiles of watershed boundaries on 3 arcsec

based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec
watershed boundaries of 10241 basins

Low-Res shapefiles of watershed boundaries on 30 arcmin

Selected 30 arcmin watershed basins based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec. Watershed boundaries of 2725 basins

Low-Res shapefiles of watershed boundaries on 5 arcmin

Selected 5 arcmin watershed basins based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec. Watershed boundaries of 6416 basins

## Python scripts

### 1_findMeritcoord.py

Use GRDC station cooordinates to find the outlet in MERIT upstream area 

Using the method of Lehner 2012 to find corrected station location

Uses data:

- grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
- grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID

Output:

- grdc_Merit_1.txt

### 2_makeshape.py

Use MERIT coordinates of upstream area to create shapefiles

uses upstream area of MERIT (UPA) and GRDC station data

to create shapefiles  from rivernetwork MERIT data

-> adding manual corrected station location before

Uses data:

-	makeshape_all10349.txt  station with new location fitted to merit UPA  and manual corrected station location
-	grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
-	grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
-	Shapefiles of basins from Lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
-	Shapefiles of basins from Lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_

Output:

-  grdc_Merit_2.txt
-  Shapefiles of basins : ashape_merit/grdc_basin_merit_


### 3_makeshape.py

Use MERIT coordinates of upstream area create smoothen shapefiles

uses upstream area of MERIT (UPA) and GRDC station data

to create smooothen shapefiles  from rivernetwork MERIT data

Uses data

-	grdc_MERIT_1.txt  station with new location fitted to merit UPA (output from  1_findMeritcoord.py)
-	grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
-	grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
-	Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
-	Shapefiles of basins : ashape_merit/grdc_basin_merit_

output

-  grdc_shape_all.txt
-  Shapefiles of basins : ashape2_merit/grdc_basin_merit_smooth_

### 4_basincomp.py

Creates shapefiles and a list in lower resolution 5 arcmin and  30arcmin

used:

- Shapefiles of basins : ashape2_merit/grdc_basin_merit_smooth_
- grdc_shape_all.txt with 10249 stations
- River network as pcraster format LDD: ldd.tif for 5 arcmin or 30arcmin

output either as 5 arcmin or 30 arcmin:

- shapefiles: /ashape30min/grdc_basin_30min_
- shapefiles: /ashape5min/grdc_basin_5min_
- results/basins_30min.txt"
- results/basins_5min.txt"

### 5_stations4calib.py

Create list of stations to calibrate or validate

using the watershed boundary shapepfile on 30min (or 5min) and the station metadata

to dismiss stations which are not qualify and to select vetween calibration and validation stations

Validation station are to close to calibration stations, but have a lower score

Either for 5arcmin or 30arcmin

input:

- shapefile of 30min basins: ashape30min/grdc_basin_30min_basin_
- station information in: basins_30min4.txt

Output:

Station information in: results/grdc_station_sel30min.txt"

-	if a station can be used for calibration it gets a C as indicator
-	if a station passes the selection but is to close to another station which is used for calibration it gets a V
-	if a station does not pass the selection criteria it gets a N


Selection criteria:

-	similarityMin: Minimum similarity between 3 arcsec basin and low-res basin: here: 0.7 (1 = identically, 0 = nothing in common)
-	areaaccorMin: Minimim area error: area provided from GRDC and area from low-res shapefile should have at least (e.g 0.4) 
-	area = 0.3 x areaGRDC is dismissed , area = 0.5 x areaGRDC is used
-	t_yrsMin = 5:  5 minimum years of record
-	t_endMin = 1985: last year should be 1985 or younger
-	d_endMin = 1985


Scoring points:

if 2 (or more) station are too close to each other the one is taken with the most scoring points

-	simiScore = 2  # every 2% scoring +1
-	areaScore = 2   # every 2% scoring +1
-	tyrsScore = 5  # every 5 years a score
-	tyrsScoreMax = 10 # max score for end years
-	tendScore = 3  # every 3 years a point
-	tendScoremax = 100 # maximum score from this can only be 100
-	dayscore = 5  # points if daily timeseries
-	MissingScore = 5  # neg. points for missing values

## Results from python scripts

### grdc_Merit_1.txt

output from 1_findMeritcoord.py

-	No: Number from 1 ...
-	GRDC_No: GRDC number
-	lat: original latitude from GRDC metafile
-	lon: original longituted from GRDC metafile
-	newlat: corrected latitude based on MERIT UPA dataset
-	newlon: corrected longitute based on MERIT UPA dataset
-	area; provided basin area from GRDC metafile
-	newarea: basin area based on MERIT UPA dataset
-	UPS_Indicator:  min error in % from MERIT UPA to provided basin area
-	dist_Indicator: distance to original pour point in [unit:100m]
-	Indicator:  ranking criteria: UPS_Indicator + 2 x dist_indicator

### grdc_Merit_2.txt

output from 2_makeshape1.py

-	No: Number from 1 ...
-	GRDC_No: GRDC number
-	lat: original latitude from GRDC metafile
-	lon: original longituted from GRDC metafile
-	area; provided basin area from GRDC metafile
-	newlat: corrected latitude based on MERIT UPA dataset
-	newlon: corrected longitute based on MERIT UPA dataset
-	newarea: basin area based on MERIT UPA dataset
-	shapearea: area of the shape calculated with geopandas directly from shape
-	Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
-	lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
-	Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar


### grdc_shape_allend_1.txt

output from 3_makeshape.py

-	GRDC_No: GRDC number
-	lat: original latitude from GRDC metafile
-	lon: original longituted from GRDC metafile
-	newlat: corrected latitude based on MERIT UPA dataset
-	newlon: corrected longitute based on MERIT UPA dataset
-	GRDCarea; provided basin area from GRDC metafile
-	area: basin area based on MERIT UPA dataset
-	shapearea: area of the shape calculated with geopandas directly from shape
-	Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
-	lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
-	Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar


### basins_30min.txt

output from 4_basincomp.py

-	No: Number from 1
-	GRDC_No: GRDC number
-	similarity: similarity of 30min low-res shape with high-res shape from 3arcsec
-	areaGRDC: area provided by GRDC
-	area: area from high-res UPA MERIT at pour point
-	lat: coorected latitude on high-res
-	lon: corrected longitude on high-res
-	area30min : area on pour point of best fitting 30arcmin grid cell
-	lat30min: latitude on 30arcmin
-	lon30min: longitude on 30arcmin
-	lat30move: latitude moved to grid cell centre
-	lon30move: longitude moved to grid cell centre
-	indloc: locatation of indice which is taken as nbest fit: a 5x5 gridcell surrounding is taken numbered 0-24
-	upsloc: location (index of 0-24) of best fitting upsstream area
-	ind:    scoring result for best upstream fitting based on upstream area fitting and shape similarity
-	indups: scoring on upstream area fitting between low-res area and high-res area
-	ups1:   low-res upstream area of best upstream fitting
-	indshape: scoring on similarity between low-res shape and high res shape
-	shapeloc: location (index of 0-24) of best fitting shape similarity
-	ind:    scoring result for best upstream fitting based on upstream area fitting and shape similarity
-	indups: scoring on upstream area fitting between low-res area and high-res area
-	ups1:   low-res upstream area of best shape fitting
-	indshape: scoring on similarity between low-res shape and high res shape


### grdc_station_sel30min.txt

output for 5_stations4calib.py

-	No: Number from 0 ...
-	GRDC_No: GRDC number
-	latGRDC: original latitude from GRDC metafile
-	lonGRDC: original longituted from GRDC metafile
-	latcorrected: coorected latitude on high-res
-	loncorrected: corrected longitude on high-res
-	lat30min: latitude moved to grid cell centre
-	lon30min: longitude moved to grid cell centre
-	areaGRDC: area provided by GRDC
-	area: area from high-res UPA MERIT at pour point
-	area30min : area on pour point of best fitting 30arcmin grid cell
-	similarity: similarity of 30min low-res shape with high-res shape from 3arcsec
-	area_accordance> comparison of area30min with areaGRDC
-	calibrate_orNot> can be use for calibration (C), validation (V) or is dismissed (N)
-	scoring: scoring based on the scoring criteria
-	grdc_scored_better: grdc number with scored better
-	better_score: better score of the grdc no which is used for calibration



## Description of CSV files in folder evaluation_csv

	
### Content

1	Summary from originally 10701 stations to 10241 selected station with corrected location and basin shapefile
2	all 10701 station with 6 categories (1 smaller than 10km2, 2 no basin found, 3 area error > 50%,4 area >= 9000km2,5 area between 1000-9000km2,6 smaller than 1000km2)
3	station quality (area error - comparisaon between upstream area from MERIT UPA and provided upstream area, comparison between provided location and new location based on MERIT UPA)
4	GRDC station selected for hydrological modelling on 30 arcmin
5	GRDC station selected for hydrological modelling on 5 arcmin


### SUMMARY

No	Overview	Number of stations
1	All station	10701
2	station <10km	-124
3	station area =>10km or -999	10577
4	station not found	-228
5	stations found location	10349
6	Manuall corrected	169
7	station area error>60%	-108
8	stations area error =<60%	10241
		
	station in this selection with no area 49	-49
	station with area	10192
		
	Stations altogeher with no area (-999)	327
	GRDC ID changed	1359




MERIT area vs. GRDC area		
		
Station area error	Sum number of stations	% stations
< 1	4332	41
< 5	7920	75
< 10	8980	85
< 15	9446	90
< 25	9862	94
< 50	10174	97
< 40	10300	98
not found	10528	100
no provided area	10577	
		
Stations not found		228
Stations not found or an error >50%		354
Stations not found or an error >25%		665
		437
Distance of new location from provided location (10241 stations)		
Distance from provided loc [km]	Sum number of station	
<= 1		8544
> 1		1697
> 5		528
> 10		240
> 50		75


### Descriptions of csv

Summary.csv
-----------

column 1: row number
column 2: description 
column 3: Variable

rows	Description
1-3		(1) 10701 station GRDC
4-13	(2) 124 station < 10km2
14-23	(3) 10577 stations area >= 10km2 or -999
24-34	(4) 228 Station with no location 
35-47	(5) 10349 with location on MERIT in(ashape_merit)  (49 with no area in reported area)
48-60	(6) 169 stations manuall checked
61-74	(7) 108 Station more than 60% area error
75-86	(8) 10241 Station less equal than 60% area error (49 without provided area)  
87-96	327 stations with no area
97-99	Former ID changed

allstations.csv
---------------

column 1: row number
column 2: description 
column 3: Variable

Categorie	
1	smaller than 10 km2
2	no basin found
3	Area error >0.5
4	>= 9000 km2
5	between 1000-9000km2
6	smaller than 1000 km2

rows	Description
1-6		(1) 10701 station GRDC with the categories described above

stationquality.csv
------------------

column 1: row number
column 2: description 
column 3: Variable

rows	Description
1-10	10241 Station less eqal than 60% area error (49 without provided area)
11-20	10241 station sorted by area error
21-30	10241 station sorted by distance for provided location
31-41	228 Station with no location
42-49	108 Station less than 60% area eror

GRDC_station_select_for_30min.csv
---------------------------------

Information
Calibration (C,V,N)
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


column 1: row number
column 2: description 
column 3: Variable

rows	Description
1-17	GRDC station select for 30min
18-32	For Calibration


GRDC_station_select_for_5min.csv
--------------------------------

column 1: row number
column 2: description 
column 3: Variable

rows	Description
1-17	GRDC station select for 5min
18-32	For Calibration 5arcmin
