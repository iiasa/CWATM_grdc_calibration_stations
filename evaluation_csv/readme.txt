The use of GRDC gauging stations for calibrating large-scale hydrological models	
	
Peter Burek,1 Mikhail Smilovic1	
1 International Institute for Applied Systems Analysis, Laxenburg, Austria	
	
Abstract. The Global Runoff Data Centre provides time series of observed discharges that are very valuable for calibrating and validating the results of hydrological models. We address a common issue in large-scale hydrology which, though  investigated several times, has not been satisfactorily solved. Grid-based hydrological models need to fit the reported station location to the river network depending on the resolution, to compare simulated discharge with observed discharge. We introduce an Intersection over Union ratio approach to selected station locations on a coarser grid scale, reducing the errors in assigning stations to the wrong basin. We update the 10-year-old database of watershed boundaries with additional stations based on a high-resolution (3 arcseconds) river network, and we provide source codes and high- and low-resolution watershed boundaries.	
	
The first step toward generating a high-resolution collection of watershed shapefiles was to update the work of Lehner (2012) to include more basins (10,241 stations vs. 7,163), based on a higher resolution river network database (3’’ MERIT Hydro from Yamazaki et al. (2019) vs. 15’’ the HydroSHEDS from Lehner et al. (2008), including the changed GRDC IDs from September 2021.	
The second step, of generating a low-resolution collection of watershed shapefiles based on Intersection over Union ratio, was inspired by the ideas of Rezatofighi et al. (2019) and Munier and Decharme (2021). It is a better approach than selecting a station location on low-resolution river network systems based only on the upstream area and distance to the original location. Here, we provide the low-resolution watershed boundaries on 30’and 5’ and the source code to produce results for different resolutions and projection systems.	
The third step, of selecting suitable stations for calibration and validation, was also based on the Intersection over Union ratio. This selection of stations can now be used in a more effective way to calibrate grid-based hydrological models at different resolutions. 	
	
	
Content:	
1	Summary from originally 10701 stations to 10241 selected station with corrected location and basin shapefile
2	all 10701 station with 6 categories (1 smaller than 10km2, 2 no basin found, 3 area error > 50%,4 area >= 9000km2,5 area between 1000-9000km2,6 smaller than 1000km2)
3	station quality (area error - comparisaon between upstream area from MERIT UPA and provided upstream area, comparison between provided location and new location based on MERIT UPA)
4	GRDC station selected for hydrological modelling on 30 arcmin
5	GRDC station selected for hydrological modelling on 5 arcmin


-----------------------------
SUMMARY


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

-----------------------------------
Descriptions of csv
===================

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