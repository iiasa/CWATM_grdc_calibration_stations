# The use of GRDC gauging stations for calibrating large-scale hydrological models


## Abstract

The Global Runoff Data Centre provides time series of observed discharges that are very valuable for calibrating and validating the results of hydrological models. We address a common issue in large-scale hydrology which, though  investigated several times, has not been satisfactorily solved. Grid-based hydrological models need to fit the reported station location to the river network depending on the resolution, to compare simulated discharge with observed discharge. We introduce an Intersection over Union ratio approach to selected station locations on a coarser grid scale, reducing the errors in assigning stations to the wrong basin. We update the 10-year-old database of watershed boundaries with additional stations based on a high-resolution (3 arcseconds) river network, and we provide source codes and high- and low-resolution watershed boundaries.

## Authors
Peter Burek, Mikhail Smilovic
International Institute for Applied Systems Analysis, Laxenburg, Austria


##Data used:

### GRDC data

Global Runoff Data Centre (GRDC) database 
“A global hydrological database is essential for research and application-oriented hydrological and climatological studies at global, regional, and basin scales. The Global Runoff Database is a unique collection of river discharge data on a global scale. It contains time series of daily and monthly river discharge data of well over 10,000 stations worldwide. This adds up to around 470,000 station-years with an average record length of 45 years” (GRDC, 2020).
GRDC: https://www.bafg.de/GRDC/EN/Home/homepage_node.html, last access: 15/05.

### Shapefile of basins from Lehner, B (2012)

Shapefiles from Lehner (2012) 
where used for comparison
https://www.bafg.de/GRDC/EN/02_srvcs/22_gslrs/222_WSB/methodology_Lehner.html
https://www.bafg.de/GRDC/EN/02_srvcs/22_gslrs/222_WSB/watershedBoundaries_node.html

Lehner, B.: Derivation of watershed boundaries for GRDC gauging stations based on the HydroSHEDS drainage network - Technical Report prepared for the GRDC, Global Runoff Data Centre, Koblenz, Germany, 18, 2012.

## GSIM
The Global Streamflow Indices and Metadata Archive (GSIM) for comparison

Do, H.; Gudmundsson, L.; Leonard, M.; Westra, S.; Senevirante, S.I (2017): The Global Streamflow Indices and Metadata Archive (GSIM) - Part 1: The production of daily streamflow archive and metadata

### MERIT database (Yamazaki et al. 2019)

http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/
MERIT Hydro is a global hydrography datasets, developed based on the MERIT DEM and multiple inland water maps. It contains flow direction, flow accumulation, hydrologically adjusted elevations, and river channel width.
High-resolution raster hydrography maps are a fundamental data source for many geoscience applications. Here we introduce MERIT Hydro, a new global flow direction map at 3 arc-second resolution (~90 m at the equator) derived from the latest elevation data (MERIT DEM) and water body datasets (G1WBM, GSWO, and OpenStreetMap). We developed a new algorithm to extract river networks near-automatically by separating actual inland basins from dummy depressions caused by the errors in input elevation data. After a minimum amount of hand-editing, the constructed hydrography map shows good agreement with existing quality-controlled river network datasets in terms of flow accumulation area and river basin shape. The location of river streamlines was realistically aligned with existing satellite-based global river channel data. Relative error in the drainage area was smaller than 0.05 for 90% of GRDC gauges, confirming the accuracy of the delineated global river networks. Discrepancies in flow accumulation area were found mostly in arid river basins containing depressions that are occasionally connected at high water levels and thus resulting in uncertain watershed boundaries. MERIT Hydro improves on existing global hydrography datasets in terms of spatial coverage (between N90 and S60) and representation of small streams, mainly due to increased availability of high-quality baseline geospatial datasets. The new flow direction and flow accumulation maps, along with accompanying supplementary layers on hydrologically adjusted elevation and channel width, will advance geoscience studies related to river hydrology at both global and local scales.

- river network LDD in 5degx5deg chunks
- upa_  (upstream area in km2 in 5x5deg)

Yamazaki D., D. Ikeshima, J. Sosa, P.D. Bates, G.H. Allen, T.M. Pavelsky, MERIT Hydro: A high-resolution global hydrography map based on latest topography datasets
Water Resources Research, vol.55, pp.5053-5073, 2019, doi: 10.1029/2019WR024873


### Python libraries

Usual libraries like geopandas, numpy, rasterio, osgeo
and

Python library flwdir from Eilander et al. (2021) 
Eilander, D., van Verseveld, W., Yamazaki, D., Weerts, A., Winsemius, H. C., and Ward, P. J.: A hydrography upscaling method for scale-invariant parametrization of distributed hydrological models, Hydrol. Earth Syst. Sci., 25, 5287-5313, 10.5194/hess-25-5287-2021, 2021.

## Output:

### Excel file

- summary from originally 10701 stations to 10241 selected station with corrected location and basin shapefile
- all 10701 station with 6 categories (1 smaller than 10km2, 2 no basin found, 3 area error > 50%,4 area >= 9000km2,5 area between 1000-9000km2,6 smaller than 1000km2)
- station quality (area error - comparisaon between upstream area from MERIT UPA and provided upstream area, comparison between provided location and new location based on MERIT UPA)
- GRDC station selected for hydrological modelling on 30 arcmin
- GRDC station selected for hydrological modelling on 5 arcmin

### Smoothen Hig-res shapefiles of watershed boundaries on 3 arcsec

based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec
watershed boundaries of 10241 basins

### Low-res shapefiles of watershed boundaries on 30 arcmin

Selected 30 arcmin watershed basins based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec 
watershed boundaries of 2725 basins

### Low-res shapefiles of watershed boundaries on 5 arcmin

Selected 5 arcmin watershed basins based on MERIT (Yamazaki et al. 2019) river network on 3 arcsec 
watershed boundaries of 6416 basins

## Python scrips

### 1_findMeritcoord.py

Use GRDC station cooordinates to find the outlet in MERIT upstream area 
Using the method of Lehner 2012 to find corrected station location

Uses data:

- grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
- grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID


Output:

- grdc_Merit_1.txt

### 2_makeshape1.py

Use MERIT coordinates of upstream area to create shapefiles
uses upstream area of MERIT (UPA) and GRDC station data
to create shapefiles  from rivernetwork MERIT data

Uses data:

- makeshape_all10349.txt  station with new location fitted to merit UPA  and manual corrected station location
- grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
- grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
- Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_

output:

-  grdc_Merit_2.txt
-  Shapefiles of basins : ashape_merit/grdc_basin_merit_


### 3_makeshape

Use MERIT coordinates of upstream area create smoothen shapefiles
uses upstream area of MERIT (UPA) and GRDC station data
to create smooothen shapefiles  from rivernetwork MERIT data
-> adding manual coorected station location before

Uses data:

- grdc_MERIT_1.txt  station with new location fitted to merit UPA (output from  1_findMeritcoord.py)
- grdc_2022_10701.txt: all data from GRDC metafile 2022 with 10701 stations
- grdc2022_lehner2012.txt  - translation file new GRDC ID - old GRDC ID
- Shapefiles of basins from lehner 2012: grdc_station_upstream/grdc_basins_smoothed_md_no_
- Shapefiles of basins : ashape_merit/grdc_basin_merit_

output:

-  grdc_shape_all.txt
-  Shapefiles of basins : ashape2_merit/grdc_basin_merit_smooth_

### 4

Creates shapefiles and a list in lower resolution 5 arcmin and  30arcmin











## Results from python scripts

### grdc_Merit_1.txt

output from 1_findMeritcoord.py

- No: Number from 1 ...
- GRDC_No: GRDC number
- lat: original latitude from GRDC metafile
- lon: original longituted from GRDC metafile
- newlat: corrected latitude based on MERIT UPA dataset
- newlon: corrected longitute based on MERIT UPA dataset
- area; provided basin area from GRDC metafile
- newarea: basin area based on MERIT UPA dataset
- UPS_Indicator:  min error in % from MERIT UPA to provided basin area
- dist_Indicator: distance to original pour point in [unit:100m]
- Indicator:  ranking criteria: UPS_Indicator + 2 x dist_indicator

### grdc_Merit_2

output from 2_makeshape1.py

- No: Number from 1 ...
- GRDC_No: GRDC number
- lat: original latitude from GRDC metafile
- lon: original longituted from GRDC metafile
- area; provided basin area from GRDC metafile
- newlat: corrected latitude based on MERIT UPA dataset
- newlon: corrected longitute based on MERIT UPA dataset
- newarea: basin area based on MERIT UPA dataset
- shapearea: area of the shape calculated with geopandas directly from shape
- Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
- lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
- Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar