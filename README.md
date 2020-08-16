# GMD_hydrostreamer
Code and (part of) the data needed to replicate the figures in the hydrostreamer 1.0 model description paper.
The code reproduces figures 3-5, table 4 and appendix tables A1 to A4.

The repository contains all data needed except for monitoring data, for which our license does not cover distribution.
The data is available from https://portal.mrcmekong.org. *Note that the code does not run without the observations.*

Other data included are:
* runoff datasets obtained from the ISIMIP archive (https://esg.pik-potsdam.de/search/isimip/), and cropped to the area-of-interest
* discharge data obtained from the ISIMIP archive, and cropped to the area-of-interest
* GRADES data for the 10 monitoring stations included, aggregated to monthly means. https://doi.org/10.1029/2019WR025287
* GLOFAS data for the 10 monitoring stations included, aggregated to monthly means. https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical
* HydroSHEDS river network for the area-of-interest. https://hydrosheds.org
* River segment specific catchment areas delineated using a Voronoi Diagram, and using HydroSHEDS drainage direction raster.
* Station locations for the 0.5 degree data, including the HydroSHEDS river segmeng ID for each station.

The ISIMIP data is distributed with CC-BY 4.0 or CC-NY-NC 4.0 license. Refer to the following link for which license applies to which model. https://www.isimip.org/gettingstarted/terms-of-use/licenses-publicly-available-isimip-data/?query=license#licences-for-isimip2a
GRADES data can be only used for research purposes.
