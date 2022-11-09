# Contrasting effects of shooting disturbance on the movement and behaviour of sympatric wildfowl species
Movement and behavioral effects of shooting on Barnacle Goose and White-fronted Geese on Islay, Scotland. 

## _Authors_
- Luke Ozsanlav-Harris
- Aimee Mcintosh

## Code description
- `R code/1_Combine tracking data + label with landcover.R`: Prepare GWfG tracking data for analysis
- `R code/2_Shooting proximity effects on Distance travelled.R`: Models to determine the immediate movement impact of shooting
- `R code/4_Daily behaviour on shooting vs non-shooting days.R`: Models to determine the behavioral changes in response to shooting in GWfG
- `R code/5_Daily ODBA on shooting vs non-shooting days.R`; Models to determine the change in daily ODBA due to shooting
- `R code/Useful Functions` folder containing useful code for this project, it is not used in the main work flow

## Data description
- `Landcover Data/Islay landcover data`: habitat data from Islay from the years 2015, 2017, 2018, 2019 amd 2020 as rasters. This was cropped from larger UK raster using `1_Crop RAW landcover data.R`
- `Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY` shapefile of all the agricultural field boundaries on Islay
- `MetaData/Tagged bird summary data new.csv` metadata file that contains the sex, ringing location and deployment dates of all tagged GWfG
- `Shooting logs/All_logs_cleaned.csv`: All cleaned shooting logs from Islay containing the information on each shooting event
- `Derived data/All_winter_GPS_with_habitat.RDS` R data file containing all of the winter Islay GPS fixes for GWfG with the habitats appended from `Landcover Data/Islay landcover data`

