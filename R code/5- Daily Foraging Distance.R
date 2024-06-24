## Aimée McIntosh

## 19/05/2023

## Determining whether GPS tracking data occurs on shooting disturbed days 
# Examining the effect that shooting has on movement, in particular the total daily distance travelled on days disturbed by shooting vs those undisturbed by shooting. 
# Idea here is to model disturbed vs undisturbed daily travel distances
# This code uses data files that can be created using the "Identifying_Shooting_Disturbance" code, files are also provided in the data folder 
# Code is initially split for GBG and GWfG to assign sampling rate and disturbance days before being combined for analysis 

## ## Packages required
pacman::p_load(tidyr,plyr,dplyr,zoo,data.table,move,ggplot2,purrr,readr,
               lubridate,tidyverse,readr,amt,geosphere,sf,DHARMa,lme4,
               MuMIn,effects,lattice,performance,MASS,glmmTMB,bbmle, emmeans)

# ------------------------------------ #
####     Prepare data for models    ####
# ------------------------------------ #

#### 1. GBG data ####
## 1.1 GPS data ##----
# GBG - data prep and distance metric calculation #

gps_gbg_all <- read.csv("tracking data/GBG_GPS_CLEAN_COMPLETE.csv")
head(gps_gbg_all)

# GBG list of shooting disturbance days, fields and IDs
gbg_shoot_id <- read.csv("Derived data/GBG_INDIV_SHOOT_DATE_797m.csv")
head(gbg_shoot_id)

# Prepare the data
gbg_shoot_id <- dplyr::select(gbg_shoot_id, c("Date_ID", "spatial_overlap", "Raw_Farm_Name", "fix_gps", "Shoot_field"))

# Create date ID column for gps data
# NB! Columns match the shooting data records
gps_gbg_all <- gps_gbg_all %>% dplyr::mutate(Date_ID = paste(Device_ID, UTC_Date, sep = "_"))
head(gps_gbg_all)

## 1.2 Combine shooting data with gps data ## ----
# - join by Date_ID
gbg_shoot <- dplyr::left_join(gps_gbg_all, gbg_shoot_id, by = "Date_ID")
head(gbg_shoot)

# NA in "spatial overlap" are non-shooting, change to 0 
gbg_shoot$spatial_overlap[is.na(gbg_shoot$spatial_overlap)] <- 0

# Streamline data
gbg_shoot2 <- dplyr::select(gbg_shoot, -c("X"))
head(gbg_shoot2)

# Filter to daytime only
gbg_shoot_day <- dplyr::filter(gbg_shoot2, tod_ == "day")
gbg_shoot_day$Timestamp <- as.POSIXct(gbg_shoot_day$Timestamp)

# Add field code to GPS fixes
islay.field <- readOGR("Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY/AS_ISLAY_GMS_FIELD_BOUNDARY.shp")
islay.field <- spTransform(islay.field, (CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")))
crs(islay.field)

gbg_winter_trk <- make_track(gbg_shoot_day, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gbg_winter_trk)
gbg_trk <- transform_coords(gbg_winter_trk, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 
str(gbg_trk)
crs(gbg_trk)
coordinates(gbg_trk) <- coordXY # extract coordinates from day data
proj4string(gbg_trk) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gbg_gps_shoot_field <- raster::extract(islay.field, gbg_trk)
gbg_gps_field <- cbind(gbg_shoot_day, gbg_gps_shoot_field)
head(gbg_gps_field)

# Add data for most frequently visited farm from GPS data - this gives an indication of the overal area utilised by geese in a given day 
farm <- read.csv("Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY/FARM_ID.csv")
farm <- dplyr::rename(farm, FIELD_ID = "Field_code")

# Add farm fix data
gbg_farm <- dplyr::left_join(gbg_gps_field, farm, by = "FIELD_ID")
head(gbg_farm)

gbg_farm_total <- gbg_farm %>% group_by(Date_ID) %>% summarize (mode_hld_farm =names(which.max(table(hld_farm))))

# Bind data frame summary with original GBG points data frame
gbg_joined_df <- dplyr::full_join(gbg_farm, gbg_farm_total, by = c("Date_ID"))
head(gbg_joined_df)

gbg_shoot_day <- dplyr::select(gbg_joined_df, -c("X", "point.ID", "poly.ID", "OID_T", "SHAPE_AREA", "SHAPE_LEN"))
head(gbg_shoot_day)

## 1.3 Re-sample GPS fixes for consistent sampling effort ## ----
# Convert Timestamp to posixct
gbg_shoot_day$Timestamp <- as.POSIXct(gbg_shoot_day$Timestamp)

# Split data by winter
gbg_19 <- dplyr::filter(gbg_shoot_day, Winter == 2019)
gbg_20 <- dplyr::filter(gbg_shoot_day, Winter == 2020)

# Sub-sample for consistent sampling rate to calculate daily distance traveled without being influenced by variable sampling rate
# NB! Using sub-sampling code from Liam Langley

trackSubSamp <- function(TD, dt=1, unit='hours'){
  
  ## order the timestamps
  TD <- TD[order(TD$Timestamp),] 
  
  # breakdown to datasets per bird
  unid = unique(TD$Device_ID) # the unique tag IDs in your data set
  nrid = length(unid)
  TDall = list(nrid)  
  TDred = list(nrid)
  timestep = paste(dt,unit)
  
  # create time sequence from min to max time with step sizes that are defined at the start of the function
  dati_start = min(TD$Timestamp,na.rm=T)
  dati_end = max(TD$Timestamp,na.rm=T)
  datiseq = seq(from=dati_start, to=dati_end, by=timestep)
  
  for (i in 1:nrid)
  {
    Dtemp = TD[TD$Device_ID == unid[i],]
    idx = sapply(datiseq, function( x) which.min( abs( difftime( Dtemp$Timestamp, x, units='mins') ) ) ) # finds closest time in data to your created time series
    TDall[[i]] = Dtemp
    TDred[[i]] = unique( Dtemp[idx,] ) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
  }
  
  TDred # return this list
}

# Sub-sample the data for each winter separately
### 2019 ###
list_gbg_19 <- trackSubSamp(gbg_19)

### 2020 ###
list_gbg_20 <- trackSubSamp(gbg_20)

# Convert list to data frame
### 2019 ###
df_gbg_19_ss <- do.call("rbind", list_gbg_19)
head(df_gbg_19_ss)

### 2020 ###
df_gbg_20_ss <- do.call("rbind", list_gbg_20)
head(df_gbg_20_ss)

# Combine the data frames
gbg_ss <- rbind(df_gbg_19_ss, df_gbg_20_ss)
str(gbg_ss)

## 1.4 Identify goose days ## ----
# Need to filter to minimum number of fixes in a day

# check number of fixes per day
gbg_ss$sum <- 1

sun_dates <- seq.Date(as.Date(min(gbg_ss$Timestamp)), as.Date(max(gbg_ss$Timestamp)), by = 1)
sunt <-
  getSunlightTimes(
    date = sun_dates,
    keep = c("dusk", "dawn"),
    lat = 55.738,
    lon = 6.246,
    tz = "UTC"
  )

## Add a column for the date
gbg_ss$date <- as.Date(gbg_ss$Timestamp)

## Now join on the susnet and sunrise times
sunt$lat <- NULL
sunt$lon <- NULL
gbg_ss2 <- left_join(gbg_ss, sunt,  by = "date")
gbg_ss2$DayNight <- ifelse(gbg_ss2$Timestamp < gbg_ss2$dawn | gbg_ss2$Timestamp > gbg_ss2$dusk,
                           "night", "day")
head(gbg_ss2)
gbg_ss3 <- dplyr::filter(gbg_ss2, DayNight == "day")
gbg_ss3

gbg_ss_summary <- with(gbg_ss3, aggregate(sum, by = list(Date_ID), "sum"))

names(gbg_ss_summary)[1] <- "Date_ID"
names(gbg_ss_summary)[2] <- "Winter"
names(gbg_ss_summary)[2] <- "Total_fixes"

# Bind data frame summary with original gpg points data frame
gbg_joined_df <- dplyr::full_join(gbg_ss3, gbg_ss_summary, by = c("Date_ID"))
head(gbg_joined_df)

# Subset the original gps points based on the number of fixes
gbg_df_points_cleaned <- dplyr::filter(gbg_joined_df, gbg_joined_df$Total_fixes >= 5)
gbg_df_points_cleaned

gbg_df_points_cleaned$Species <- "GBG"
head(gbg_df_points_cleaned)
min(gbg_df_points_cleaned$Total_fixes)
max(gbg_df_points_cleaned$Total_fixes)



#### 2. GWfG data ####
## 2.1 GPS data ## ----
gps_gwf_all <- read.csv("tracking data/GWFG_GPS_CLEAN_COMPLETE_hab.csv")

# GBG list of shooting disturbance days, fields and IDs
gwf_shoot_id <- read.csv("Derived data/GWFG_INDIV_SHOOT_DATE_797m.csv")

gwf_shoot_id <- dplyr::select(gwf_shoot_id, c("Date_ID", "spatial_overlap", "Raw_Farm_Name", "fix_gps", "Shoot_field"))

# Create date ID column for gps data
# Note: Columns match the shooting data records
gps_gwf_all$Timestamp <- as.POSIXct(gps_gwf_all$Timestamp, format = "%Y-%m-%d %H:%M:%S")

gps_gwf_all <- mutate(gps_gwf_all, UTC_Date = format(gps_gwf_all$Timestamp, '%Y-%m-%d'))

gps_gwf_all <- gps_gwf_all %>% dplyr::mutate(Date_ID = paste(Device_ID, UTC_Date, sep = "_"))

## 2.2 Combine shooting data with gps data ## ----
# join by Date_ID
gwf_shoot <- dplyr::left_join(gps_gwf_all, gwf_shoot_id, by = "Date_ID")

# NA in "spatial overlap" are non-shooting, change to 0 
gwf_shoot$spatial_overlap[is.na(gwf_shoot$spatial_overlap)] <- 0

# Streamline data
gwf_shoot2 <- dplyr::select(gwf_shoot, -c( "X"))

# Filter to daytime only
gwf_shoot_day <- dplyr::filter(gwf_shoot2, tod_ == "day")
gwf_shoot_day$Timestamp <- as.POSIXct(gwf_shoot_day$Timestamp)

# Split into each winter as too many fixes causes to crash
gwf_1213 <- dplyr::filter(gwf_shoot_day, Timestamp < "2014-10-01")
gwf_1415 <- dplyr::filter(gwf_shoot_day, Timestamp >= "2014-10-01", Timestamp < "2016-10-01")
gwf_1617 <- dplyr::filter(gwf_shoot_day, Timestamp >= "2016-10-01", Timestamp < "2018-10-01")
gwf_18oct <- dplyr::filter(gwf_shoot_day, Timestamp >= "2018-10-01", Timestamp < "2019-01-01")
gwf_18jan <- dplyr::filter(gwf_shoot_day, Timestamp >= "2019-01-01", Timestamp < "2019-03-01")
gwf_18apr <- dplyr::filter(gwf_shoot_day, Timestamp >= "2019-03-01", Timestamp < "2019-05-01")
gwf_19oct <- dplyr::filter(gwf_shoot_day, Timestamp >= "2019-10-01", Timestamp < "2020-01-01")
gwf_19jan <- dplyr::filter(gwf_shoot_day, Timestamp >= "2020-01-01", Timestamp < "2020-03-01")
gwf_19apr <- dplyr::filter(gwf_shoot_day, Timestamp >= "2020-03-01", Timestamp < "2020-05-01")
gwf_2021 <- dplyr::filter(gwf_shoot_day, Timestamp >= "2020-10-01", Timestamp < "2021-10-01")
gwf_2122 <- dplyr::filter(gwf_shoot_day, Timestamp >= "2021-10-01", Timestamp < "2022-10-01")

# Add field code to GPS fixes
islay.field <- readOGR("Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY/AS_ISLAY_GMS_FIELD_BOUNDARY.shp")
islay.field <- spTransform(islay.field, (CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")))
crs(islay.field)

# 2012/2013
gwf_winter_trk_1213 <- make_track(gwf_1213, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_1213)

gwf_trk_1213 <- transform_coords(gwf_winter_trk_1213, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_1213) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_1213) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_1213 <- extract(islay.field, gwf_trk_1213)

gwf_gps_field_1213<- cbind(gwf_1213, gwf_gps_shoot_field_1213)

# 2014/2015
gwf_winter_trk_1415 <- make_track(gwf_1415, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_1415)

gwf_trk_1415 <- transform_coords(gwf_winter_trk_1415, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_1415) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_1415) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_1415 <- extract(islay.field, gwf_trk_1415)

gwf_gps_field_1415<- cbind(gwf_1415, gwf_gps_shoot_field_1415)

# 2016/2017
gwf_winter_trk_1617 <- make_track(gwf_1617, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_1617)

gwf_trk_1617 <- transform_coords(gwf_winter_trk_1617, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_1617) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_1617) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_1617 <- extract(islay.field, gwf_trk_1617)

gwf_gps_field_1617<- cbind(gwf_1617, gwf_gps_shoot_field_1617)

# 2018 October
gwf_winter_trk_18oct <- make_track(gwf_18oct, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_18oct)

gwf_trk_18oct <- transform_coords(gwf_winter_trk_18oct, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_18oct) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_18oct) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_18oct <- extract(islay.field, gwf_trk_18oct)

gwf_gps_field_18oct<- cbind(gwf_18oct, gwf_gps_shoot_field_18oct)

# 2018 January
gwf_winter_trk_18jan <- make_track(gwf_18jan, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_18jan)

gwf_trk_18jan <- transform_coords(gwf_winter_trk_18jan, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_18jan) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_18jan) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_18jan <- extract(islay.field, gwf_trk_18jan)

gwf_gps_field_18jan<- cbind(gwf_18jan, gwf_gps_shoot_field_18jan)

# 2018 April
gwf_winter_trk_18apr <- make_track(gwf_18apr, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_18apr)

gwf_trk_18apr <- transform_coords(gwf_winter_trk_18apr, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_18apr) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_18apr) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_18apr <- extract(islay.field, gwf_trk_18apr)

gwf_gps_field_18apr<- cbind(gwf_18apr, gwf_gps_shoot_field_18apr)

# 2019 October
gwf_winter_trk_19oct <- make_track(gwf_19oct, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_19oct)

gwf_trk_19oct <- transform_coords(gwf_winter_trk_19oct, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_19oct) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_19oct) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_19oct <- extract(islay.field, gwf_trk_19oct)

gwf_gps_field_19oct<- cbind(gwf_19oct, gwf_gps_shoot_field_19oct)

# 2019 January
gwf_winter_trk_19jan <- make_track(gwf_19jan, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_19jan)

gwf_trk_19jan <- transform_coords(gwf_winter_trk_19jan, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_19jan) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_19jan) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_19jan <- extract(islay.field, gwf_trk_19jan)

gwf_gps_field_19jan<- cbind(gwf_19jan, gwf_gps_shoot_field_19jan)

# 2019 April
gwf_winter_trk_19apr <- make_track(gwf_19apr, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_19apr)

gwf_trk_19apr <- transform_coords(gwf_winter_trk_19apr, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_19apr) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_19apr) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_19apr <- extract(islay.field, gwf_trk_19apr)

gwf_gps_field_19apr<- cbind(gwf_19apr, gwf_gps_shoot_field_19apr)

# 2020/2021
gwf_winter_trk_2021 <- make_track(gwf_2021, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_2021)

gwf_trk_2021 <- transform_coords(gwf_winter_trk_2021, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_2021) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_2021) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_2021 <- extract(islay.field, gwf_trk_2021)

gwf_gps_field_2021<- cbind(gwf_2021, gwf_gps_shoot_field_2021)

# 2021/2022
gwf_winter_trk_2122 <- make_track(gwf_2122, .x = Longitude, .y = Latitude, .t = Timestamp, id = Device_ID, crs = CRS("+proj=longlat +dataum=WGS84"), all_cols = TRUE)
str(gwf_winter_trk_2122)

gwf_trk_2122 <- transform_coords(gwf_winter_trk_2122, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

coordXY <- c("x_", "y_")  # make spatial points data frame 

coordinates(gwf_trk_2122) <- coordXY # extract coordinates from day data
proj4string(gwf_trk_2122) <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

gwf_gps_shoot_field_2122 <- extract(islay.field, gwf_trk_2122)

gwf_gps_field_2122<- cbind(gwf_2122, gwf_gps_shoot_field_2122)

# Combine all together to sing GWfG file
gwf_gps_field <- rbind(gwf_gps_field_1213,gwf_gps_field_1415,gwf_gps_field_1617,gwf_gps_field_18oct,gwf_gps_field_18jan,gwf_gps_field_18apr,gwf_gps_field_19oct,gwf_gps_field_19jan,gwf_gps_field_19apr, gwf_gps_field_2021,gwf_gps_field_2122)
head(gwf_gps_field)

# Add farm data 
farm <- read.csv("Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY/FARM_ID.csv")
head(farm)
gwf_gps_field <- dplyr::rename(gwf_gps_field, Field_code = "FIELD_ID")
head(gwf_gps_field)

# Add farm fix data
gwf_farm <- dplyr::left_join(gwf_gps_field, farm, by = "Field_code")

gwf_farm_total <- gwf_farm %>% group_by(Date_ID) %>% summarize (mode_hld_farm =names(which.max(table(hld_farm))))

# Bind data frame summary with original gpg points data frame
gwf_joined_df <- dplyr::full_join(gwf_farm, gwf_farm_total, by = c("Date_ID"))
gwf_shoot_day <- dplyr::select(gwf_joined_df, -c("X", "point.ID", "poly.ID", "OID_T", "SHAPE_AREA", "SHAPE_LEN"))

# Convert Timestamp to posixct
gwf_shoot_day$Timestamp <- as.POSIXct(gwf_shoot_day$Timestamp)

# Remove GPS fixes outisde of Islay
# Filter out for Long and Lat within Islay area
dist_Islay <- filter(gwf_shoot_day, between(Latitude, 55.57, 55.94))

gwf_shoot_day <- filter(dist_Islay, between(Longitude, -6.5, -6.02))

# Split data by winter
gwf_12 <- dplyr::filter(gwf_shoot_day, Winter == 2012)
gwf_13 <- dplyr::filter(gwf_shoot_day, Winter == 2013)
gwf_14 <- dplyr::filter(gwf_shoot_day, Winter == 2014)
gwf_15 <- dplyr::filter(gwf_shoot_day, Winter == 2015)
gwf_16 <- dplyr::filter(gwf_shoot_day, Winter == 2016)
gwf_17 <- dplyr::filter(gwf_shoot_day, Winter == 2017)
gwf_18 <- dplyr::filter(gwf_shoot_day, Winter == 2018)
gwf_19 <- dplyr::filter(gwf_shoot_day, Winter == 2019)
gwf_20 <- dplyr::filter(gwf_shoot_day, Winter == 2020)
gwf_21 <- dplyr::filter(gwf_shoot_day, Winter == 2021)

## 2.3 Sub-sample for consistent sampling rate ## ----
# to calculate daily distance traveled without being influenced by variable sampling rate
# NB! Using sub-sampling code from Liam Langley

trackSubSamp <- function(TD, dt=1, unit='hours'){
  
  ## order the timestamps
  TD <- TD[order(TD$Timestamp),] 
  
  # breakdown to datasets per bird
  unid = unique(TD$Device_ID) # the unique tag IDs in your data set
  nrid = length(unid)
  TDall = list(nrid)  
  TDred = list(nrid)
  timestep = paste(dt,unit)
  
  # create time sequence from min to max time with step sizes that are defined at the start of the function
  dati_start = min(TD$Timestamp,na.rm=T)
  dati_end = max(TD$Timestamp,na.rm=T)
  datiseq = seq(from=dati_start, to=dati_end, by=timestep)
  
  for (i in 1:nrid)
  {
    Dtemp = TD[TD$Device_ID == unid[i],]
    idx = sapply(datiseq, function( x) which.min( abs( difftime( Dtemp$Timestamp, x, units='mins') ) ) ) # finds closest time in data to your created time series
    TDall[[i]] = Dtemp
    TDred[[i]] = unique( Dtemp[idx,] ) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
  }
  
  TDred # return this list
}

# Sub-sample the data for each winter separately
list_gwf_12 <- trackSubSamp(gwf_12)
list_gwf_13 <- trackSubSamp(gwf_13)
list_gwf_14 <- trackSubSamp(gwf_14)
list_gwf_15 <- trackSubSamp(gwf_15)
list_gwf_16 <- trackSubSamp(gwf_16)
list_gwf_17 <- trackSubSamp(gwf_17)
list_gwf_18 <- trackSubSamp(gwf_18)
list_gwf_19 <- trackSubSamp(gwf_19)
list_gwf_20 <- trackSubSamp(gwf_20)
list_gwf_21 <- trackSubSamp(gwf_21)

# Convert list to data frame
df_gwf_12_ss <- do.call("rbind", list_gwf_12)

df_gwf_13_ss <- do.call("rbind", list_gwf_13)

df_gwf_14_ss <- do.call("rbind", list_gwf_14)

df_gwf_15_ss <- do.call("rbind", list_gwf_15)

df_gwf_16_ss <- do.call("rbind", list_gwf_16)

df_gwf_17_ss <- do.call("rbind", list_gwf_17)

df_gwf_18_ss <- do.call("rbind", list_gwf_18)

df_gwf_19_ss <- do.call("rbind", list_gwf_19)

df_gwf_20_ss <- do.call("rbind", list_gwf_20)

df_gwf_21_ss <- do.call("rbind", list_gwf_21)

# Combine the data frames
gwf_ss <- rbind(df_gwf_12_ss, df_gwf_13_ss, df_gwf_14_ss, df_gwf_15_ss, df_gwf_16_ss, df_gwf_17_ss, df_gwf_18_ss, df_gwf_19_ss, df_gwf_20_ss, df_gwf_21_ss)
str(gwf_ss)

## 2.4 Identify goose days ## ----

sun_dates2 <- seq.Date(as.Date(min(gwf_ss$Timestamp)), as.Date(max(gwf_ss$Timestamp)), by = 1)
sunt2 <-
  getSunlightTimes(
    date = sun_dates2,
    keep = c("dusk", "dawn"),
    lat = 55.738,
    lon = 6.246,
    tz = "UTC"
  )

## Add a column for the date
gwf_ss$date <- as.Date(gwf_ss$Timestamp)

## Now join on the sunset and sunrise times
sunt$lat <- NULL
sunt$lon <- NULL
gwf_ss2 <- left_join(gwf_ss, sunt2,  by = "date")
gwf_ss2$DayNight <- ifelse(gwf_ss2$Timestamp < gwf_ss2$dawn | gwf_ss2$Timestamp > gwf_ss2$dusk,
                           "night", "day")
head(gwf_ss2)
gwf_ss3 <- dplyr::filter(gwf_ss2, DayNight == "day")
gwf_ss3

# check number of fixes per day
gwf_ss3$sum <- 1

# Summary table 
# Use Date_ID to summarise, each ID = date and individual

gwf_ss_summary <- with(gwf_ss3, aggregate(sum, by = list(Date_ID), "sum"))

names(gwf_ss_summary)[1] <- "Date_ID"
names(gwf_ss_summary)[2] <- "Winter"
names(gwf_ss_summary)[2] <- "Total_fixes"

# Bind data frame summary with original gpg points data frame
gwf_joined_df <- dplyr::full_join(gwf_ss3, gwf_ss_summary, by = c("Date_ID"))
head(gwf_joined_df)

# Subset the original gps points based on the number of fixes
# NB! for now this is set as 5, look at what else could use/need for this later # Updated 24/01/2022
min(gwf_joined_df$Total_fixes)
max(gwf_joined_df$Total_fixes)

gwf_df_points_cleaned <- dplyr::filter(gwf_joined_df, gwf_joined_df$Total_fixes >= 5)
gwf_df_points_cleaned

gwf_df_points_cleaned$Species <- "GWFG"
min(gwf_df_points_cleaned$Total_fixes)
max(gwf_df_points_cleaned$Total_fixes)




#### 3. Calculate daily distance metrics ####
# For both species gps fixes have been sub-sampled to 1 hr with a minimum of 5 fixes per individual per day
# Include data from October - April
## 3.1 Clean and combine species GPS data together: ----
gwf_df_points_cleaned2 <- dplyr::select(gwf_df_points_cleaned, -c("lon", "lat"))
colnames(gwf_df_points_cleaned2)

gbg_df_points_cleaned$Timestamp <- as.POSIXct(gbg_df_points_cleaned$Timestamp)
gwf_df_points_cleaned$Timestamp <- as.POSIXct(gwf_df_points_cleaned$Timestamp)

# Streamline data tables to match columns
gbg_df_points_cleaned <- dplyr::select(gbg_df_points_cleaned, -c(Darvic, Sex, Site, Tag))

# Add year_day to GBG
gbg_df_points_cleaned$year_day <- yday(gbg_df_points_cleaned$Timestamp)

gbg_df_points_cleaned$UTC_Date <- as.POSIXct(gbg_df_points_cleaned$UTC_Date, format = "%d/%m/%Y")
gwf_df_points_cleaned$UTC_Date <- as.POSIXct(gwf_df_points_cleaned$UTC_Date, format = "%Y-%m-%d")

colnames(gbg_df_points_cleaned)
colnames(gwf_df_points_cleaned)
gbg_df_points_cleaned <- dplyr::rename(gbg_df_points_cleaned, Field_code = "FIELD_ID")
gwf_df_points_cleaned <- dplyr::rename(gwf_df_points_cleaned, FIELD_ID = "Field_code")
gbg_gwf_1hr <- rbind(gbg_df_points_cleaned, gwf_df_points_cleaned2)
head(gbg_gwf_1hr)

# Change NA in Shoot_field to "No"
gbg_gwf_1hr$Shoot_field[is.na(gbg_gwf_1hr$Shoot_field)] <- "No"

# Ensure timestamp in correct format
gbg_gwf_1hr$Timestamp <- as.POSIXct(gbg_gwf_1hr$Timestamp, format = "%Y-%m-%d %H:%M:%S")

## 3.2 Calculate cumulative total distance ## ----
# Between sequential GPS points by day 
# Note Use unique ID for each individual for each day

# Create a data frame to store the trip metrics
bird_trips <- unique(gbg_gwf_1hr$Date_ID) # Identify unique bird trips (Date_ID)
metrics <- as.data.frame(bird_trips) # Create blank data frame of unique birds trips to store metrics in later

# Create blank columns within data frame to fill
# Note! Not essential for the script to run, but if do not do the columns will fill with value of first trip, makes it hard to spot errors

# Bird info
metrics$birdID <- NA

# Total distance travelled
metrics$totaldist <- NA
metrics$daydist <- NA

# Check data without running loop j=1

for (j in 1:length(bird_trips)){ # Loop through each trip at a site
  trip <- subset(gbg_gwf_1hr,gbg_gwf_1hr$Date_ID == bird_trips[j])
  
  # Store the additional data about each trip in the metrics dataframe
  metrics$birdID[j]    <- as.character(trip$Device_ID[1])
  
  # Calculate the total distance travelled and maximum distance from the colony (km)
  trip$dist[2:nrow(trip)] <- pointDistance(
    matrix(c(trip$Longitude[1:(nrow(trip)-1)], trip$Latitude[1:(nrow(trip)-1)]),ncol=2),
    matrix(c(trip$Longitude[2:nrow(trip)], trip$Latitude[2:nrow(trip)]),ncol=2),longlat=TRUE)/1000
  metrics$totaldist[j] <- sum(trip$dist, na.rm = T)
  
  print(length(bird_trips)-j)
}

# Add other data into data frame
day_data <- dplyr::select(gbg_gwf_1hr, Device_ID, UTC_Date, Species, Month, Year, Winter, Date_ID, spatial_overlap, Total_fixes, year_day, fix_gps, Raw_Farm_Name, Field_code, Shoot_field, mode_hld_farm)
day_data <- distinct(day_data)
day_data <- dplyr::rename(day_data, bird_trips = Date_ID)

trip_data_full <- dplyr::left_join(metrics, day_data, by = "bird_trips")



#### 4.  Add winter shooting metrics to data #### 
# Note: Include 2 metrics of winter shooting disturbance:
# 1. Total number of shooting events in a winter
# 2. Frequency of shooting events - number of events/day

## 4.1 Add and summarise shooting data ## ----
# Load in shooting data
Logs <- read.csv("Shooting logs/All_logs_cleaned.csv")
Logs$Date_Cl <- as.POSIXct(Logs$Date_Cl)

# Add column for winter
Logs <- Logs %>% mutate(Month = month(Date_Cl), Year = year(Date_Cl))

Logs_total <- Logs %>% mutate(Logs, Winter = ifelse(Month >5, Year, Logs$Year-1))

# Tally number of shooting events per winter in new dataframe
# Note Group by winter

# Add column for shooting events
Logs_total <- mutate(Logs_total, shooting_event = ifelse(Shots_fired_Cl >= 1, 1, 0))

# Change NA into zero 
Logs_total$shooting_event[is.na(Logs_total$shooting_event)] <- 0
Logs_total

winter_shoot_total <- Logs_total %>% group_by(Winter) %>% tally(shooting_event)
winter_shoot_total <- as.data.frame(winter_shoot_total)
winter_shoot_total <- dplyr::rename(winter_shoot_total, events_total = n)

# Calculate number of events per day 
# NB! Number of shooting events/number of shooting days within shooting season - taken as 1st Nov to 1st April
winter_shoot_total2 <- mutate(winter_shoot_total, shoot_freq = events_total/180)
winter_shoot_total2$Winter <- as.factor(winter_shoot_total2$Winter)

# Filter out data to within shooting period (i.e from 1st November)
dist_data_islay <- dplyr::filter(trip_data_full, Month != 10)
unique(dist_data_islay$Month)
dist_data_islay <- dplyr::rename(dist_data_islay, Disturbance = spatial_overlap)

# Add day since start of shooting winter
dist_islay_shoot <- mutate(dist_data_islay, shoot_date = ifelse(dist_data_islay$year_day >= 305, dist_data_islay$year_day - 305, dist_data_islay$year_day + 60))
dist_islay_shoot$totaldist_log <- log(dist_islay_shoot$totaldist)
dist_islay_shoot$Winter <- as.factor(dist_islay_shoot$Winter)

winter_shoot_total2$Winter <- as.factor(winter_shoot_total2$Winter)

islay_distance <- dplyr::left_join(dist_islay_shoot, winter_shoot_total2, by = "Winter")

islay_distance$Disturbance <- as.numeric(islay_distance$Disturbance)
indiv_winter_tot <- islay_distance %>% group_by(Winter, Device_ID) %>% tally(Disturbance)
head(islay_distance)

# Calculate cumulative number of disturbance events for each individual for each winter
islay_distance1 <- dplyr::select(islay_distance, -c("Field_code"))
islay_distance3 <- distinct(islay_distance1)

islay_distance <- islay_distance3 %>% group_by(Winter, Device_ID) %>% mutate(indiv_dist = cumsum(Disturbance)) 
islay_distance <- as.data.frame(islay_distance)
max(islay_distance$indiv_dist)

largest_dist <- dplyr::filter(islay_distance, indiv_dist == 43)
largest_dist
largest_dsit2 <- dplyr::filter(islay_distance, Device_ID == "BLO27" & Winter == "2014")

# recode Disturbance to disturbed and undisturbed
islay_distance$Disturbance2 <- islay_distance$Disturbance
islay_distance$Disturbance2 <- dplyr::recode(islay_distance$Disturbance2, "0" = "No_Shoot", "1" = "Shoot")

# Note: Limit to >8 fixes per day 
islay_distance <- dplyr::filter(islay_distance, Total_fixes >= 8)

## 4.2 Format variables ## ----
# Scale continuous variables
islay_distance$shoot_date_sc <- as.numeric(scale(islay_distance$shoot_date))
islay_distance$events_total_sc <- as.numeric(scale(islay_distance$events_total))
islay_distance$shoot_freq_sc <- as.numeric(scale(islay_distance$shoot_freq))
islay_distance$indiv_dist_sc <- as.numeric(scale(islay_distance$indiv_dist))

# Ensure categorical variables are factors
islay_distance$Species <- as.factor(islay_distance$Species)
islay_distance$Device_ID <- as.factor(islay_distance$Device_ID)
islay_distance$Disturbance2 <- as.factor(islay_distance$Disturbance2)
islay_distance$Winter <- as.factor(islay_distance$Winter)
islay_distance$Shoot_field <- as.factor(islay_distance$Shoot_field)
non_shot <- dplyr::filter(islay_distance, Disturbance2 == "No_Shoot")
shot <- dplyr::filter(islay_distance, Disturbance2 == "Shoot")

str(islay_distance)

## 4.3 Add Sex to data ## ----
min(islay_distance$Total_fixes)
islay_distance_sub <- dplyr::filter(islay_distance, Total_fixes >= 8)

sex <- read.csv("tracking data/GPS_TAG_SEX.csv")
islay_distance_sub2 <- dplyr::left_join(islay_distance_sub, sex, by = "Device_ID")

## 4.4 Prep data for model  ## ----
islay_distance_sub2 <- islay_distance_sub2 %>% mutate(
  Tag_Winter = paste(Device_ID, Winter, sep = "_")
)

islay_distance_sub2 <- dplyr::mutate(islay_distance_sub2, Disturbance_shoot = ifelse(Shoot_field == "Yes", "shoot_dist", islay_distance_sub2$Disturbance2)
)

#write.csv(islay_distance_sub2, "Foraging_Distance_Model_Data.csv") # Data to use specifically tp run models




# ------------------------------------ 
####      Daily Distance Models     ####
# ------------------------------------ 

#### 1. Read in and prep model data ####
islay_distance_sub2 <- read.csv("Derived data/Foraging_Distance_Model_Data.csv")
islay_distance_sub2$shoot_date <- as.factor(islay_distance_sub2$shoot_date)
islay_distance_sub2$Tag_Winter <- as.factor(islay_distance_sub2$Tag_Winter)
islay_distance_sub2$group_id <- as.numeric(islay_distance_sub2$Tag_Winter)
islay_distance_sub2$Winter <- as.factor(islay_distance_sub2$Winter)





#### 2. Models with 2-level disturbance classification ####
# re-level factor to GWfG as base line
#islay_distance_sub2$Species <- relevel(as.factor(islay_distance_sub2$Species),"GBG")
dist_glob <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                       Disturbance2*Species + Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) + 
                       Disturbance2*Species*scale(indiv_dist) + ar1(shoot_date + 0|group_id) +
                       (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                     data = islay_distance_sub2,
                     family = Gamma("log"), na.action = "na.fail")

dist_globb <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) + 
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globc <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Disturbance2*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globd <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globe <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globf <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"))

dist_globg <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globh <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globi <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globj <- glmmTMB(totaldist ~ Disturbance2 + Species + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globk <- glmmTMB(totaldist ~ Species + Sex + scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

dist_globl <- glmmTMB(totaldist ~ Species + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")





#### 3. Models with 3-level disturbance classification ####
shoot_glob <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                        Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) + 
                        Disturbance_shoot*Species*scale(indiv_dist) + ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = islay_distance_sub2,
                      family = Gamma("log"), na.action = "na.fail")

shoot_globb <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) + 
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globc <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globd <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globe <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globf <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"))

shoot_globg <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globh <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globi <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globj <- glmmTMB(totaldist ~ Disturbance_shoot + Species + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globk <- glmmTMB(totaldist ~ Species + Sex + scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globl <- glmmTMB(totaldist ~ Species + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = islay_distance_sub2,
                       family = Gamma("log"), na.action = "na.fail")





#### 4. Model Checks ####
## 4.1 Calculate and compare models with AICc ## ----
ICtab(dist_glob,dist_globb,dist_globc,dist_globd,dist_globe,dist_globf,dist_globg,dist_globh,
      dist_globi, dist_globj, dist_globk, dist_globl,shoot_glob,
      shoot_globb,shoot_globc,shoot_globd,shoot_globe,shoot_globf,
      shoot_globg,shoot_globh,shoot_globi, shoot_globj, shoot_globk, shoot_globl,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

## 4.2 Dhama model checks  ## ----
# Simulate residuals and plot - this is for best model can check for all models

simtop1 <- simulateResiduals(dist_globf)
plot(simtop1)

## 4.3 Calculate r-squared ## ----
r.squaredGLMM(dist_glob)
r.squaredGLMM(dist_globb)
r.squaredGLMM(dist_globc)
r.squaredGLMM(dist_globd)
r.squaredGLMM(dist_globe)
r.squaredGLMM(dist_globf)
r.squaredGLMM(dist_globg)
r.squaredGLMM(dist_globh)
r.squaredGLMM(dist_globi)
r.squaredGLMM(dist_globj)
r.squaredGLMM(dist_globk)
r.squaredGLMM(dist_globl)




#### 5. Model Plots ####
## 5.1 Extract effects from models ## ----
df_dist_top <- as.data.frame(effect("Disturbance2:Species",dist_globf)) # Specific interaction terms
effects_top <- as.data.frame(allEffects(dist_globd)) # All effects

#df_dist_top2 <- emmeans(dist_globf, specs = ~ Disturbance2:Species, by = "Species", type = "response")
#df_dist_top2
#summary(dist_globf)
confint(dist_globf) # Confidence intervals of top model

# Visualise all effects in top model
plot(allEffects(dist_globf))

## 5.2 Plot overall top model ## ----
# Top model = dist_globf -> disturbance model with two-level factor of disturbance
# Plots of main model outputs of the effect of disturbance on total daily distance traveled
Disturbance_labels <- c("No Shooting", "Shooting")
pd <- position_dodge(0.43)

disturbance_plot_topmod <- ggplot(df_dist_top, aes(Disturbance2, fit)) +
  geom_point(aes(colour = Species), position = pd, size = 3) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Species), width=.4, position = pd) +
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(x = " Disturbance Day", y = "Daily distance travelled (km)", colour = "Species")+
  scale_x_discrete(labels = Disturbance_labels) +
  scale_y_continuous(limits = c(1,6), breaks = seq(1,6, by = 0.5)) +
  theme_classic() +
  theme(text = element_text(size = global_size))
disturbance_plot_topmod

ggsave("Plots/Script 5) plots/Daily travel distance - 2level disturbance - top model.jpeg", disturbance_plot_topmod)

## 5.3 Plot model of 3-level disturbance factor ## ----
# supplementary material
Disturbance_labels <- c("No Shooting", "Shooting nearby", "Shooting in field")
pd <- position_dodge(0.43)

field_disturbance_plot <- ggplot(emmean_dist, aes(Disturbance_shoot, response)) +
  geom_point(aes(colour = Species), position = pd, size = 3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, color = Species), width=.4, position = pd) +
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(x = " Disturbance Day", y = "Daily distance travelled (km)", colour = "Species")+
  scale_x_discrete(labels = Disturbance_labels) +
  scale_y_continuous(limits = c(1,6), breaks = seq(1,6, by = 0.5)) +
  theme_classic()
field_disturbance_plot


ggsave("Plots/Script 5) plots/Daily travel distance - 3level disturbance - top.jpeg", field_disturbance_plot)



#### End ####