## Luke Ozsanlav-Harris
## 21/04/22

## Compare behavior days that are ans aren't shot at using classified acc data from Ornitela tags

## Packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(svMisc)
library(suncalc)

library(geosphere)
library(sf)
library(amt)

library(lme4)
library(DHARMa)
library(MuMIn)
library(effects)
library(performance)
library(glmmTMB)



## TO DO:
## Properly filter the ACC data so no burst from off Islay are used
## Filter out to only have birds that were on Islay for a certain number of days



##
#### 1. Read in data sets ####
##


##
#### 1.1 Read in classified ACC data ####
##

## Data Cleaning, extracting just the GPS data
orn_behav <- readRDS("tracking data/Classified_ACC.RDS")
table(orn_behav$device_id)

## read in meta data so i can just keep the Scotland birds
Meta <- fread("MetaData/Tagged bird summary data new.csv")
Meta <- filter(Meta, Ringing.location %in% c("ISLA", "WEST", "LOKE", "DYFI", "KINT") | S.N == 17795)

## filter the gps data for the device ids left
orn_behav <-filter(orn_behav, device_id %in% unique(Meta$S.N))

## parse timestamp
orn_behav$UTC_datetime <- ymd_hms(orn_behav$UTC_datetime)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind <- orn_behav %>% 
       dplyr::select(UTC_datetime, device_id) %>%
       duplicated()
sum(ind)
orn_behav <- orn_behav %>% filter(!ind)


## Streamline the dataset for a first pass, first just have subset of the columns
setnames(orn_behav, old = "device_id", new = "Tag_ID")

## Now filter for the winter, create a year day column first
orn_behav$year_day <- yday(orn_behav$UTC_datetime)
winter_orn_behav <- filter(orn_behav, year_day >= 305 | year_day <= 92)
table(winter_orn_behav$Tag_ID)




##
#### 1.2 Read in Islay winter GPS data ####
##

GPS_data <- readRDS("Derived data/All_winter_GPS_with_habitat.RDS")

## only use the GPS data from birds with ACC data
ACC_tags <- unique(winter_orn_behav$Tag_ID)
GPS_set <- filter(GPS_data, Tag_ID %in% ACC_tags)

## give each fix its own unique ID number
GPS_set$ID <- 1:nrow(GPS_set)





##                                                 
#### 1.3 Read in shooting data and field boundaries ####
##                                                 


## read in the shooting logs
Logs <- fread("Shooting logs/All_logs_cleaned.csv")

## create a timestamp column for the shooting logs 
## think the ones that fail to parse are the ones were the time is missing and it is just a date
Logs$timestamp <- ymd_hms(paste0(Logs$Date_Cl, " ", Logs$Time_Cl))

## Now filter the logs for the same time period as the tracking data
winter_logs <- filter(Logs, timestamp > min(GPS_data$UTC_datetime) & timestamp < max(GPS_data$UTC_datetime))

## Now filter out just the shooting events
winter_shots <- filter(winter_logs, Shots_fired_Cl > 0)



## read in the shapefiles of the Islay field boundaries
Fields <- st_read(dsn = "Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY",
                  layer = "AS_ISLAY_GMS_FIELD_BOUNDARY")
st_crs(Fields)

## remove any duplicates in the field data
ind3 <- duplicated(Fields$FIELD_ID)
sum(ind3) 
Fields <- Fields %>% filter(!ind3)




##                                                   
#### 2. Combine shooting data and spatial field data ####
##                                                   


## calculate the centroid of each field
Field_centres <- st_centroid(Fields)

## now create a data set with the shooting events but where each one has the center for that field
setnames(winter_shots, old = "Field_code_Cl", new = "FIELD_ID")
winter_shots_cent <- left_join(winter_shots, Field_centres, by = "FIELD_ID")

## for now filter out the shooting events were we don't know the field code
winter_shots_cent2 <- filter(winter_shots_cent, is.na(SHAPE_LEN)==F)

## Give each shooting event a unique ID number
winter_shots_cent2$shot_ID <- 1:nrow(winter_shots_cent2)

## Make winter_shots_cent2 an sf object again
winter_shots_cent_sf <- st_as_sf(winter_shots_cent2)
winter_shots_cent_sf <- st_transform(winter_shots_cent_sf, crs = st_crs(Field_centres))







##
#### 3. Subset ACC data to just the extent of Islay ####
##

## The GPS_set data set contains just GPS data from the winter on Islay
## Therefore if i restrict the ACC data to the same tag periods
## then i can constrain the ACC data to just come from Islay winter

## Add a tag date column to both data sets
GPS_set$tag_date <- paste0(GPS_set$Tag_ID, "_", as.Date(GPS_set$UTC_datetime))
winter_orn_behav$tag_date <- paste0(winter_orn_behav$Tag_ID, "_", as.Date(winter_orn_behav$UTC_datetime))

## Now filter out only the tag dates in the GPS data set
UnQs <- unique(GPS_set$tag_date)
Islay_orn_behav <- filter(winter_orn_behav, tag_date %in% UnQs)






##                                                   
#### 4. Label fixes associated with each shooting event ####
##                                                     


## Define function that converts timestamp immediately before shooting event the smallest value
## y needs to be the timestamps from the tracking data and x is just a single value from the shooting data sequence
difftime2 <- function(x=x, y=y){
  
  time_difs <- difftime(y, x, units='mins')*-1
  time_difs[time_difs<0] <- 1000000000
  return(as.vector(time_difs))
}

## Define function that will loop through each tag and find the nearest fix before each shooting event while the tag was active
## The fix prior to the shooting event has to be within a certain spatial and temporal buffer to the shooting event
## TD = tracking data. Need tag ID (Tag_ID), lat/long (Longitude/Latitude), timestamp (UTC_datetime) and unique ID for each for row (ID)
## SD = the shooting data. Need the time of the shot (timestamp), the field center (as an sf geometry) and unique ID for each for row (shot_ID)
## time_thresh = the max time in minutes that a fix can be before a shooting event
## dist_thresh = the max distance in meters the fix before a shooting event can be from the location of the shooting event
## field_crs = an sf object that has the crs of the field centers in SD

Shoot_proximity <- function(TD = TD, SD = SD, time_tresh = time_tresh, dist_thresh = dist_thresh, field_crs = field_crs){
  
  ##  make sure that the data set is ordered by time
  TD <- TD[order(TD$UTC_datetime),] 
  
  ## breakdown to datasets per bird
  unid <- unique(TD$Tag_ID) # list of unique tag IDs
  nrid <- length(unid) # number of unique birds
  
  ## Create lists to put each tag ID in
  Prox_list <- vector(mode = "list", length = nrid)
  
  ## Loop through each tag
  for (i in 1:nrid) { ## ** CHANGE 15 BACK TO nrid **###
    
    message("Tag ", i, " out of ", nrid)
    
    
    ## DATA PREP ##
    ## filter out each tag one at a time
    Dtemp <- TD[TD$Tag_ID == unid[i],]
    
    ## create sequence of shooting data while the tag was active
    Start_time <- min(Dtemp$UTC_datetime)
    End_time <- max(Dtemp$UTC_datetime)
    SD_sub <- filter(SD , timestamp > Start_time, timestamp < End_time)
    shootseq <- SD_sub$timestamp
    
    if(length(shootseq)==0){next}
    
    
    ## FINDING NEAREST PRIOR FIX ##
    # find the first fix before the shooting event
    idx = sapply(shootseq, function(x) which.min(difftime2(y = Dtemp$UTC_datetime, x = x)))
    
    ## combine these fixes immediately before the shooting event with the shooting data
    aa <- Dtemp[idx,]
    bb <- cbind(aa , SD_sub)
    
    ##** COULD POTENTIALLY TAKE THIS SECTION OUT OF THE LOOP **##  
    
    ## CALCUALTE TIME DIFFERENCE AND DISTANCE ##
    ## calculate time differences between fix and shooting event
    bb$time_dif <- as.duration(ymd_hms(bb$timestamp) %--% ymd_hms(bb$UTC_datetime))
    bb$time_dif <- as.numeric(bb$time_dif, "minutes")
    
    ## Now calculate distances between the fix and the shooting event
    ## Make the fixes an sf object
    gps_sf <- st_as_sf(aa, coords = c("Longitude", "Latitude"), 
                       crs = 4326, agr = "constant")
    ## now convert to the same crs as the field data
    gps_sf <- st_transform(gps_sf, crs = st_crs(field_crs))
    
    ## make the shooting events an sf object
    shootseq_sf <- st_as_sf(SD_sub)
    shootseq_sf <- st_transform(shootseq_sf, crs = st_crs(field_crs))
    
    ## calculate the distance between the two sets of points
    distances <- st_distance(gps_sf, shootseq_sf, by_element = TRUE)
    
    ## put the distances back into the shooting proximity data set
    bb$dist_to_shot <- as.numeric(distances)
    
    
    ## RULES TO KEEP FIXES IN SPATIO TEMPORAL BUFFER ##
    ## Remove any fixes that were greater then 120 before the shooting event or were after the shooting event
    bb2 <- filter(bb, time_dif <= 0, time_dif >= -time_tresh)
    
    ## Remove any fixes that were greater than a certain number of meters from the shooting event
    bb3 <- filter(bb2, dist_to_shot < dist_thresh)
    
    ## Put this data set into a list that was initialized outside of the loop
    Prox_list[[i]] <- bb3
    
  }
  
  Prox_list
  
}



## run the function
## Used 2000m as the distance threshold, they ony tested up to c1,500m in the Islay study and at this distance the effect became small
## Added c500m on to allow a bit of wiggle room
system.time(Prox_list <- Shoot_proximity(TD = GPS_set, SD = winter_shots_cent2, time_tresh = 60, dist_thresh = 797, field_crs = Field_centres))

## bind the lists together into one data frame
All_prox <- plyr::ldply(Prox_list)








##                                                                     
#### 5. Summarize on what days and how many times birds were shot at ####
##                                                                     

## create a date column for the grouping
All_prox$date <- as.Date(All_prox$UTC_datetime)

## create a dummy column to sum rows
All_prox$dummy <- 1

## now summarise the  number of shooting events by each tag day
DistSum <- All_prox %>% 
           group_by(date, Tag_ID) %>% 
           summarise(n_shot = sum(dummy))





##
#### 6. Summarize daily behavior ####
##


## First label times as either day or night
## "sunrise" and "sunset" gave short days, i.e. the geese wouls still be out on daytime feeding sites
## Either use  "dusk", "dawn" or "nauticalDusk", "nauticalDawn" (nautical versions give the longer day length)
sun_dates <- seq.Date(as.Date(min(Islay_orn_behav$UTC_datetime)), as.Date(max(Islay_orn_behav$UTC_datetime)), by = 1)
sunt <-
  getSunlightTimes(
    date = sun_dates,
    keep = c("dusk", "dawn"),
    lat = 55.738,
    lon = 6.246,
    tz = "UTC"
  )


## Add a column for the date
Islay_orn_behav$date <- as.Date(Islay_orn_behav$UTC_datetime)

## Now join on the sunet and sunrise times
sunt$lat <- NULL
sunt$lon <- NULL
Islay_orn_behav2 <- left_join(Islay_orn_behav, sunt,  by = "date")
Islay_orn_behav2$DayNight <- ifelse(Islay_orn_behav2$UTC_datetime < Islay_orn_behav2$dawn | Islay_orn_behav2$UTC_datetime > Islay_orn_behav2$dusk,
                                     "night", "day")

## Add binary columns that indicate which behavior
Islay_orn_behav2$Flight <- ifelse(Islay_orn_behav2$Behaviour == "flight", 1, 0)
Islay_orn_behav2$Walk <- ifelse(Islay_orn_behav2$Behaviour == "walk", 1, 0)
Islay_orn_behav2$Graze <- ifelse(Islay_orn_behav2$Behaviour == "graze", 1, 0)
Islay_orn_behav2$Stat <- ifelse(Islay_orn_behav2$Behaviour == "stationary", 1, 0)
Islay_orn_behav2$Active <- ifelse(Islay_orn_behav2$Behaviour == "graze" | Islay_orn_behav2$Behaviour == "walk", 1, 0)



## Create summary of behaviors before and after shooting
BehSum <- Islay_orn_behav2 %>% 
          group_by(date, Tag_ID, DayNight) %>% 
          summarise(Flight = sum(Flight), 
                    Walk = sum(Walk),
                    Graze = sum(Graze),
                    Stat = sum(Stat),
                    Active = sum(Active),
                    Total = length(Behaviour))



##
#### 7. Add on number of shooing events to the dataset ####
##


## join number of shots and daily behavioral summary
BehSum$Tag_ID <- as.character(BehSum$Tag_ID)
BehSum2 <- left_join(BehSum, DistSum, by = c("date", "Tag_ID"))
stopifnot(nrow(BehSum2)==nrow(BehSum2))


## change all the NAs in the n_shot column to 0s
BehSum2$n_shot <- ifelse(is.na(BehSum2$n_shot)==T, 0, BehSum2$n_shot)







##
#### 8. Model the proportion of time spent of certain behaviors varies ####
##


## Add extra columns to put in to the model
BehSum2$year <- year(BehSum2$date) # year column
BehSum2$year_day <- yday(BehSum2$date) 
BehSum2$from_nov1 <- ifelse(BehSum2$year_day > 250, 
                            BehSum2$year_day - 305, 
                            BehSum2$year_day + 60) # day since November 1st
BehSum2$shot <- ifelse(BehSum2$n_shot >  0, "Shot", "Not_Shot") # whether there was shooting or not
BehSum2$winter <- ifelse(BehSum2$year_day < 150, paste0((BehSum2$year-1), "-", BehSum2$year), paste0(BehSum2$year, "-", (BehSum2$year+1)))
BehSum2$tag_winter <- paste0(BehSum2$Tag_ID, "_", BehSum2$winter) # tag winter identifier
BehSum2$prev_shot <- ifelse(lag(BehSum2$shot, n = 2L) == "Shot", "Prev_shot", "Not")
BehSum2$prev_shot <- ifelse(is.na(BehSum2$prev_shot)==T, "Not", BehSum2$prev_shot)


## Create indicator columns, e.g. fro not grazing, not walking etc
BehSum2$Not_Graze <- BehSum2$Total - BehSum2$Graze
BehSum2$Not_Walk <- BehSum2$Total - BehSum2$Walk
BehSum2$Not_Flight <- BehSum2$Total - BehSum2$Flight
BehSum2$Not_Stat <- BehSum2$Total - BehSum2$Stat


## set variables as correct class
BehSum2$year <- as.factor(BehSum2$year)
BehSum2$from_nov1 <- as.numeric(BehSum2$from_nov1)
BehSum2$DayNight <- as.factor(BehSum2$DayNight)
BehSum2$Tag_ID <- as.factor(BehSum2$Tag_ID)
BehSum2$shot <- as.factor(BehSum2$shot)
BehSum2$from_nov1Factor <- as.factor(BehSum2$from_nov1)


BehSum2$obs <- 1:nrow(BehSum2)
BehSum3 <- filter(BehSum2, DayNight == "day")




##
#### 8.0 Minimum data requirement checks ####
##

## remove data points based on whether they had a certain number of data points per day
BehSum2 <- filter(BehSum2, Total >= 6) ##**CHANGED THIS FROM 10 TO 6 (12/12/22). Now aligned with the daily ODBA models


## summarize number of data point for each tag winter for day and night
BehSum2$dummy <- 1
ID_sum <- BehSum2 %>% 
          group_by(Tag_ID, winter) %>% 
          summarise(n_days = sum(dummy), 
                    minNov =min(from_nov1), 
                    maxNov1 =max(from_nov1))

## now filter out instances with 5 or less days
ID_sum <- filter(ID_sum, n_days > 5)
ID_sum$n_days <- NULL

## Use inner join to remove the tag winters with 5 or less days
BehSum3 <- inner_join(BehSum2, ID_sum, by = c("Tag_ID", "winter"))
## remove the first and last day fo the winter for each inivdiual as these likely contain periods of migration
BehSum3 <- filter(BehSum3, !from_nov1 == minNov)
BehSum3 <- filter(BehSum3, !from_nov1 == maxNov1)








##
#### 8.1 Walk model ####
##

## Shoot term not significant, model checks OKAY

hist(sqrt(BehSum2$Walk))
Walk_mod <- glmmTMB(formula = cbind(Walk, Not_Walk) ~ shot*DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                     data = BehSum3,
                     family = betabinomial)


summary(Walk_mod)
drop1(Walk_mod, test = "Chi")

qqnorm(residuals(Walk_mod))

## Check model performance and assumptions
#model_performance(Walk_mod)
#check_model(Walk_mod)


## Model checks with DHARMa

##---- use DHARMa to get QQ plot and resids vs fitted ----##
Resids_Walk <- simulateResiduals(Walk_mod)
plot(Resids_Walk)

## plot reiduals vs predicted values
par(mfrow=c(1,3))
plotResiduals(Resids_Walk[["scaledResiduals"]], form = BehSum2$DayNight)
plotResiduals(Resids_Walk[["scaledResiduals"]], form = BehSum2$shot)
plotResiduals(Resids_Walk[["scaledResiduals"]], form = BehSum2$from_nov1)

## Check dispersion
testDispersion(Resids_Walk)
testZeroInflation(Resids_Walk)

## check autocorrelation
New_resids <- recalculateResiduals(Resids_Walk, group = BehSum2$from_nov1)
testTemporalAutocorrelation(Resids_Walk, time = unique(BehSum2$from_nov1))


## create all candidate models using dredge (trace shows progress bar)
dredge_set <- MuMIn::dredge(Walk_mod, trace = 2)
nested_set <- subset(dredge_set, !MuMIn::nested(.), recalc.weights=T)
delta6_set <- subset(nested_set, delta<=6, recalc.weights=T)


## plot the model effects
top_mod_effects <- predictorEffects(Walk_mod)
plot(top_mod_effects)







##
#### 8.2 Stat model ####
##

## Compare models to determine effect of shooting variable
Stat_modInt <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ shot*DayNight + winter + poly(from_nov1,2) + 
                                                         (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                        data = BehSum3,
                        family = betabinomial)
Stat_modfix <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ shot + DayNight + winter + poly(from_nov1,2) + 
                                                         (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                       data = BehSum3,
                       family = betabinomial)
Stat_modNo <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ DayNight + winter + poly(from_nov1,2) + 
                                                        (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                       data = BehSum3,
                       family = betabinomial)


AIC(Stat_modInt, Stat_modfix, Stat_modNo)
summary(Stat_modNo)
drop1(Stat_mod, test = "Chi")
emmeans::emmeans(Stat_mod, ~shot*DayNight, type = "response")
confint(Stat_mod)


## Check model performance and assumptions
#model_performance(Stat_mod)
#check_model(Stat_mod)


## Model checks with DHARMa

##---- use DHARMa to get QQ plot and resids vs fitted ----##
Resids_Stat <- simulateResiduals(Stat_modInt)
plot(Resids_Stat)

## plot reiduals vs predicted values
par(mfrow=c(1,3))
plotResiduals(Resids_Stat[["scaledResiduals"]], form = BehSum2$DayNight)
plotResiduals(Resids_Stat[["scaledResiduals"]], form = BehSum2$shot)
plotResiduals(Resids_Stat[["scaledResiduals"]], form = BehSum2$from_nov1)

## Check dispersion
testDispersion(Resids_Stat)
testZeroInflation(Resids_Stat)

## check autocorrelation
New_resids <- recalculateResiduals(Resids_Stat, group = BehSum2$from_nov1)
testTemporalAutocorrelation(Resids_Stat, time = unique(BehSum2$from_nov1))


## create all candidate models using dredge (trace shows progress bar)
dredge_set <- MuMIn::dredge(Stat_mod, trace = 2)
nested_set <- subset(dredge_set, !MuMIn::nested(.), recalc.weights=T)
delta6_set <- subset(nested_set, delta<=6, recalc.weights=T)


## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effects <- predictorEffects(Stat_mod)
plot(top_mod_effects)
## now extract the fits for the first variable and bind them together
effects1 <- top_mod_effects["shot"]
fit1 <- as.data.frame(cbind(effects1[["shot"]][["fit"]], effects1[["shot"]][["lower"]], 
                            effects1[["shot"]][["upper"]], effects1[["shot"]][["x"]][["shot"]]))
## change the names to something meaningful
setnames(fit1, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "level"))
## transform variables back to proportion scale
fit1$fit <- boot::inv.logit(fit1$fit)
fit1$lower <- boot::inv.logit(fit1$lower)
fit1$upper <- boot::inv.logit(fit1$upper)
fit1$Shooting <- c("No Shooting", "Shooting", "No Shooting", "Shooting")
fit1$Time_of_Day <- c("Day", "Day", "Night", "Night")

## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit1, aes(x= Shooting, ymin = lower, ymax = upper, colour = Time_of_Day), width = 0.4, size = 0.8) +
  geom_point(data=fit1, aes(x= Shooting, y = fit, colour = Time_of_Day), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of bursts Stationary") +
  theme_bw() +
  scale_colour_manual(values=c("#000000", "#cc0000")) +
  labs(colour="Time of Day") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        #legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =12), 
        axis.title.x = element_text(size =14),
        axis.title.y = element_text(size =14), 
        strip.text.x = element_text(size =12))


## Save a plot
# ggsave("Plots/Script 4) plots/Proportion_Stationary.png", 
#        width = 19, height = 22, units = "cm")






##
#### 8.3 Flight model ####
##

## compare three models
hist(sqrt(BehSum2$Flight))
Flight_modInt <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ shot*DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                    data = BehSum3,
                    family = betabinomial)
Flight_modFix <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ shot + DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                         data = BehSum3,
                         family = betabinomial)
Flight_modNo <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                         data = BehSum3,
                         family = betabinomial)

AIC(Flight_modInt, Flight_modFix, Flight_modNo)
summary(Flight_modNo)
drop1(Flight_modNo, test = "Chi")

## Check model performance and assumptions
#model_performance(Flight_mod)
#check_model(Flight_mod)
#check_model(Flight_mod, check ="reqq") # echek the random effects strucutre


## Model checks with DHARMa

##---- use DHARMa to get QQ plot and resids vs fitted ----##
Resids_Flight <- simulateResiduals(Flight_mod)
plot(Resids_Flight)

## plot reiduals vs predicted values
par(mfrow=c(1,3))
plotResiduals(Resids_Flight[["scaledResiduals"]], form = BehSum2$DayNight)
plotResiduals(Resids_Flight[["scaledResiduals"]], form = BehSum2$shot)
plotResiduals(Resids_Flight[["scaledResiduals"]], form = BehSum2$tag_winter)
unique(BehSum2$tag_winter)
## Check dispersion
testDispersion(Resids_Flight)
testZeroInflation(Resids_Flight)

## check autocorrelation
New_resids <- recalculateResiduals(Resids_Flight, group = BehSum2$from_nov1)
testTemporalAutocorrelation(Resids_Flight, time = unique(BehSum2$from_nov1))


## create all candidate models using dredge (trace shows progress bar)
dredge_set <- MuMIn::dredge(Flight_mod, trace = 2)
nested_set <- subset(dredge_set, !MuMIn::nested(.), recalc.weights=T)
delta6_set <- subset(nested_set, delta<=6, recalc.weights=T)


## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effectsF <- predictorEffects(Flight_mod)
plot(top_mod_effectsF)
## now extract the fits for the first variable and bind them together
effects2 <- top_mod_effectsF["shot"]
fit2 <- as.data.frame(cbind(effects2[["shot"]][["fit"]], effects2[["shot"]][["lower"]], 
                            effects2[["shot"]][["upper"]], effects2[["shot"]][["x"]][["shot"]]))
## change the names to something meaningful
setnames(fit2, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "level"))
## transform variables back to proportion scale
fit2$fit <- boot::inv.logit(fit2$fit)
fit2$lower <- boot::inv.logit(fit2$lower)
fit2$upper <- boot::inv.logit(fit2$upper)
fit2$Shooting <- c("No Shooting", "Shooting", "No Shooting", "Shooting")
fit2$Time_of_Day <- c("Day", "Day", "Night", "Night")

## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit2, aes(x= Shooting, ymin = lower, ymax = upper, colour = Time_of_Day), width = 0.4, size = 0.8) +
  geom_point(data=fit2, aes(x= Shooting, y = fit, colour = Time_of_Day), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of bursts in Flight") +
  theme_bw() +
  scale_colour_manual(values=c("#000000", "#cc0000")) +
  labs(colour="Time of Day") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        #legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =12), 
        axis.title.x = element_text(size =14),
        axis.title.y = element_text(size =14), 
        strip.text.x = element_text(size =12))


## Save a plot
# ggsave("Plots/Script 4) plots/Proportion_Flight.png", 
#        width = 19, height = 22, units = "cm")










##
#### 8.4 Graze model ####
##

## Again removed the three rogue individuals like in the flight model make the model checks acceptable
## When i changes day to between dawn and dusk instead of sunrise and sunset any significant effect of shooting disappeared
## When day was classified in between sunrise and sunset the was significantly less feeding at night following shooting in the day, with no effect on day time feeding

hist(sqrt(BehSum2$Graze))
BehSum2F$DayNight <- relevel(BehSum2F$DayNight, "night")
Graze_modInt <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ shot*DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                      data = BehSum3,
                      family = betabinomial)
Graze_modFix <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ shot + DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                     data = BehSum3,
                     family = betabinomial)
Graze_modNo <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ DayNight + winter + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                     data = BehSum3,
                     family = betabinomial)

AIC(Graze_modInt, Graze_modFix, Graze_modNo)
summary(Graze_modInt)
drop1(Graze_modFix, test = "Chi")
confint(Graze_modInt)
MuMIn::r.squaredGLMM(Graze_modFix)

## Check model performance and assumptions
#model_performance(Graze_mod)
#check_model(Graze_mod)


## Model checks with DHARMa

##---- use DHARMa to get QQ plot and resids vs fitted ----##
Resids_Graze <- simulateResiduals(Graze_mod)
plot(Resids_Graze)

## plot reiduals vs predicted values
par(mfrow=c(1,3))
plotResiduals(Resids_Graze[["scaledResiduals"]], form = BehSum2$DayNight)
plotResiduals(Resids_Graze[["scaledResiduals"]], form = BehSum2$shot)
plotResiduals(Resids_Graze[["scaledResiduals"]], form = BehSum2$from_nov1)

## Check dispersion
testDispersion(Resids_Graze)
testZeroInflation(Resids_Graze)

## check autocorrelation
New_resids <- recalculateResiduals(Resids_Graze, group = BehSum2$from_nov1)
testTemporalAutocorrelation(Resids_Graze, time = unique(BehSum2$from_nov1))


## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effectsg <- predictorEffects(Graze_modInt)
plot(top_mod_effectsg)
## now extract the fits for the first variable and bind them together
effects3 <- top_mod_effectsg["shot"]
fit3 <- as.data.frame(cbind(effects3[["shot"]][["fit"]], effects3[["shot"]][["lower"]], 
                            effects3[["shot"]][["upper"]], effects3[["shot"]][["x"]][["shot"]]))
## change the names to something meaningful
setnames(fit3, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "level"))
## transform variables back to proportion scale
fit3$fit <- boot::inv.logit(fit3$fit)
fit3$lower <- boot::inv.logit(fit3$lower)
fit3$upper <- boot::inv.logit(fit3$upper)
fit3$level <- c("No Shooting", "Shooting")


## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit3, aes(x= level, ymin = lower, ymax = upper), width = 0.4, size = 0.8) +
  geom_point(data=fit3, aes(x= level, y = fit), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of bursts Grazing") +
  theme_bw() +
  #scale_colour_manual(values=c("#000000", "#cc0000")) +
  labs(colour="Time of Day") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        #legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =12), 
        axis.title.x = element_text(size =14),
        axis.title.y = element_text(size =14), 
        strip.text.x = element_text(size =12))

## Save a plot
# ggsave("Plots/Script 4) plots/Proportion_Grazing.png", 
#        width = 19, height = 22, units = "cm")


##
## Plot the effect from the best grazing model ##
##

summary(Graze_modFix)
## sums below are just the 95% CIs for the parameter
-7.630e-02 + (1.96*3.907e-02)
-7.630e-02 - (1.96*3.907e-02)
## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effectsg2 <- predictorEffects(Graze_modFix)
plot(top_mod_effectsg2)
## now extract the fits for the first variable and bind them together
effects4 <- top_mod_effectsg2["shot"]
fit4 <- as.data.frame(cbind(effects4[["shot"]][["fit"]], effects4[["shot"]][["lower"]], 
                            effects4[["shot"]][["upper"]], effects4[["shot"]][["x"]][["shot"]]))
## change the names to something meaningful
setnames(fit4, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "level"))
## transform variables back to proportion scale
fit4$fit <- boot::inv.logit(fit4$fit)
fit4$lower <- boot::inv.logit(fit4$lower)
fit4$upper <- boot::inv.logit(fit4$upper)
fit4$Shooting <- c("No Shooting", "Shooting")


## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit4, aes(x= Shooting, ymin = lower, ymax = upper), width = 0.4, size = 0.8) +
  geom_point(data=fit4, aes(x= Shooting, y = fit), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of bursts Grazing") +
  theme_bw() +
  ylim(0, 0.5) +
  scale_colour_manual(values=c("#000000", "#cc0000")) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        #legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =12), 
        axis.title.x = element_text(size =14),
        axis.title.y = element_text(size =14), 
        strip.text.x = element_text(size =12))

## Save a plot
# ggsave("Plots/Script 4) plots/Proportion_Grazing_NoInteraction.png", 
#        width = 19, height = 22, units = "cm")






##
#### 8.5 Plot all three behaviors together ####
##

## combine all of the data set
fit1$Behaviour <- "Stationary"
fit2$Behaviour <- "Flying"
fit3$Behaviour <- "Grazing"
fitAll <- rbind(fit1, fit2, fit3)

## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fitAll, aes(x= Shooting, ymin = lower, ymax = upper, colour = Time_of_Day), width = 0.4, size = 0.8) +
  geom_point(data= fitAll, aes(x= Shooting, y = fit, colour = Time_of_Day), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of accelerometer bursts") +
  facet_wrap(~Behaviour) +
  theme_bw() +
  scale_colour_manual(values=c("#000000", "#cc0000")) +
  labs(colour="Time of Day") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=15,), 
        panel.grid.minor.x = element_blank(), 
        #legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =12), 
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size =16), 
        strip.text.x = element_text(size =12))

## Save a plot
# ggsave("Plots/Script 4) plots/Proportion_Allbehaviours.png", 
#        width = 24, height = 22, units = "cm")

