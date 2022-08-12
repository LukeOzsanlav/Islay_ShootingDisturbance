## Luke Ozsanlav-Harris


##                                                      ##
#### 4. Label fixes associated with each shooting event ####
##                                                      ##


## Define function that converts timestamp immediately before shooting event to the smallest positive value
## y needs to be the timestamps from the tracking data and x is just a single value from the shooting data sequence
## Tis is needed to find the fix immediatley before the shooting event instead of just the closest which mght be forward in time
difftime2 <- function(x=x, y=y){
  
  time_difs <- difftime(y, x, units='mins')*-1
  time_difs[time_difs<0] <- 1000000000
  return(as.vector(time_difs))
}

## Define function that will loop through each tag and find the nearest fix before each shooting event while the tag was active
## The fix prior to the shooting event has to be within a certain spatial and temporal buffer to the shooting event
## TD is the tracking data. Need tag ID (Tag_ID), lat/long (Longitude/Latitude), timestamp (UTC_datetime) and unique ID for each for row (ID)
## SD is the shooting data. Need the time of the shot (timestamp), the field center (as an sf geometry) and unique ID for each for row (shot_ID)
## time_thresh is the max time in minutes that a fix can be before a shooting event
## dist_thresh is the max distance in meters the fix before a shooting event can be from the location of the shooting event
## field_crs is an sf object that has the crs of the field centers in SD

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
    
    
    ## FINDING NEAREST PRIOR FIX ##
    # find the first fix before the shooting event
    idx = sapply(shootseq, function(x) which.min(difftime2(y = Dtemp$UTC_datetime, x = x)))
    
    ## combine these fixes immediately before the shooting event with the shooting data
    aa <- (Dtemp[idx,])
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
system.time(Prox_list <- Shoot_proximity(TD = winter_gps, SD = winter_shots_cent2, time_tresh = 120, dist_thresh = 2000, field_crs = Field_centres))

## bind the lists together into one data frame
All_prox <- plyr::ldply(Prox_list)
