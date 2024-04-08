##---------------------------------------------------------##
## Luke Ozsanlav-Harris
## Created: 29/04/2022

## Aims:
## Model if ODBA changes on shooting and non-shooting days
##
## Prediction 5: 
## energy expenditure (measured using accelerometery) will be greater on days individuals 
## are exposed to shooting due to increased flight and compensatory feeding.  
##---------------------------------------------------------##


## ## Packages required
pacman::p_load(tidyverse, lubridate, data.table, zoo, svMisc, suncalc,
               geosphere, sf, sfheaders, DHARMa, MuMIn, effects, 
               performance, glmmTMB, emmeans, raster)



#--------------------------------#
#### 1. Read in the data sets ####
#--------------------------------#

## Read in all of the ODBA data for GWfG and GBG on Islay. This gives the calcualted ODBA value for each accelerometer burst and the timestamp
All_winter_ODBA <- readRDS("Biologging Data/Script4_BiologgingData.RDS")

## Read in a data sets that contains all the instances that birds were close spatially and temporally to a shooting event
## while the device was collecting accelerometer data. This can be used to wokr ot which days birds were disturbed on
All_prox <- readRDS("Biologging Data/Script4_ShootingProximity.RDS")



#-----------------------------------------------------------------------#                                                                    
#### 2. Summarize on what days and how many times birds were shot at ####
#-----------------------------------------------------------------------#                                                                  

## create a date column for the grouping
All_prox$date <- as.Date(All_prox$UTC_datetime)

## create a dummy column to sum rows
All_prox$dummy <- 1

## now summarise the  number of shooting events by each tag day
DistSum <- All_prox %>% 
            group_by(date, Tag_ID) %>% 
            summarise(n_shot = sum(dummy))




#---------------------------------------#                                                                   
#### 3. Summarize average daily ODBA ####
#---------------------------------------# 

## First label times as either day or night
## "sunrise" and "sunset" gave short days, i.e. the geese would still be out on daytime feeding sites
## Either use  "dusk", "dawn" or "nauticalDusk", "nauticalDawn" (nautical versions give the longer day length)
sun_dates <- seq.Date(as.Date(min(All_winter_ODBA$UTC_datetime)), as.Date(max(All_winter_ODBA$UTC_datetime)), by = 1)
sunt <-
  getSunlightTimes(
    date = sun_dates,
    keep = c("dusk", "dawn"),
    lat = 55.738,
    lon = 6.246,
    tz = "UTC"
  )


## Add a column for the date
All_winter_ODBA$date <- as.Date(All_winter_ODBA$UTC_datetime)

## Now join on the susnet and sunrise times
sunt$lat <- NULL
sunt$lon <- NULL
All_winter_ODBA2 <- left_join(All_winter_ODBA, sunt,  by = "date")
All_winter_ODBA2$DayNight <- ifelse(All_winter_ODBA2$UTC_datetime < All_winter_ODBA2$dawn | All_winter_ODBA2$UTC_datetime > All_winter_ODBA2$dusk,
                                    "night", "day")


## filter out ODBA values of 0
All_winter_ODBA2 <- filter(All_winter_ODBA2, !ODBA == 0)


## Create summary of behaviors before and after shooting
All_winter_ODBA2$dummy <- 1
ODBASum <- All_winter_ODBA2 %>% 
           group_by(date, device_id, DayNight, Sex, Species) %>% 
           summarise(Avg_ODBA = mean(ODBA),
                     Var_ODBA = var(ODBA),
                     Sum_ODBA = sum(ODBA),
                     n_readings = sum(dummy),
                     StEr_ODBA = sqrt(Var_ODBA / n_readings))




#-------------------------------------------------------# 
#### 4. Add on number of shooting events to the data ####
#-------------------------------------------------------# 

## join number of shots and daily behavioral summary
setnames(ODBASum, old = "device_id", new = "Tag_ID")
ODBASum2 <- left_join(ODBASum, DistSum, by = c("date", "Tag_ID"))
stopifnot(nrow(ODBASum2)==nrow(ODBASum2))


## change all the NAs in the n_shot column to 0s
ODBASum2$n_shot <- ifelse(is.na(ODBASum2$n_shot)==T, 0, ODBASum2$n_shot)




#---------------------------------------------------------------------# 
#### 5. Model Average daily ODBA on shooting vs non-shooting days ####
#---------------------------------------------------------------------#

#-----------------------------------------#
#### 5.1 Add columns for use in model ####
#-----------------------------------------#

## Add extra columns to put in to the model
ODBASum2$year <- year(ODBASum2$date) # year column
ODBASum2$year_day <- yday(ODBASum2$date) 
ODBASum2$from_nov1 <- ifelse(ODBASum2$year_day > 250, 
                            ODBASum2$year_day - 305, 
                            ODBASum2$year_day + 60) # day since November 1st
ODBASum2$from_nov1Factor <- as.factor(ODBASum2$from_nov1)
ODBASum2$shot <- ifelse(ODBASum2$n_shot >  0, "Shot", "Not_Shot") # whether there was shooting or not
ODBASum2$winter <- ifelse(ODBASum2$year_day < 150, paste0((ODBASum2$year-1), "-", ODBASum2$year), paste0(ODBASum2$year, "-", (ODBASum2$year+1)))
ODBASum2$tag_winter <- paste0(ODBASum2$Tag_ID, "_", ODBASum2$winter) # tag winter identifier
ODBASum2$prev_shot <- ifelse(lag(ODBASum2$shot, n = 2L) == "Shot", "Prev_shot", "Not")
ODBASum2$prev_shot <- ifelse(is.na(ODBASum2$prev_shot)==T, "Not", ODBASum2$prev_shot)

## set variables as correct class
ODBASum2 <- ODBASum2 %>% ungroup() %>%  mutate(across(c(year, DayNight, Tag_ID, shot, from_nov1), as.factor))
ODBASum2$from_nov1 <- as.numeric(ODBASum2$from_nov1)

## Create another column for the inverse of the varience and standard errors
ODBASum2$InVar_ODBA <- 1/ODBASum2$Var_ODBA
ODBASum2$InStEr_ODBA <- 1/ODBASum2$StEr_ODBA




#-------------------------#
#### 5.2 Prelim plots ####
#-------------------------#

## plot by species and Tag ID
ggplot(data = ODBASum2) +
  geom_boxplot(aes(Avg_ODBA, Tag_ID, colour = Species)) + 
  theme_bw()

## plot by GWfG Tag ID
ODBASum2 %>% 
  filter(Species == "GWfG") %>% 
  ggplot() +
    geom_boxplot(aes(Avg_ODBA, Tag_ID)) + 
    theme_bw()

## plot by GWfG Tag ID
ODBASum2 %>% 
  filter(Species == "GBG") %>% 
  ggplot() +
  geom_boxplot(aes(Avg_ODBA, Tag_ID)) + 
  theme_bw()

## historgam of Avg_ODBA
hist(sqrt(ODBASum2$Avg_ODBA))




#--------------------------------------------#
#### 5.3 Minimum data requirement checks ####
#--------------------------------------------#

## remove data points based on whether they had a certain number of data points per day
ODBASum3 <- filter(ODBASum2, n_readings >= 6)

## now filter out instances where an individual has 5 or less day/night worth of data
ID_sum <- table(ODBASum3$tag_winter) %>% data.frame()
ID_sum <- filter(ID_sum, Freq > 10)


## Use inner join to remove the tag winters with 5 or less days
ODBASum4 <- filter(ODBASum3, tag_winter %in% ID_sum$Var1)

## create summary to see how much data I have
datasize <- ODBASum4 %>% group_by(Species, DayNight, n_shot) %>% summarise(samp = n())
datasize

## Low samples size of shooting disturbed at night for GBG, so filtering out night for both species
ODBASum4 <- filter(ODBASum4, DayNight == "day")




#--------------------------------------------------------------#
#### 5.4 Test multiple models with different terms in them ####
#--------------------------------------------------------------#

## Start with the interaction term but then gradually drop the terms and see which has the lowest AIC

Mod_Int <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot*Species + Sex + (1|winter) + poly(from_nov1,2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_NoInt <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot + Species + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_Species <-  glmmTMB(formula = Sum_ODBA ~ n_readings + Species + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_Shot <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                      data = ODBASum4,
                      family = Gamma(link = "log"),
  
                      REML = FALSE,
                      control=glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS")))

Mod_None <-  glmmTMB(formula = Sum_ODBA ~ n_readings + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

## compare models with AIC
AIC(Mod_Int, Mod_NoInt, Mod_Species, Mod_Shot, Mod_None)

## Now get the estimates from the main mode
summary(Mod_None) # quick summary
confint(Mod_None)[1:10,] #confidence intervals
MuMIn::r.squaredGLMM(Mod_None)
confint(Mod_Int)

## model checks
ResidsMod_ReadingsWeight <- simulateResiduals(Mod_Int)
plot(ResidsMod_ReadingsWeight)
par(mfrow=c(1,3))
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$DayNight)
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$shot)
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$from_nov1)


## the below computes all of the pairwise comparisons
## and then calculates the confidence intervals for these
## create all the pairwise comparisons within the two way interaction
grid1 <- emmeans(Mod_Int, ~shot*Species)
PairComp <- as.data.frame(pairs(grid1))

## Now estimate the CIs for each of the estimates
PairComp$lowCI <- PairComp$estimate - (1.96*PairComp$SE)
PairComp$highCI <- PairComp$estimate + (1.96*PairComp$SE)




#------------------------------------------------#
#### 5.5 Plot the output from the best model ####
#------------------------------------------------#

## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effects <- predictorEffects(Mod_None)
plot(top_mod_effects)

## now extract the fits for the first variable and bind them together
effects3 <- top_mod_effects["shot"]
fit3 <- as.data.frame(cbind(as.numeric(effects3[["shot"]][["fit"]]), as.numeric(effects3[["shot"]][["lower"]]), 
                            as.numeric(effects3[["shot"]][["upper"]]), as.character(effects3[["shot"]][["x"]][["shot"]]),
                            as.character(effects3[["shot"]][["x"]][["Species"]])))
## change the names to something meaningful
setnames(fit3, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "shot", "Species"))
## transform variables back to proportion scale
fit3$fit <- exp(as.numeric(as.character(fit3$fit)))
fit3$lower <- exp(as.numeric(as.character(fit3$lower)))
fit3$upper <- exp(as.numeric(as.character(fit3$upper)))

## Change the names of some of the variables for plotting
fit3$shot <- ifelse(fit3$shot == "Shot", "Shooting", "No Shooting")
ODBASum5 <- ODBASum4
ODBASum5$shot <- ifelse(ODBASum5$shot == "Shot", "Shooting", "No Shooting")

## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit3, aes(x= shot, ymin = lower, ymax = upper, colour = Species), width = 0.4, size = 0.8) +
  geom_point(data=fit3, aes(x= shot, y = fit, colour = Species), size = 2, shape = 15)  +
  #geom_point(data = ODBASum5, aes(x =shot, y= Sum_ODBA, colour = Species), size = 0.75, alpha = 0.1) +
  xlab("Shooting Exposure") + ylab("Average daily ODBA") +
  theme_bw() +
  facet_wrap(~Species)+
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  labs(colour="Species") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =14), 
        axis.title.x = element_text(size =18),
        axis.title.y = element_text(size =18), 
        strip.text.x = element_text(size =14), 
        legend.title = element_text(size =14),
        legend.text = element_text(size =14))

## Save a plot
# ggsave("Plots/Script 4) plots/Average daily ODBA by species.png", 
#        width = 18, height = 20, units = "cm")

