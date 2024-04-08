##----------------------------------------##
## Luke Ozsanlav-Harris
## Created: 13/01/2022

## Aim: 
## Examining the effect that shooting has on movement
## In particular comparing the step length immediately before and during shooting disturbance
##
## Prediction 1: 
## shooting disturbance will cause flight initiation 
## but this effect will decay with increasing distance from the shooting event. 
##
##----------------------------------------##

## Packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(lme4)
library(DHARMa)
library(effects)
library(performance)
library(lmerTest)
library(ggplot2)



##                                                               ##
#### 0.  Read in biologging data ####
##                                                               ##

## Read in the data set for modelling
## This contains the GPS tracking data for each individual and identifies any fixes where the bird
## was disturbed by shooting before the next fix was taken, along with some data about the shooting event. 
## Lat/longs have been rounded but the step lengths were already calculated using more precise location. 
Biolog <- readRDS("Biologging Data/Script2_BiologgingData.RDS")



##                                                               ##
#### 1. Identify suitable paired observations within individuals ####
##                                                               ##

##
#### 1.1 Shoot and t-1 as control ####
##

## Compare the disturbed fix to both the fix immediately before and after
## Copy the main data set just to find suitable pairs for now
ppairs <- Biolog

## Add columns to the dataset for the time difs, shoot status and ID of the previous row....
ppairs$TD_control <- lag(ppairs$time_dif)
ppairs$shoot_control <- ifelse(is.na(lag(ppairs$Shoot_ID1))==T, 0, 1)
ppairs$ID_control <- lag(ppairs$ID)
ppairs$SL_control <- lag(ppairs$Step_len)

## Now just filter out the rows that experienced a shooting event, to streamline data set a bit now
ppairs2 <- ppairs %>% filter(is.na(Shoot_ID1) == F) 
## filter out instances were the disturbance step length is calculated from a fix interval greater than 130 mins
ppairs2 <- ppairs2 %>% filter(time_dif < 130)
## filter out instances when the fix interval of the proposed pair is greater than 10 minutes different
ppairs2 <- ppairs2 %>% filter(TD_control < time_dif + 10 & TD_control > time_dif - 10)
## filter out instances where the paired observation is also a disturbed
ppairs2 <- ppairs2 %>% filter(shoot_control == 0)


## Now extract the pairs and give each pair and their time dif and step length
Dist_fix <- subset(ppairs2, select = c("ID", "time_dif", "Step_len"))
Pair_fix <- subset(ppairs2, select = c("ID_control", "TD_control", "SL_control"))

## label each pair
Dist_fix$Pair_ID <- 1:nrow(ppairs2)

## calculate the speeds for both data sets, swapped these around as the 
Dist_fix$speed <- (Dist_fix$Step_len*1000)/(Dist_fix$time_dif*60) #meters/seconds
Pair_fix$prev_speed <- (Pair_fix$SL_control*1000)/(Pair_fix$TD_control*60)

## bind the data sets side by side now
final_pairs <- cbind(Dist_fix, Pair_fix)

## signify whether it is before of after the disturbance fix
final_pairs$control <- "Shoot_T-1"



##
#### 1.2 Bind all pairs with shooting info ####
##

# ## Bind the before and after paris together
# final_pairs3 <- plyr::rbind.fill(final_pairs, final_pairs2, final_pairsX)
# 
# ## for t-1/t+1 we had to already assign the ID if the shooting event fix as it isnt either of the fixes in the sample
# ## for the other tow samples the shooting event ID is just the ID
# final_pairs3$ID_Shot <- ifelse(is.na(final_pairs3$ID_Shot)==T, final_pairs3$ID, final_pairs3$ID_Shot)

## remove the columns that are already present in final_pairs
Biolog <- subset(Biolog, select = c(-Step_len, -time_dif))

## bind to final_pairs data set
final_pairs2 <- left_join(final_pairs, Biolog, by = "ID")
stopifnot(nrow(final_pairs2) == nrow(final_pairs))

## Just going to run this analysis with pairs were the fix length is roughly one hour
final_pairs2 <- filter(final_pairs2, time_dif > 50 & time_dif < 70)





##                                                                   
#### 2. Prepare data for models ####
##                                                                    

## Some disturbance fixes are exposed to multiple shooting events
## Therefore going to pick the minimum distance for each for
## First pick out the column with the distance to shot data
Shot_dists <- subset(final_pairs2, select =c("Shoot_ID1_dist_to", "Shoot_ID2_dist_to", "Shoot_ID3_dist_to", "Shoot_ID4_dist_to"))
final_pairs2$min_shot_dist <- apply(Shot_dists, 1, FUN = min, na.rm = TRUE)
summary(final_pairs2$min_shot_dist)

## If I modeled these as differences in speeds then this would standardize for fix length
final_pairs2$speed_dif <- final_pairs2$speed - final_pairs2$prev_speed
final_pairs2$Dist_difm <- (final_pairs2$Step_len*1000) - (final_pairs2$SL_control*1000)


## Add other variables to the data set that will go in the model
## year and column
final_pairs2$year <- year(final_pairs2$UTC_datetime)
## shooting time of day column
final_pairs2$TOD_ <- as.numeric(hour(final_pairs2$Shoot_ID1_time)) + as.numeric(minute(final_pairs2$Shoot_ID1_time)/60)
## days Since November 1st column
final_pairs2$year_day <- yday(final_pairs2$UTC_datetime)
final_pairs2$from_nov1 <- ifelse(final_pairs2$year_day > 250, 
                           final_pairs2$year_day - 305, 
                           final_pairs2$year_day + 60)
## add winter column
final_pairs2$year <- as.numeric(as.character(final_pairs2$year))
final_pairs2$winter <- ifelse(final_pairs2$year_day < 150, paste0((final_pairs2$year-1), "-", final_pairs2$year), paste0(final_pairs2$year, "-", (final_pairs2$year+1)))


## scale continuous variables
final_pairs2$TOD_sc <- scale(final_pairs2$TOD_)
final_pairs2$from_nov1_sc <- scale(final_pairs2$from_nov1)
final_pairs2$min_shot_dist_sc <- scale(final_pairs2$min_shot_dist)

## set variables as the right classes
final_pairs2$Tag_ID <- as.factor(final_pairs2$Tag_ID)
final_pairs2$year <- as.factor(final_pairs2$year)
final_pairs2$species <- as.factor(final_pairs2$species)

## check the distribution of the response
hist(final_pairs2$speed_dif)
hist(final_pairs2$Dist_difm)

  


##
##### 3. Model Step Length Change as Function of Distance to Shot ####
##

## 
#### 3.1 Prepare Data sets ####
##

## Pick which data sets I want to put into the models
Consec_pair <- dplyr::filter(final_pairs2, control == "Shoot_T-1")
Consec_pair <- dplyr::filter(Consec_pair, is.na(species) == F) # removes any rows were the data for the shooting event is missing

## filter the data set for each species
GBG_pairs <- dplyr::filter(Consec_pair, species == "GBG")
GWfG_pairs <- dplyr::filter(Consec_pair, species == "GWfG")

## remove one big value, in this instance the bird actually started migration so capturing a different behaviour to this analysis
GWfG_pairs <- GWfG_pairs %>%  filter(Dist_difm > -10000)



## 
#### 3.2 GBG Model ####
##

## Model for Barnacle Geese
model_GB <- lmer(Dist_difm ~ log(min_shot_dist) + poly(TOD_,2) + winter + (1|Tag_ID), 
                    data= GBG_pairs,
                    REML=FALSE)

## Investigate model
summary(model_GB) # model summary
confint(model_GB) # model CIs
performance::check_model(model_GB) # model check


## Use the "effects" package to create a data set that can be used to plot the effects of distance to shot
## Extract estimates across full range of distances to shot
GBG_mod_effects <- predictorEffects(model_GB, focal.levels = 4000)
plot(GBG_mod_effects[1]) ## check the plot

## Extract just the distance to shot estimates
effectsGBG <- GBG_mod_effects[1]

## Manipulate data frame for plotting
fitGBG <- as.data.frame(cbind(effectsGBG[["min_shot_dist"]][["fit"]], effectsGBG[["min_shot_dist"]][["lower"]], 
                              effectsGBG[["min_shot_dist"]][["upper"]], effectsGBG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGBG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
fitGBG$Species <- "Barnacle" # add on species column



## 
#### 3.3 GWfG Model ####
##

## Model for GWfG
## Winter now a random intercept term as there are many more levels in this variable
model_GW <- lmer(Dist_difm ~ log(min_shot_dist) + poly(TOD_,2) + (1|winter) + (1|Tag_ID),
                    REML=FALSE,
                    data= GWfG_pairs)

## Investigate model
summary(model_GW) # model summary
confint(model_GW) # model CIs
performance::check_model(model_GW) # model check


## Use the "effects" package to create a data set that can be used to plot the effects of distance to shot
## Extract estimates across full range of distances to shot
GWfG_mod_effects <- predictorEffects(model_GW, focal.levels = 4000)
plot(GWfG_mod_effects[1]) ## check the plot

## Extract just the distance to shot estimates
effectsGWfG <- GWfG_mod_effects[1]

## Manipulate data frame for plotting
fitGWfG <- as.data.frame(cbind(effectsGWfG[["min_shot_dist"]][["fit"]], effectsGWfG[["min_shot_dist"]][["lower"]], 
                               effectsGWfG[["min_shot_dist"]][["upper"]], effectsGWfG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGWfG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
fitGWfG$Species <- "White-fronted" # add on species column



## 
#### 3.4 Identify point regression lines cross y=0  ####
##

## **Note using this data set it is possible to find the points that the 95% CI of the regression line
##   first cross the line y=0
CrossPointGWfG <- filter(fitGWfG, lower <0 & lag(lower)>0)
CrossPointGBG <- filter(fitGBG, lower <0 & lag(lower)>0)



## 
#### 3.5 Plot GBG and GWfG model outputs ####
##

## Join together the model estimates for GBG and GWfG
Bestfits <- rbind(fitGBG, fitGWfG) 

## Add this species column to keep naming consisten with the MS
Bestfits$Species <- ifelse(Bestfits$Species == "White-fronted", "GWfG", "GBG")

## Make the plot
ggplot() + 
  geom_ribbon(data = Bestfits, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Species, colour = Species), 
              alpha = 0.3, colour = NA, fill = "grey") +
  geom_line(data=Bestfits, aes(x= dist_to_shot, y = fit, group = Species, colour = Species), size = 1.25)  +
  ylab("Difference in distance (m) [t-(t-1)]") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title=element_text(size=12,), 
        axis.text = element_text(size =14), 
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size =16), 
        strip.text.x = element_text(size =14), 
        legend.title = element_text(size =14),
        legend.text = element_text(size =14))

## save the plot 
ggsave("Plots/Script 2) plots/Displacement [t-(t-1)] Distance version.png", 
       width = 25, height = 15, units = "cm")


