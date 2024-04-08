##---------------------------------------------------##
## Luke Ozsanlav-Harris
## Created: 21/04/22

## Aim: 
## Compare behavior on days where individuals are disturbed by shooting vs
## not disturbed by shooting during the day. 
##
## Prediction 4: 
## shooting disturbance will increase flight as individuals are displaced, 
## increase compensatory feeding and increase vigilance in response to other stimuli post shooting. 

##---------------------------------------------------##

## Packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(suncalc)
library(DHARMa)
library(effects)
library(performance)
library(glmmTMB)




##                                                                     
#### 1. Read in the data sets ####
##

## Read in the classified accelerometer data set from Islay, Scotland
Islay_orn_behav <- read_rds("Biologging Data/Script3_BiologgingData.RDS")

## Read in a data sets that contains all the instances that birds were close spatially and temporally to a shooting event
## while the device was collecting accelerometer data. This can be used to wokr ot which days birds were disturbed on
All_prox <- read_rds("Biologging Data/Script3_ShootingProximity.RDS") %>% select(-geometry)




##                                                                     
#### 2. Summarize on what days and how many times birds were shot at ####
##                                                                     

## create a date column for the grouping
All_prox$date <- as.Date(All_prox$UTC_datetime)

## create a dummy column to sum rows
All_prox$dummy <- 1

## now summarise the  number of shooting events by each tag day
DistSum <- All_prox %>% 
           as.data.frame() %>% 
           group_by(date, Tag_ID) %>% 
           summarise(n_shot = sum(dummy))




##
#### 3. Summarize daily behavior ####
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
#### 4. Add on number of shooing events to the dataset ####
##

## join number of shots and daily behavioral summary
BehSum$Tag_ID <- as.character(BehSum$Tag_ID)
BehSum2 <- left_join(BehSum, DistSum, by = c("date", "Tag_ID"))
stopifnot(nrow(BehSum2)==nrow(BehSum2))

## change all the NAs in the n_shot column to 0s
BehSum2$n_shot <- ifelse(is.na(BehSum2$n_shot)==T, 0, BehSum2$n_shot)




##
#### 5. Model the proportion of time spent of certain behaviors varies ####
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




##
#### 6. Minimum data requirement checks ####
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
#### 6.1 Stationary model ####
##

## Compare models to determine effect of shooting variable
Stat_modInt <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ shot*DayNight + (1|winter) + poly(from_nov1,2) + 
                                                         (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                        data = BehSum3,
                        family = betabinomial)
Stat_modfix <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ shot + DayNight + (1|winter) + poly(from_nov1,2) + 
                                                         (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                       data = BehSum3,
                       family = betabinomial)
Stat_modNo <- glmmTMB(formula = cbind(Stat, Not_Stat) ~ DayNight + (1|winter) + poly(from_nov1,2) + 
                                                        (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                       data = BehSum3,
                       family = betabinomial)

## compare models
AIC(Stat_modInt, Stat_modfix, Stat_modNo)

## return model parameters and CIs
summary(Stat_modNo)
emmeans::emmeans(Stat_mod, ~shot*DayNight, type = "response")
confint(Stat_modInt)[1:10,]


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
# ggsave("Plots/Script 3) plots/Proportion_Stationary.png", 
#        width = 19, height = 22, units = "cm")




##
#### 6.2 Flight model ####
##

## compare three models
hist(sqrt(BehSum2$Flight))
BehSum3$DayNight <- relevel(BehSum3$DayNight, "day")

Flight_modInt <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ shot*DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                    data = BehSum3,
                    family = betabinomial)
Flight_modFix <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ shot + DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                         data = BehSum3,
                         family = betabinomial)
Flight_modNo <- glmmTMB(formula = cbind(Flight, Not_Flight) ~ DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                         data = BehSum3,
                         family = betabinomial)

## Compare models
AIC(Flight_modInt, Flight_modFix, Flight_modNo)

## return model parameters and CIs
summary(Flight_modInt)
drop1(Flight_modNo, test = "Chi")
confint(Flight_modFix)[1:10,]


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
top_mod_effectsF <- predictorEffects(Flight_modInt)
plot(top_mod_effectsF)
## now extract the fits for the first variable and bind them together
effects2 <- top_mod_effectsF["shot"]
fit2 <- as.data.frame(cbind(effects2[["shot"]][["fit"]], effects2[["shot"]][["lower"]], 
                            effects2[["shot"]][["upper"]], effects2[["shot"]][["x"]][["shot"]],
                            effects2[["shot"]][["x"]][["DayNight"]]))
## change the names to something meaningful
setnames(fit2, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "shotlevel", "Daylevel"))
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
# ggsave("Plots/Script 3) plots/Proportion_Flight.png", 
#        width = 19, height = 22, units = "cm")




##
#### 6.3 Graze model ####
##

## Again removed the three rogue individuals like in the flight model make the model checks acceptable
## When i changes day to between dawn and dusk instead of sunrise and sunset any significant effect of shooting disappeared
## When day was classified in between sunrise and sunset the was significantly less feeding at night following shooting in the day, with no effect on day time feeding

hist(sqrt(BehSum2$Graze))

Graze_modInt <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ shot*DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                      data = BehSum3,
                      family = betabinomial)
Graze_modFix <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ shot + DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                     data = BehSum3,
                     family = betabinomial)
Graze_modNo <- glmmTMB(formula = cbind(Graze, Not_Graze) ~ DayNight + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                     data = BehSum3,
                     family = betabinomial)

## Compare models
AIC(Graze_modInt, Graze_modFix, Graze_modNo)

## return model parameters and CIs
summary(Graze_modFix)
confint(Graze_modInt)[1:10,]
MuMIn::r.squaredGLMM(Graze_modFix)



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
top_mod_effectsg <- predictorEffects(Graze_modFix)
plot(top_mod_effectsg)
## now extract the fits for the first variable and bind them together
effects3 <- top_mod_effectsg["shot"]
fit3 <- as.data.frame(cbind(effects3[["shot"]][["fit"]], effects3[["shot"]][["lower"]], 
                            effects3[["shot"]][["upper"]], effects3[["shot"]][["x"]][["shot"]]))
## change the names to something meaningful
setnames(fit3, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "shotlevel"))
## transform variables back to proportion scale
fit3$fit <- boot::inv.logit(fit3$fit)
fit3$lower <- boot::inv.logit(fit3$lower)
fit3$upper <- boot::inv.logit(fit3$upper)
fit3$shotlevel <- c("No Shooting", "Shooting")



## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit3, aes(x= shotlevel, ymin = lower, ymax = upper), width = 0.4, size = 0.8) +
  geom_point(data=fit3, aes(x= shotlevel, y = fit), size = 1.5)  +
  xlab("Shooting Exposure") + ylab("Proportion of bursts Grazing") +
  theme_bw() +
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
# ggsave("Plots/Script 3) plots/Proportion_Grazing.png", 
#        width = 19, height = 22, units = "cm")

