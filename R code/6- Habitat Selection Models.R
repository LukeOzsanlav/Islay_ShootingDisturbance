##---------------------------------------------------------##
## Aimée McIntosh
## Created 12/09/23

## Aims:
## Model whether foraging habitat selection differs between and within shooting disturbance days
## 
## Prediction 6: 
## Individuals disturbed by shooting should show weaker selection for preferred improved grassland fields 
## and/or select sites further from shooting and other forms of anthropogenic disturbance (i.e. further from roads)
##---------------------------------------------------------##

## Packages required
pacman::p_load(tidyr,plyr,dplyr,zoo,data.table,move,ggplot2,patchwork, purrr,readr,sfheaders,
               lubridate,tidyverse,readr,amt,raster, geosphere,sf,DHARMa,lme4,pROC,
               MuMIn,effects,lattice,performance,MASS,glmmTMB,bbmle,emmeans, MetBrewer)



#--------------------------------#
####  1. Read and prep data   ####
#--------------------------------#
# This includes the full data used in final models. GPS fixes and specific farm locations have been removed/anonamized for GDPR. 

## 1.1 Between disturbance days data ----
gbg_dd_dat <- read.csv("Derived Data/GBG_1200m_Between_Model_Data.csv") # GBG data
gwf_dd_dat <- read.csv("Derived Data/GWFG_644m_Between_Model_Data.csv") # GWfG data

gbg_dd_dat$Species <- "GBG"
gwf_dd_dat$Species <- "GWfG"

gbg_dd_dat$Winter <- as.factor(gbg_dd_dat$Winter)
gbg_dd_dat$Disturbance_Day <- as.factor(gbg_dd_dat$Disturbance_Day)

gwf_dd_dat$Winter <- as.factor(gwf_dd_dat$Winter)
gwf_dd_dat$Disturbance_Day <- as.factor(gwf_dd_dat$Disturbance_Day)

## 1.2 Within disturbance days data ----
gbg_PP <- read.csv("Derived Data/GBG_1200m_Within_Model_Data.csv")
gwf_PP <- read.csv("Derived Data/GWFG_644m_Within_Model_Data.csv")



#-------------------------------------------#
####  2. Between disturbance GBG models  ####
#-------------------------------------------#

## 2.1 Set weights and re-level model ----
gbg_dd_dat$weight2 <- ifelse(gbg_dd_dat$case_ == "TRUE", 10, 1)
gbg_dd_dat$ceh_broad <- relevel(as.factor(gbg_dd_dat$ceh_broad),"other") # set 'other' habitat as reference level in the model

## 2.2 Run models for habitat selection ----
# Model with disturbance-habitat interaction 
gbg_dd_rsfc1 <- glm(Case ~ ceh_broad + Disturbance_Day + ceh_broad*Disturbance_Day, family=binomial(), data = gbg_dd_dat, weights = weight2)
summary(gbg_dd_rsfc1)

# Model without interaction 
gbg_dd_rsfc2 <- glm(Case ~ ceh_broad + Disturbance_Day, family=binomial(), data = gbg_dd_dat, weights = weight2)
summary(gbg_dd_rsfc2)

## 2.3 Model checks ----
# AIC comparison table
ICtab(gbg_dd_rsfc1, gbg_dd_rsfc2,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

# Calculate auc for each model
p <- as.numeric(predict(gbg_dd_rsfc1, type="response")> 0.5) # Set for best performing model, change model to compare 

auc(roc(gbg_dd_dat$Case ~ p))

# Confusion matrix for correct classification etc
# Compare model fit to raw data
gbg_dd_dat$Case.factor <- as.factor(gbg_dd_dat$Case)

p2 <- as.factor(p)

confusionMatrix(p2, gbg_dd_dat$Case.factor, positive = "1")

# Create a data frame of model variables for the interaction terms in the model 
gbgSHOOT.HISTa <- as.data.frame(effect("ceh_broad:Disturbance_Day",gbg_dd_rsfc1)) 
head(gbgSHOOT.HISTa)



#--------------------------------------------#
####  3. Between disturbance GWfG models  ####
#--------------------------------------------#

## 3.1 Set weights and re-level model ----
gwf_dd_dat$weight2 <- ifelse(gwf_dd_dat$case_ == "TRUE", 10, 1)
gwf_dd_dat$ceh_broad <- relevel(as.factor(gwf_dd_dat$ceh_broad),"other") # set 'other' habitat as reference level in the model 

## 3.2 Run models for habitat selection ----

gwf_dd_rsfc1 <- glm(Case ~ ceh_broad + Disturbance_Day + ceh_broad*Disturbance_Day, family=binomial(), data = gwf_dd_dat, weights = weight2)
summary(gwf_dd_rsfc1)

gwf_dd_rsfc2 <- glm(Case ~ ceh_broad + Disturbance_Day, family=binomial(), data = gwf_dd_dat, weights = weight2)

## 3.3 Model checks ----
# AIC comparison table
ICtab(gwf_dd_rsfc1, gwf_dd_rsfc2,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

# Calculate AUC for each model 
p <- as.numeric(predict(gwf_dd_rsfc1, type="response")> 0.5) # Set for best performing model, change model to compare

auc(roc(gwf_dd_dat$Case ~ p))

# Confusion matrix for correct classification etc
# Compare model fit to raw data
gwf_dd_dat$Case.factor <- as.factor(gwf_dd_dat$Case)

p2 <- as.factor(p)

confusionMatrix(p2, gwf_dd_dat$Case.factor, positive = "1")

# Create a data frame of model variables for the interaction terms in the model 
gwfSHOOT.HISTa <- as.data.frame(effect("ceh_broad:Disturbance_Day",gwf_dd_rsfc1)) #100 to get nice line



#--------------------------------------------#
####          4. Create Figure 5          ####
#--------------------------------------------#

## 4.1 Create data frame with model variables for both GBG and GWfG ----

gbgSHOOT.HISTa$Species <- "GBG"
gbgSHOOT.HISTa <- dplyr::filter(gbgSHOOT.HISTa, ceh_broad != c("other"))
gbgSHOOT.HISTa$ceh_broad <- recode(gbgSHOOT.HISTa$ceh_broad, "coastal_saltmarsh" = "saltmarsh_coastal")

gwfSHOOT.HISTa$Species <- "GWfG"
goose_2way <- rbind(gbgSHOOT.HISTa, gwfSHOOT.HISTa)

## 4.2 Create Figure 5a ----
# Re-code Habitat variables
goose_2way$ceh_broad <- dplyr::recode(goose_2way$ceh_broad, arable_hort = "Other Arable", bog = "Bog",
                                      improved_grassland = "Improved grassland",other_grassland = "Other grassland",
                                      saltmarsh_coastal = "Saltmarsh & coastal", freshwater = "Freshwater",
                                      other = "Other")

# Re-code Shooting fix
goose_2way$Disturbance_Day <- recode(goose_2way$Disturbance_Day, "1" = "Shooting", "0" = "No Shooting" )

# goose_3way$PP_Shoot_f <- factor(goose_3way$PP_Shoot, levels = c("Pre-shooting", "Post-shooting"))

pd <- position_dodge(0.5) # Dodge point position to avoid complete overlap in the figure

Fig_5a <- ggplot(goose_2way, aes(ceh_broad, fit))+
  geom_hline(yintercept = 0.5, colour ="black", linetype = "dashed") +
  geom_point(aes(color = Disturbance_Day), position = pd, size = 3)+
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Disturbance_Day), width=1.5, position = pd) +
  scale_colour_manual(values = c("black", "red")) +
  labs(x = "Habitat", y = "Relative probability of habitat use", colour = "Disturbance")+
  facet_grid(~Species) +
  #ggtitle("a)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", size = 0.7),
        axis.text.x = element_text(size=11, angle  = 30, hjust = 1), 
        axis.title=element_text(size=15),
        text = element_text(size = 15)) +
  #  scale_y_continuous(limits = c(0,6000))
  scale_y_continuous(limits = c(0,1))

Fig_5a

# ggsave("Goose_Total_Habitat_Selection_1200m.jpeg",Fig_5a) # Save plot


## 4.3 Prep data to create Figure 5b ----

# Create plot of proportion of fixes in each habitat category by day and species using real GPS fixes
all_goose_dd <- rbind(gbg_dd_dat, gwf_dd_dat)

# Create a column of 1s to aggregate by
all_goose_dd$sum <- 1

# Create a Breakdown for all birds by site and case for real points and pseudoabsences
# Create a data frame of fixes in each habitat type by species by shooting day
goosefixdf <- with(all_goose_dd, aggregate(sum, by = list(Species, Disturbance_Day, Used,
                                                          ceh_broad), "sum"))

names(goosefixdf)[1] <- "Species"
names(goosefixdf)[2] <- "Disturbance_Day"
names(goosefixdf)[3] <- "Used"
names(goosefixdf)[4] <- "ceh_broad"
names(goosefixdf)[5] <- "number_of_fixes"

# Create a data frame of total fixes per year
goosefixtot <- with(all_goose_dd, aggregate(sum, by = list(Species, Disturbance_Day, Used), "sum"))

names(goosefixtot)[1] <- "Species"
names(goosefixtot)[2] <- "Disturbance_Day"
names(goosefixtot)[3] <- "Used"
names(goosefixtot)[4] <- "total_fixes"

# Bind data frames together using a full_join
alldf_all <- full_join(goosefixdf, goosefixtot, by = c("Species", "Disturbance_Day", "Used"))
head(alldf_all)

# Order data frame by site (Species)
alldf_all <- alldf_all[order(alldf_all$Species), ]
head(alldf_all)

# Create a column for the proportion of fixes 
alldf_all$ppn_fixes <- alldf_all$number_of_fixes/alldf_all$total_fixes

# Convert this to a percentage and create as a new column
alldf_all$percentage_fixes <- alldf_all$ppn_fixes*100

# Remove temporary data frames
rm(allfixdf, allfixtot)

## 4.4 Create Figure 5b ----

# Re-code habitat classifications for figure
alldf_all$ceh_broad <- dplyr::recode(alldf_all$ceh_broad, arable_hort = "Arable", bog = "Bog", improved_grassland = "Improved grassland", other_grassland = "Other grassland", saltmarsh_coastal = "Saltmarsh & coastal", freshwater = "Freshwater", other = "Other")

# Re-code Shooting day
alldf_all$Disturbance_Day <- recode(alldf_all$Disturbance_Day, "0" = "No Shooting", "1" = "Shooting" )

# Re-code GPS fix
alldf_all$Used <- recode(alldf_all$Used, "Avail" = "Pseudo", "Used" = "Used" )

# Create colour palette for the figure
x <- met.brewer("Lakota", n = 9) # using Metbrewer colour palettes
x[1];x[2];x[3];x[4];x[5];x[6];x[7];x[8];x[9] # Get a list of the codes for each colour 
x # View the colours

# Create plot 
Fig_5b <- ggplot(alldf_all, aes(x = Used, y = percentage_fixes)) +
  geom_col(aes(fill = ceh_broad)) +
  #scale_fill_manual()) +
  facet_grid(Species ~ Disturbance_Day) +
  labs(x = "Point Type", y = "Proportion of locations (%)", fill = "Habitat Type") +
  theme_bw() +
  scale_fill_manual(values = c("#DABB48","#CA8740","#346565", "#20235B", "#6FAF82", "#B85D32", "#04A3BD")) +
  ggtitle("b)") +
  theme(panel.grid.major = element_blank(), 
        ##panel.grid.minor = element_blank(), 
        ##panel.border = element_blank(), 
        ##axis.line = element_line(colour = "black", size=0.7),
        ##axis.text = element_text(size=12, angle = 90, hjust = 1), 
        axis.title=element_text(size=15),
        text = element_text(size = 15),
        strip.text.x = element_text(size = 15))
Fig_5b



# 4.5 Combine plots to create completed Figure 5 ----
Figure_5 <- Fig_5a / Fig_5b
Figure_5
ggsave(filename = "Plots/Script 6) plots/Figure_5.png", plot = Figure_5, dpi = 320, units = "in", height = 10.7, width = 8.59) # Save Figure 5




#-----------------------------------------------#
####  5. Within disturbance days GBG models  ####
#-----------------------------------------------#

## 5.1 Set weights, prep data and re-level factors for models ----
gbg_PP <- mutate(gbg_PP, ceh_simp = ifelse(ceh_broad != "improved_grassland", "other", "improved_grassland"))
gbg_PP <- dplyr::mutate(gbg_PP, dist_sc = scale(gbg_dist_points))
gbg_PP$ceh_simp <- relevel(as.factor(gbg_PP$ceh_simp),"other")
gbg_PP$ceh_simp_group <- gbg_PP$ceh_broad # create column for the broader shooting pre- vs post disturbance

# recode ceh_broad into categories according to whether there are 50 + fixes in each pre and post category for each habitat
gbg_PP$ceh_simp_group <- recode(gbg_PP$ceh_simp_group, "freshwater" = "other", "bog" = "other", "other_grassland" = "other" )

gbg_PP$ceh_simp_group <- relevel(as.factor(gbg_PP$ceh_simp_group),"other")

## 5.2 Run models for habitat selection ----
gbg_pp_rsfc1 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc + ceh_simp_group*dist_sc +
                      ceh_simp_group*PP_Shoot*dist_sc, family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc2 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc + PP_Shoot*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc3 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc4 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + PP_Shoot*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc5 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*dist_sc + PP_Shoot*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc6 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc7 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc8 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      PP_Shoot*dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

gbg_pp_rsfc9 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc,
                    family=binomial(), data = gbg_PP, weights = weight2)

## 5.3 Model checks ----
# AIC comparison tables 
ICtab(gbg_pp_rsfc1, gbg_pp_rsfc2,gbg_pp_rsfc3,gbg_pp_rsfc4,gbg_pp_rsfc5,gbg_pp_rsfc6,gbg_pp_rsfc7,gbg_pp_rsfc8,gbg_pp_rsfc9,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

# Calculate AUC for each model 
p <- as.numeric(predict(gbg_pp_rsfc1, type="response")> 0.5) # Set for best performing model, change model to compare

auc(roc(gbg_PP$Case ~ p))

# Confusion matrix for correct classification etc
# Compare model fit to raw data
gbg_PP$Case.factor <- as.factor(gbg_PP$Case)

p2 <- as.factor(p)

confusionMatrix(p2, gbg_PP$Case.factor, positive = "1")

# Create data frame of model estimates for interaction terms - best performing model included 3-way interaction term 
gbgSHOOT.HISTa <- as.data.frame(effect("ceh_simp_group:PP_Shoot:dist_sc",gbg_pp_rsfc1, xlevels = 100)) #100 to get nice line

# Transform scaled and centered distance metrics
gbgSHOOT.HISTa$Road_Distance <- (gbgSHOOT.HISTa$dist_sc*sd(gbg_PP$gbg_dist_points)) + mean(gbg_PP$gbg_dist_points)

# drop freshwater and other from plot as there are no real points for these categories
#gbgSHOOT.HISTa <- dplyr::filter(gbgSHOOT.HISTa, ceh_broad != c("freshwater", "other"))


#------------------------------------------------#
####  6. Within disturbance days GWfG models  ####
#------------------------------------------------#

## 6.1 Set weights, prep data and re-level factors for models ----
gwf_PP$ceh_simp_group <- gwf_PP$ceh_broad

# recode ceh_broad into categories according to whether there are 50 + fixes in each pre and post category for each habitat
gwf_PP$ceh_simp_group <- recode(gwf_PP$ceh_simp_group, "freshwater" = "other", "bog" = "other", "saltmarsh_coastal" = "other" )

gwf_PP$dist_sc <- scale(gwf_PP$gwf_dist_points)
gwf_PP$Case <- as.factor(gwf_PP$Case)
gwf_PP$Winter <- as.factor(gwf_PP$Winter)

gwf_PP$weight2 <- ifelse(gwf_PP$case_ == "TRUE", 10, 1)

gwf_PP$ceh_simp_group <- relevel(as.factor(gwf_PP$ceh_simp_group),"other")

## 6.2 Run models for habitat selection ----

gwf_pp_rsfc1 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc + PP_Shoot*dist_sc +
                      ceh_simp_group*PP_Shoot*dist_sc, family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc2 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc + PP_Shoot*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc3 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + ceh_simp_group*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc4 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot + PP_Shoot*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc5 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*dist_sc + PP_Shoot*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc6 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*PP_Shoot,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc7 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      ceh_simp_group*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc8 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc +
                      PP_Shoot*dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

gwf_pp_rsfc9 <- glm(Case ~ ceh_simp_group + PP_Shoot + dist_sc,
                    family=binomial(), data = gwf_PP, weights = weight2)

## 6.3 Model checks ----
# AIC comparison tables 
ICtab(gwf_pp_rsfc1, gwf_pp_rsfc2,gwf_pp_rsfc3,gwf_pp_rsfc4,gwf_pp_rsfc5,gwf_pp_rsfc6,gwf_pp_rsfc7,gwf_pp_rsfc8,gwf_pp_rsfc9,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

# Calculate AUC for each model 
p <- as.numeric(predict(gwf_pp_rsfc1, type="response")> 0.5) # Set for best performing model, change model to compare

auc(roc(gwf_PP$Case ~ p))

# Confusion matrix for correct classification etc
# Compare model fit to raw data
gwf_PP$Case.factor <- as.factor(gwf_PP$Case)

p2 <- as.factor(p)

confusionMatrix(p2, gwf_PP$Case.factor, positive = "1")

# Create data frame of model estimates for interaction terms - best performing model included 3-way interaction term 
gwfSHOOT.HISTa <- as.data.frame(effect("ceh_simp_group:PP_Shoot:dist_sc",gwf_pp_rsfc1, xlevels = 100)) #100 to get nice line

# Transform scaled and centred distance metrics
gwfSHOOT.HISTa$Road_Distance <- (gwfSHOOT.HISTa$dist_sc*sd(gwf_PP$gwf_dist_points)) + mean(gwf_PP$gwf_dist_points)

## 6.4 Create Figure 6a ----
# Combine model estiamte data frames for each species
gbgSHOOT.HISTa$Species <- "GBG"
gwfSHOOT.HISTa$Species <- "GWfG"
goose_3way <- rbind(gbgSHOOT.HISTa, gwfSHOOT.HISTa)

# Re-code Habitat variables
goose_3way$ceh_simp_group <- dplyr::recode(goose_3way$ceh_simp_group, arable_hort = "Arable", bog = "Bog", improved_grassland = "Improved grassland", other_grassland = "Other grassland", saltmarsh_coastal = "Saltmarsh & coastal", other = "Other")

# Re-code Shooting fix
goose_3way$PP_Shoot <- recode(goose_3way$PP_Shoot, post = "Post-shooting", pre = "Pre-shooting" )

goose_3way$PP_Shoot_f <- factor(goose_3way$PP_Shoot, levels = c("Pre-shooting", "Post-shooting"))

Fig_6a <- ggplot(goose_3way, aes(Road_Distance, fit))+
  #geom_hline(yintercept = 0.5, colour ="black", linetype = "dashed") +
  geom_line(aes(color = ceh_simp_group), position = pd, size = 1)+
  #geom_errorbar(aes(ymin=lower, ymax=upper, color = PP_Shoot), width=.4, position = pd) +
  scale_colour_manual(values = c("#CA8740", "#6FAF82", "#DABB48", "#04A3BD", "#B85D32")) +
  labs(x = "Distance from road (m)", y = "Relative probability of habitat use", color = "Habitat Type")+
  facet_grid(Species~  PP_Shoot_f) +
  theme_bw() +
  ggtitle("a)") +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", size = 0.7),
        axis.text.x = element_text(size=12), 
        axis.title=element_text(size=15),
        text = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,6000))+
  scale_y_continuous(limits = c(0,1))

Fig_6a


## 6.5 Prep data to create figure 6b ----

# Breakdown for all birds by site and case
gbg_PP <- dplyr::select(gbg_PP, -c("ceh_simp", "X.1"))
gbg_PP <- rename(gbg_PP, "dist_points" = "gbg_dist_points")
gbg_PP$Species <- "GBG"
gwf_PP$Species <- "GWfG"
gwf_PP <- rename(gwf_PP, "dist_points" = "gwf_dist_points")
goose_PP <- rbind(gbg_PP, gwf_PP)
goose_PP$sum <- 1

# Create a data frame of fixes in each habitat type by species by shooting day
PPfixdf <- with(goose_PP, aggregate(sum, by = list(Species, PP_Shoot, Used,
                                                   ceh_simp_group), "sum"))

names(PPfixdf)[1] <- "Species"
names(PPfixdf)[2] <- "PP_Shoot"
names(PPfixdf)[3] <- "Used"
names(PPfixdf)[4] <- "ceh_simp_group"
names(PPfixdf)[5] <- "number_of_fixes"

# Create a data frame of total fixes per year
PPfixtot <- with(goose_PP, aggregate(sum, by = list(Species, PP_Shoot, Used), "sum"))

names(PPfixtot)[1] <- "Species"
names(PPfixtot)[2] <- "PP_Shoot"
names(PPfixtot)[3] <- "Used"
names(PPfixtot)[4] <- "total_fixes"

# Bind data frames together using a full_join
PPdf_all <- full_join(PPfixdf, PPfixtot, by = c("Species", "PP_Shoot", "Used"))

# Order data frame by site (Species)
PPdf_all <- PPdf_all[order(PPdf_all$Species), ]

# Create a column for the proportion of fixes 
PPdf_all$ppn_fixes <- PPdf_all$number_of_fixes/PPdf_all$total_fixes

# Convert this to a percentage and create as a new column
PPdf_all$percentage_fixes <- PPdf_all$ppn_fixes*100

## 6.6 Create Figure 6b ----
# Re-code habitat classifications for figure
PPdf_all$ceh_simp_group <- dplyr::recode(PPdf_all$ceh_simp_group, arable_hort = "Arable", improved_grassland = "Improved grassland", other_grassland = "Other grassland", saltmarsh_coastal = "Saltmarsh & coastal", other = "Other")

# Re-code Shooting day
PPdf_all$PP_Shoot <- recode(PPdf_all$PP_Shoot, "post" = "Post-shooting", "pre" = "Pre-shooting" )

# Re-code GPS fix
PPdf_all$Used <- recode(PPdf_all$Used, "Avail" = "Pseudo", "Used" = "Used" )

# Create plot 
Fig_6b <- ggplot(PPdf_all, aes(x = Used, y = percentage_fixes)) +
  geom_col(aes(fill = ceh_simp_group)) +
  #scale_fill_manual()) +
  facet_grid(Species ~ PP_Shoot) +
  #labs(x = "Point Type", y = "Proportion of locations (%)", fill = "Habitat Type") +
  labs(x = "Point Type", y = "Proportion of locations (%)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#CA8740", "#6FAF82", "#DABB48","#B85D32", "#04A3BD"))+
  ggtitle("b)") +
  theme(panel.grid.major = element_blank(), 
        ##panel.grid.minor = element_blank(), 
        ##panel.border = element_blank(), 
        ##axis.line = element_line(colour = "black", size=0.7),
        ##axis.text = element_text(size=12, angle = 90, hjust = 1), 
        axis.title=element_text(size=15),
        text = element_text(size = 15),
        strip.text.x = element_text(size = 15))
Fig_6b


# 6.7 Combine plots to create Figure 6 ----
PP_Overall_habitat <- Fig_6a / Fig_6b
ggsave("Plots/Script 6) plots/Figure_6.png", plot = PP_Overall_habitat, dpi = 320, units = "in", height = 10.7, width = 8.59)

#### END ####