##---------------------------------------------------------##
## Aimée McIntosh
## Created 27/10/23

## Aims:
## Model whether individual daily movement distances differs between shooting disturbed and non-disturbed days
##
## Prediction 3: 
## Daily movement distances will be greater on days when individuals are exposed to shooting
##---------------------------------------------------------##


## Packages required
pacman::p_load(tidyr,plyr,dplyr,zoo,data.table,move,ggplot2,purrr,readr,
               lubridate,tidyverse,readr,amt,geosphere,sf,DHARMa,lme4,
               MuMIn,effects,lattice,performance,MASS,glmmTMB,bbmle,emmeans, patchwork)



#--------------------------------#
####  1. Read and prep data   ####
#--------------------------------#
distance_data <- read.csv("Foraging_Distance_Model_Data.csv")
distance_data$shoot_date <- as.factor(distance_data$shoot_date)
distance_data$Tag_Winter <- as.factor(distance_data$Tag_Winter)
distance_data$group_id <- as.numeric(distance_data$Tag_Winter)
distance_data$Winter <- as.factor(distance_data$Winter)
head(distance_data)

#------------------------------------#
####  2.Model 2-level disturbance ####
#------------------------------------#
# re-level factor to GWfG as base line

dist_glob <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                       Disturbance2*Species + Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) + 
                       Disturbance2*Species*scale(indiv_dist) + ar1(shoot_date + 0|group_id) +
                       (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                     data = distance_data,
                     family = Gamma("log"), na.action = "na.fail")

dist_globb <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) + 
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globc <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Disturbance2*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globd <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globe <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*scale(indiv_dist) + Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globf <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*Species + ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"))

dist_globg <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Disturbance2*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globh <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +  
                        Species*scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globi <- glmmTMB(totaldist ~ Disturbance2 + Species + scale(indiv_dist) + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globj <- glmmTMB(totaldist ~ Disturbance2 + Species + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globk <- glmmTMB(totaldist ~ Species + Sex + scale(indiv_dist) +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

dist_globl <- glmmTMB(totaldist ~ Species + Sex +
                        ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

#------------------------------------#
####  3.Model 3-level disturbance ####
#------------------------------------#
shoot_glob <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                        Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) + 
                        Disturbance_shoot*Species*scale(indiv_dist) + ar1(shoot_date + 0|group_id) +
                        (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                      data = distance_data,
                      family = Gamma("log"), na.action = "na.fail")

shoot_globb <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) + 
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globc <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Disturbance_shoot*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globd <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globe <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*scale(indiv_dist) + Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globf <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*Species + ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"))

shoot_globg <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Disturbance_shoot*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globh <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +  
                         Species*scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globi <- glmmTMB(totaldist ~ Disturbance_shoot + Species + scale(indiv_dist) + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globj <- glmmTMB(totaldist ~ Disturbance_shoot + Species + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globk <- glmmTMB(totaldist ~ Species + Sex + scale(indiv_dist) +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

shoot_globl <- glmmTMB(totaldist ~ Species + Sex +
                         ar1(shoot_date + 0|group_id) +
                         (1|Winter) + (1|Device_ID) + (1|mode_hld_farm), 
                       data = distance_data,
                       family = Gamma("log"), na.action = "na.fail")

#------------------------------------#
####  4. Model Checks ####
#------------------------------------#

## 4.1 Calculate and compare models with AICc ----
ICtab(dist_glob,dist_globb,dist_globc,dist_globd,dist_globe,dist_globf,dist_globg,dist_globh,
      dist_globi, dist_globj, dist_globk, dist_globl,shoot_glob,
      shoot_globb,shoot_globc,shoot_globd,shoot_globe,shoot_globf,
      shoot_globg,shoot_globh,shoot_globi, shoot_globj, shoot_globk, shoot_globl,
      weights = TRUE, 
      delta = TRUE, 
      base = TRUE)

## 4.2 DHARMa model checks ----
# Simulate residuals and plot - this is for best model can check for all models

simtop1 <- simulateResiduals(dist_glob)
plot(simtop1)

## 4.3 Calculate r-squared ----
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
r.squaredGLMM(shoot_glob)
r.squaredGLMM(shoot_globb)
r.squaredGLMM(shoot_globc)
r.squaredGLMM(shoot_globd)
r.squaredGLMM(shoot_globe)
r.squaredGLMM(shoot_globf)
r.squaredGLMM(shoot_globg)
r.squaredGLMM(shoot_globh)
r.squaredGLMM(shoot_globi)
r.squaredGLMM(shoot_globj)
r.squaredGLMM(shoot_globk)
r.squaredGLMM(shoot_globl)

#-------------------------------------#
####         5.Model Plots         ####
#-------------------------------------#

## 5.1 Extract effects from models ----
df_dist_top <- as.data.frame(effect("Disturbance2:Species:scale(indiv_dist)",dist_glob), xlevels = 35) # Specific interaction terms
effects_top <- as.data.frame(allEffects(dist_glob)) # All effects

df_dist_top2a <- as.data.frame(effect("Disturbance2:scale(indiv_dist)",dist_glob, xlevels = 35)) # Specific interaction terms
effects_top2a <- as.data.frame(allEffects(dist_glob)) # All effects

df_dist_top2b <- as.data.frame(effect("Disturbance2:Species",dist_glob)) # Specific interaction terms
effects_top2b <- as.data.frame(allEffects(dist_glob)) # All effects

confint(dist_glob) # Confidence intervals of top model

# Visualise all effects in top model
plot(allEffects(dist_glob2))

## 5.2 Plot Figure 4a ----
# This is plotted for the best overall performing model
df_dist_top$Disturbance <- recode(df_dist_top$Disturbance2, No_Shoot = "No Shooting", Shoot = "Shooting" )

Fig_4a <- ggplot(df_dist_top, aes(indiv_dist, fit))+
  geom_line(aes(colour = Disturbance), size = 1)+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Disturbance), alpha = 0.3) +
  #scale_colour_manual(values = c("#44AA99", "#DDCC77", "#999933", "#117733", "#CC6677", "#882255", "#AA4499")) +
  labs(x = "Cumulative shooting experience", y = "Daily dsitance travelled (km)", colour = "Disturbance Day", fill = "Disturbance Day")+
  scale_colour_manual(values = c("No Shooting" = "#708090", "Shooting" = "#CD5C5C")) +
  scale_fill_manual(values = c("No Shooting" = "#708090", "Shooting" = "#CD5C5C")) +
  facet_grid(~Species) +
  theme_bw() +
  ggtitle("a) Three-way model interaction") +
  theme(plot.title = element_text(size=12)) +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", size = 0.7),
        axis.text.x = element_text(size=10), 
        axis.title=element_text(size=10),
        text = element_text(size = 10)) +
  scale_x_continuous(limits = c(0,35)) +
  scale_y_continuous(limits = c(0,17),  breaks = seq(0,17, by = 2))

Fig_4a
ggsave("Disturbance_Paper_Files/Daily travel distance_3-way - 2level disturbance - model.jpeg", Top_model_plot)

## 5.3 Plot Figure 4b ----
# This is the overall top model disturbance interaction
df_dist_top2b$Disturbance <- recode(df_dist_top2b$Disturbance2, No_Shoot = "No Shooting", Shoot = "Shooting" )

Disturbance_labels <- c("No Shooting", "Shooting")
pd <- position_dodge(0.43)

Fig_4b <- ggplot(df_dist_top2b, aes(Disturbance2, fit)) +
  geom_point(aes(colour = Species), position = pd, size = 3) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Species), width=.4, position = pd) +
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(x = " Disturbance Day", y = "Daily distance travelled (km)", colour = "Species")+
  scale_x_discrete(labels = Disturbance_labels) +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, by = 1)) +
  theme_bw() +
  ggtitle("b) Two-way model interaction") +
  theme(plot.title = element_text(size=12)) +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", size = 0.7),
        axis.text.x = element_text(size=10), 
        axis.title=element_text(size=10),
        text = element_text(size = 12)) 
Fig_4b

ggsave("Disturbance_Paper_Files/Daily travel distance_3-way - 2level disturbance - model.jpeg", Top_model_plot)

## 5.4 Plot overall Figure 4
# Combine all the plots together
Figure_5 <- Fig_4a / Fig_4b
Figure_5

ggsave("Figure_5.jpeg", Figure_5)


## 5.5 Plot of individual experience ----
df_dist_top2a$Disturbance <- recode(df_dist_top2a$Disturbance2, No_Shoot = "No Shooting", Shoot = "Shooting" )

Disturbance_labels <- c("No Shooting", "Shooting")
pd <- position_dodge(0.43)

Top_model_plot_experience <- ggplot(df_dist_top2a, aes(indiv_dist, fit))+
  geom_line(aes(colour = Disturbance), size = 1)+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Disturbance), alpha = 0.3) +
  #scale_colour_manual(values = c("#44AA99", "#DDCC77", "#999933", "#117733", "#CC6677", "#882255", "#AA4499")) +
  labs(x = "Cumulative shooting experience", y = "Daily dsitance travelled (km)", colour = "Disturbance")+
  scale_colour_manual(values = c("No Shooting" = "#708090", "Shooting" = "#CD5C5C")) +
  scale_fill_manual(values = c("No Shooting" = "#708090", "Shooting" = "#CD5C5C")) +
  #  facet_grid(~Species) +
  theme_bw() +
  ggtitle("b) Two-way model interaction") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", size = 0.7),
        axis.text.x = element_text(size=10), 
        axis.title=element_text(size=10),
        text = element_text(size = 10)) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size=12))+
  scale_x_continuous(limits = c(0,35)) +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, by = 1))

Top_model_plot_experience

#### End ####