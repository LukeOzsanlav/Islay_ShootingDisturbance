# 26/01/2021
#### Threshold models for Luke ####

## Single Threshold Model Function ##
threshold1<-function(age,T1)
{
  # T1 threshold for age
  #
  age.1 <- age.2 <- age
  age.1[age.1 > T1] <- T1
  age.2[age.2 <= T1] <- 0
  age.2[age.2 > T1] <- age.2[age.2 > T1] - T1
  cbind(age1.1 = age.1, age1.2 = age.2)
}

## Double Threshold Model Function ##
threshold2<-function(age,T1,T2)
{
  # T1 & T2 thresholds for age
  #
  age.1 <- age.2 <- age.3 <- age
  age.1[age.1 > T1] <- T1
  age.2[age.2 <= T1] <- 0
  age.2[age.2 <= T2 & age.2 > T1] <- age.2[age.2 <= T2 & age.2 > T1] - T1
  age.2[age.2 > T2] <- T2 - T1
  age.3[age.3 <= T2] <- 0
  age.3[age.3 > T2] <- age.3[age.3 > T2] - T2
  cbind(age2.1 = age.1, age2.2 = age.2, age2.3 = age.3)
}

#### glmer models ####

# Make empty outputs to store model AICs
summary(Mums) # Mum ages are 2-14
output <- data.frame(model=as.numeric(),AIC=numeric()) # This is for single threshold model 
output2 <- matrix(NA,13,7) # This is for double threshold model 
names <- matrix(NA,13,7) # This is for double threshold model


# Single threshold model 
jj<-seq(3,13,1)
for(i in 1:length(jj)) {
  model <- glmer(number.offspring ~ threshold1(mum.age,jj[i]) + SocSize + LYO + scale(ALC) + (1|mum.ID) + (1|socg) + (1|year), data=Mums3, family = "poisson", REML=FALSE)
  output[i,]<-c(jj[i],AIC(model))
}

output # AIC 1963.987 for threshold at mum age of 3 years 

# Double threshold model
# I've chosen 3-7 years and 8-13 years as the ranges that threshold 1 and threshold 2 can sit in. Use common sense/badger biology for this decision ;)
for(i in 3:7){
  for (j in 8:13){
    model <- glmer(number.offspring ~ threshold2(mum.age,i,j) + SocSize + LYO + scale(ALC) + (1|mum.ID) + (1|socg) + (1|year), data=Mums3, family = "poisson", REML=FALSE)
    output2[j,i]<-AIC(model)
    names[j,i]<-paste(i,j)
  }}

# Check AIC values to get best fitting threshold point/s
output # AIC 1967.472 age 3
output2 # AIC 1964.567 ages 3 and 8   

#### Plotting threshold models ####
# For reference Luke, I've used a double threshold model here but you'll get the idea. Also predict() starting having errors with random effects so I've cheated and just used TMB without random effects for the purposes of plotting. 


# TOP MODEL: glmer(number.offspring ~ threshold2(mum.age,3,8) + SocSize + LYO + scale(ALC) + (1|mum.ID) + (1|socg) + (1|year), family = "poisson", data=Mums3, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
library(glmmTMB)
top.model <- glmmTMB(number.offspring ~ threshold2(mum.age,3,8) + SocSize + LYO + scale(ALC), family = "poisson", data=Mums3)

summary(Mums2)
Mums$mum.ID <- as.factor(Mums$mum.ID)
Mums$socg <- as.factor(Mums$socg)

newdata <- expand.grid(mum.age  = seq(2, 14, length = 100),
                       SocSize = 9.86,
                       LYO = "0",
                       ALC = 2933) # This is the scale(meanALC) value

pred <- predict(top.model, newdata, re.form = NA, se = TRUE, type = "response")

newdata$newy <- pred[[1]]  # subset first bit in list (fit in pred)
newdata$lwrCI <- pred[[1]]-(1.96*pred[[2]]) # lower confidence interval (1.96xSE)
newdata$uprCI <- pred[[1]]+(1.96*pred[[2]]) # upper confidence interval

ggplot(data = newdata,
       aes(x = mum.age,
           y = newy))+
  geom_line(lwd=1, lty=1, colour = "red")+
  geom_jitter(data=Mums,
              aes(x=mum.age,
                  y=number.offspring), size=1.3, alpha = 0.2)+
  geom_ribbon(aes(ymin=lwrCI,
                  ymax=uprCI),
              alpha = 0.2, fill= "red")+
  scale_x_continuous(breaks = seq(2,14, by = 1))+
  xlab("Maternal age (years)")+
  ylab("Annual offspring production")+
  theme_classic()+
  theme(axis.line = element_line(size= 1, colour = "black"),
        axis.text.x = element_text(colour = "black", margin= margin(t=5)),
        axis.text.y = element_text(colour = "black", margin = margin(r=5)),
        axis.title.x = element_text(size=18,vjust = 0),
        axis.title.y = element_text(size=18, vjust = 0.5, angle = 90, margin = margin(r=10)),
        axis.text = element_text(size=18),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(2, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())