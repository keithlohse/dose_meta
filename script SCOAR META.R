#By Keith Lohse, Rehabilitation Informatics Lab, 2016-04-26

# For this analysis, you will need to install and then open the following packages:
# install.packages("metafor"); install.packages("dplyr"); install.packages("ggplot2")
library("metafor"); library("dplyr"); library("ggplot2")

# Read in the full data set
DATA<-read.table("data SCOAR TEXT OUTLIERS REMOVED.txt", header = TRUE, sep="\t") 
head(DATA)

# Alternatively, you can read the text file into R directly from my GitHub repo:
# DATA<-read.table("https://raw.github.com/keithlohse/dose_meta/master/data%20SCOAR%20TEXT%20OUTLIERS%20REMOVED.txt",
#    header = TRUE, sep="\t")
# head(DATA)

# Note that all effect sizes greater than d = 3.0 have been removed

##########################################
### Analysis of the primary extraction ###
##########################################

# In this first analysis, we took only the primary outcome or
# the first usable secondary outcome if no primary was listed.
# Below, we will conduct separate analyses for the FMA and gait
# speed (regardless of if those measures were primary outcomes).

# Create different sets of data
# 1. DATA <- contains all data for experimental and control groups
# 2. LOHSE <- contains only primary outcomes and excludes missing cases
# 3. BIGN <- contains only primary outcomes for groups with base n > 30
# 3. CTRLS <- contains only the control group data from LOHSE
# 4. EXPS <- contains on the experimental group data from LOHSE
# 5. FMA <- contains data for Fugl-Meyer Assessment outcomes excludes missing cases
# 6. SPEED <- contains data for gait speed outcomes excludes missing cases

# Creating the "LOHSE" subset
LOHSE<-subset(DATA, SCOAR_outcome == "lohse")
# Length should be 489 (total number of independent groups)
length(LOHSE$SCOAR_outcome)
# Sum should be 12,846 (total number of participants)
sum(LOHSE$base_n)

# Remove studies missing time for therapy
LOHSE<-subset(LOHSE,!time_50 == "na") 
# Remove studies missing duration of therapy
LOHSE<-subset(LOHSE,!exp_dur == "na") 
# Remove studies missing age at baseline
LOHSE<-subset(LOHSE,!age_base == "na")
# Remove studies missing days post-stroke
LOHSE<-subset(LOHSE,!days_ps == "na") 
# Remove studies missing terminal effect-sizes
LOHSE<-subset(LOHSE,!term_d == "na") 
# New length should be 303
length(LOHSE$SCOAR_outcome)
# New sum should be 6767
sum(LOHSE$base_n)

# Creating a subset of only those studies with large sample sizes
BIGN<-subset(LOHSE, base_n >= 20)
length(BIGN$SCOAR_outcome)
sum(BIGN$base_n)

# Creating a subset of control groups only
CTRLS<-subset(LOHSE, group == "ctrl")
summary(CTRLS$group)

# Creating a subset of experimental groups only
EXPS<-subset(LOHSE, group == "exp")
summary(EXPS$group)

# FMA and SPEED datasets are defined below

#################################################
# Descriptive Statistics of the included Groups #
#################################################
# average age of patients at baseline
summary(LOHSE$age_base) 
sd(LOHSE$age_base)
# average days from stroke to the beginning of therapy
summary(LOHSE$days_ps) 
sd(LOHSE$days_ps)
# duration of the intervention
summary(LOHSE$exp_dur) 
sd(LOHSE$exp_dur)
# time scheduled for therapy based on the Max Time calculation
summary(LOHSE$time_MAX) 
sd(LOHSE$time_MAX)
# time scheduled for therapy based on the 50% time calculation
summary(LOHSE$time_50) 
sd(LOHSE$time_50)
# time scheduled for therapy based on the Min time calculation
summary(LOHSE$time_MIN) 
sd(LOHSE$time_MIN)
# How was time in therapy quantified? 1 = time scheduled; 4 = repetitions
summary(as.factor(LOHSE$detailed_time)) 
# Did the authors specify an intention to treat analysis
summary(as.factor(LOHSE$itt_analysis)) 
# Number of subjects contributing to baseline means
summary(LOHSE$base_n) 
sd(LOHSE$base_n)
# Number of subjects contributing to the terminal means
summary(LOHSE$term_n) 
sd(LOHSE$term_n)
# Number of subjects contributing for follow up means
summary(LOHSE$fu_n) 
sd(LOHSE$fu_n, na.rm=TRUE) 
# Note that we did not filter the data by follow up outcomes,
# so there are missing cases that we need to exclude in this calculation

# Time from the baseline assessment to the follow-up assessment
summary(LOHSE$fu_dur) 
sd(LOHSE$fu_dur, na.rm=TRUE)
# Number of experimental versus control groups
summary(LOHSE$group) 
# Number of subjects contributing to the terminal means
summary(LOHSE$term_g) 
sd(LOHSE$term_g)


##############################
# Analysis of All Group Data #
##############################

# Random Effect model of all groups, using maximum likelihood estimation
m0<-rma(term_g, term_Vg, data=LOHSE, method="ML")
m0
confint(m0)
# Random Effect model of CONTROL groups, using maximum likelihood estimation
m0A<-rma(term_g, term_Vg, data=CTRLS, method="ML") 
m0A
# Random Effect model of EXPERIMENTAL groups, using maximum likelihood estimation
m0B<-rma(term_g, term_Vg, data=EXPS, method="ML") 
m0B


# Creating a forest plot to show the RE model of all of the data
forest(m0, slab=paste(LOHSE$author, LOHSE$year, sep=", "), cex=1.5)

# Creating a funnel plot to show potential bias in the full dataset
palette(c("orange","dodgerblue"))
funnel(m0, pch=21, bg=LOHSE$group, cex = 1.5, cex.axis=1.25, cex.lab=1.25, xlab="Terminal Outcome (g)")
# Statistical test of symmetry (see Eggers, 1997)
regtest(m0, model = "lm") 
# Significant result indicates bias in the distribution of results

# Moderator Analysis of Treatment Versus Control Groups
## First we need to center the "group" variable
summary(as.numeric(LOHSE$group))
LOHSE$exp.c<-as.numeric(LOHSE$group)-1.5
head(LOHSE)
## Then we can included this centered group variable in our model
m1<-rma(term_g~exp.c, term_Vg, data=LOHSE, method="ML")
m1
## We can also get the 95% confidence estimates for the model
confint(m1)

#Moderator Analysis of Groups and DPS
## First we center the DPS variable
summary(LOHSE$days_ps)
LOHSE$DPS.c<-LOHSE$days_ps-mean(LOHSE$days_ps)
hist(LOHSE$days_ps)

## Following our transformation, we now center the rtDPS variable
LOHSE$rtDPS<-sqrt(LOHSE$days_ps)
summary(LOHSE$rtDPS)
LOHSE$rtDPS.c<-LOHSE$rtDPS-mean(LOHSE$rtDPS)
hist(LOHSE$rtDPS.c)
## We then include the centered rtDPS variable in our model
m2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=LOHSE, method="ML")
m2
confint(m2)

#Moderator Analysis of Group and Estimated Time Scheduled
## First we mean center the TIME variable
summary(LOHSE$time_50)
hist(LOHSE$time_50)
LOHSE$time50.c<-LOHSE$time_50-mean(LOHSE$time_50)
## Following our transformation, we now mean center the rtTIME variable
LOHSE$rtTIME<-sqrt(LOHSE$time_50)
summary(LOHSE$rtTIME)
LOHSE$rtTIME.c<-LOHSE$rtTIME-mean(LOHSE$rtTIME)
hist(LOHSE$rtTIME)
## We can then add the mean centered rtTIME to your model
m3<-rma(term_g~exp.c+rtDPS.c+rtTIME.c, term_Vg, data=LOHSE, method="ML")
m3
confint(m3)
plot(m3)

# Interaction of Dose with Time
## No new variables need to be created to test this interaction
m4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=LOHSE, method="ML")
m4
confint(m4)
plot(m4)

#Moderator Analysis of Group, DPS, TIME and Age
summary(LOHSE$age_base)
hist(LOHSE$age_base)
## First we mean center the AGE variable
LOHSE$age.c<-LOHSE$age_base-mean(LOHSE$age_base)
## Then we can include mean centered age in our model
m5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=LOHSE, method="ML")
m5
confint(m5)


#Controlling for the effect of therapy duration
summary(LOHSE$exp_dur)
hist(LOHSE$exp_dur)
## First we mean center the variable of duration
LOHSE$DUR.c<-LOHSE$exp_dur-mean(LOHSE$exp_dur)
## Following our transformation, we now mean center the rtDUR variable
LOHSE$rtDUR<-sqrt(LOHSE$exp_dur)
summary(LOHSE$rtDUR)
hist(LOHSE$rtDUR)
LOHSE$rtDUR.c<-LOHSE$rtDUR-mean(LOHSE$rtDUR)
## We also want to check for potential colinearity between duration and time
plot(LOHSE$rtDUR.c,LOHSE$rtTIME.c)
cor.test(LOHSE$rtDUR.c,LOHSE$rtTIME.c)
## Finally we can add the mean centered rtDUR variable to the model
m6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=LOHSE, method="ML")
m6
confint(m6)

# Next, we export the data (including the variables we created) to a csv file
write.csv(LOHSE, file="data SCOAR R MAIN ANALYSIS.csv")

#Comparison between the different models
summary(m0)
# logLik   deviance        AIC        BIC       AICc  
# -196.4562   666.8507   396.9124   404.3398   396.9524
summary(m1)
# logLik   deviance        AIC        BIC       AICc  
# -179.4614   632.8610   364.9228   376.0640   365.0030  
summary(m2)
# logLik   deviance        AIC        BIC       AICc  
# -133.6918   541.3218   275.3836   290.2385   275.5178  
summary(m3)
# logLik   deviance        AIC        BIC       AICc  
# -127.3528   528.6440   264.7057   283.2744   264.9077
summary(m4)
# logLik   deviance        AIC        BIC       AICc  
# -125.4686   524.8755   262.9372   285.2196   263.2210 
summary(m5)
# logLik   deviance        AIC        BIC       AICc  
# -121.7344   517.4070   257.4687   283.4649   257.8484  
summary(m6)
# logLik   deviance        AIC        BIC       AICc  
# -120.4646   514.8674   256.9291   286.6390   257.4189  

# So m5 is our best fitting model...
# but lets double-check our statistical assumptions/diagnostics
# The distribution of residuals for this model is approximately normal
plot(density(resid(m5))) 
# The distribution of residuals for this model show some heteroscedasticity
plot(resid(m5)~fitted(m5)) 
# but not a significant relationship between fitted values and residuals.
cor.test(resid(m5),fitted(m5))


#########################################
# Plots for Multivariable Relationships #
#########################################
palette(c("purple","chartreuse","cyan"))

# Plot of outcomes by DOSE and DPS
plot(LOHSE$term_g~LOHSE$rtTIME, bty='n', type='p', pch=21, bg=LOHSE$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,
     ylab="Terminal Outcome (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")
## Cyan = <90 days post stroke, green = <365 days post stroke, purple = > 365 days post stroke


# Plot of outcomes by AGE
# the ntile function chops the variable into X equally sized groups.
LOHSE$AGEq <- ntile(LOHSE$age_base, 4)  
# Topocolors is just a color pallette that I like
palette(topo.colors(4)) 

plot(LOHSE$term_g~LOHSE$age_base, bty='n', xlim=c(40,90), ylim=c(-0.5,3.5), bg=LOHSE$AGEq, col = "black",
     pch=21, cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,  
     ylab="Terminal Outcome (g)", xlab="Average Age (years)")

#Plot of outcomes by DPS
plot(LOHSE$term_g~LOHSE$yrs_ps, bty='n', type='p', pch=21, bg=LOHSE$days_cat, col="black",  
     lwd=1.0, xlim=c(0,10), ylim=c(-0.5,3.5), cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,
     ylab="Terminal Outcome (g)", xlab="Time Post Stroke (yrs)")

plot(LOHSE$term_g~LOHSE$rtDPS, bty='n', type='p', pch=21, bg=LOHSE$days_cat, col="black",  
     lwd=1.0, xlim=c(0,60), ylim=c(-0.5,3.5), cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,
     ylab="Terminal Outcome (g)", xlab="Time Post Stroke (sqrt(days)")
#The square root of DPS is the actual variable in our models, so I think we should use it in the figures

#Plot of outcomes by Time Scheduled for therapy
LOHSE$TIMEq <- ntile(LOHSE$rtTIME, 4)  
palette(topo.colors(4))

plot(LOHSE$term_g~LOHSE$rtTIME, bty='n', xlim=c(0,15), ylim=c(-0.5,3.5), bg=LOHSE$TIMEq, col = "black",
     pch=21, cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,  
     ylab="Terminal Outcome (g)", xlab="Estimated Time Scheduled for Therapy (sqrt(hrs))")

plot(LOHSE$term_g~LOHSE$time_50, bty='n', xlim=c(0,250), ylim=c(-0.5,3.5), bg=LOHSE$TIMEq, col = "black",
     pch=21, cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,  
     ylab="Terminal Outcome (g)", xlab="Estimated Time Scheduled for Therapy (hrs))")



#################################################
### Analysis of the "Big" groups only, n > 20 ###
#################################################
# Descriptive Statistics for the BIGN data
summary(BIGN$age_base)
sd(BIGN$age_base, na.rm=TRUE)

summary(BIGN$days_ps)
sd(BIGN$days_ps, na.rm=TRUE)

summary(BIGN$exp_dur)
sd(BIGN$exp_dur, na.rm=TRUE)

summary(BIGN$time_50)
sd(BIGN$time_50, na.rm=TRUE)

summary(BIGN$group)

summary(BIGN$base_n)
sd(BIGN$base_n, na.rm=TRUE)

summary(BIGN$term_g)
sd(BIGN$term_g, na.rm=TRUE)

# Overall analysis of all trials
b0<-rma(term_g, term_Vg, data=BIGN, method="ML")
b0
confint(b0)

# Creating a forest plot to show the RE model of all of the data
forest(b0, slab=paste(BIGN$author, BIGN$year, sep=", "), cex=1.5)

# Creating a funnel plot to show potential bias in the full dataset
palette(c("orange","dodgerblue"))
funnel(b0, pch=21, bg=BIGN$group, cex = 1.5, cex.axis=1.25, cex.lab=1.25, xlab="Terminal Outcome (g)")
# Statistical test of symmetry (see Eggers, 1997)
regtest(m0, model = "lm") 
# Significant result indicates bias in the distribution of results

# Moderator Analysis of Treatment Versus Control Groups
summary(as.numeric(BIGN$group))
BIGN$exp.c<-as.numeric(BIGN$group)-1.5
head(BIGN)
b1<-rma(term_g~exp.c, term_Vg, data=BIGN, method="ML")
b1
confint(b1)

# Moderator Analysis of Groups and DPS
summary(BIGN$days_ps)
BIGN$DPS.c<-BIGN$days_ps-mean(BIGN$days_ps)
hist(BIGN$days_ps)

BIGN$rtDPS<-sqrt(BIGN$days_ps)
summary(BIGN$rtDPS)
BIGN$rtDPS.c<-BIGN$rtDPS-mean(BIGN$rtDPS)
hist(BIGN$rtDPS.c)

b2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=BIGN, method="ML")
b2
confint(b2)

# Moderator Analysis of Group and Estimated Time Scheduled
summary(BIGN$time_50)
hist(BIGN$time_50)
BIGN$time50.c<-BIGN$time_50-mean(BIGN$time_50)

BIGN$rtTIME<-sqrt(BIGN$time_50)
summary(BIGN$rtTIME)
BIGN$rtTIME.c<-BIGN$rtTIME-mean(BIGN$rtTIME)
hist(BIGN$rtTIME)

b3<-rma(term_g~exp.c+rtDPS.c+rtTIME.c, term_Vg, data=BIGN, method="ML")
b3
confint(b3)
plot(b3)

# Adding the interaction of DPS and TIME
b4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=BIGN, method="ML")
b4
confint(b4)
plot(fitted(b4),resid(b4))
plot(density(resid(b4)))

# Moderator Analysis of Group, DPS, TIME and Age
summary(BIGN$age_base)
hist(BIGN$age_base)
BIGN$age.c<-BIGN$age_base-mean(BIGN$age_base)

b5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=BIGN, method="ML")
b5
confint(b5)


# Controlling for the effect of therapy duration
summary(BIGN$exp_dur)
hist(BIGN$exp_dur)
BIGN$DUR.c<-BIGN$exp_dur-mean(BIGN$exp_dur)

BIGN$rtDUR<-sqrt(BIGN$exp_dur)
summary(BIGN$rtDUR)
BIGN$rtDUR.c<-BIGN$rtDUR-mean(BIGN$rtDUR)
hist(BIGN$rtDUR.c)
## Testing for colinearity between duration and time
plot(BIGN$rtDUR.c,BIGN$rtTIME.c)
cor.test(BIGN$rtDUR.c,BIGN$rtTIME.c)
# Controlling for the duration of therapy
b6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=BIGN, method="ML")
b6
confint(b6)

# Comparing the fit of the different models
summary(b0)
# logLik  deviance       AIC       BIC      AICc  
# -70.2411  325.1712  144.4821  150.0902  144.5830  
summary(b1)
# logLik  deviance       AIC       BIC      AICc  
# -65.2092  315.1075  136.4185  144.8306  136.6219 
summary(b2)
# logLik  deviance       AIC       BIC      AICc  
# -39.8544  264.3979   87.7088   98.9249   88.0507 
summary(b3)
# logLik  deviance       AIC       BIC      AICc  
# -36.3532  257.3954   82.7064   96.7265   83.2236 
summary(b4)
# logLik  deviance       AIC       BIC      AICc  
# -34.4962  253.6814   80.9924   97.8165   81.7228  
summary(b5)
# logLik  deviance       AIC       BIC      AICc  
# -33.5010  251.6911   81.0020  100.6302   81.9845
summary(b6)
# logLik  deviance       AIC       BIC      AICc  
# -32.1270  248.9430   80.2540  102.6861   81.5283 

# So b4 is our best fitting model...
# but lets double-check our statistical assumptions/diagnostics
# The distribution of residuals for this model is approximately normal
plot(density(resid(b4))) 
# The distribution of residuals for this model show some heteroscedasticity
plot(resid(b4)~fitted(b4)) 
# but not a significant relationship between fitted values and residuals.
cor.test(resid(b4),fitted(b4))

# We can also export the BIGN data (and the variables we created) to a csv
write.csv(BIGN, file="data SCOAR BIG N ANALYSIS.csv")

# Plot of outcomes by AGE
# the ntile function chops the variable into X equally sized groups.
BIGN$AGEq <- ntile(BIGN$age_base, 4)  
# Topocolors is just a color pallette that I like
palette(topo.colors(4)) 

plot(BIGN$term_g~BIGN$age_base, bty='n', xlim=c(40,90), ylim=c(-0.5,3.5), bg=BIGN$AGEq, col = "black",
     pch=21, cex=sqrt(BIGN$base_n)/2, cex.axis=1.25, cex.lab=1.25,  
     ylab="Terminal Outcome (g)", xlab="Average Age (years)")

#Plot of outcomes by DPS
plot(BIGN$term_g~BIGN$rtDPS, bty='n', type='p', pch=21, bg=BIGN$days_cat, col="black",  
     lwd=1.0, xlim=c(0,60), ylim=c(-0.5,3.5), cex=sqrt(BIGN$base_n)/2, cex.axis=1.25, cex.lab=1.25,
     ylab="Terminal Outcome (g)", xlab="Time Post Stroke (sqrt(days)")
#The square root of DPS is the actual variable in our models, so I think we should use it in the figures

#Plot of outcomes by Time Scheduled for therapy
BIGN$TIMEq <- ntile(BIGN$rtTIME, 4)  
palette(topo.colors(4))

plot(BIGN$term_g~BIGN$rtTIME, bty='n', xlim=c(0,15), ylim=c(-0.5,3.5), bg=BIGN$TIMEq, col = "black",
     pch=21, cex=sqrt(BIGN$base_n)/2, cex.axis=1.25, cex.lab=1.25,  
     ylab="Terminal Outcome (g)", xlab="Estimated Time Scheduled for Therapy (sqrt(hrs))")




###########################################
### Analysis of the FMA for the UE only ###
###########################################

# All outcomes are fma, fma-ue, fma-se, or fma-wh
# fma-le had been excluded
# only one outcome per study
FMA<-subset(DATA, fma_outcome == "fma")
# Remove studies missing time for therapy
FMA<-subset(FMA,!time_50 == "na") 
# Remove studies missing duration of therapy
FMA<-subset(FMA,!exp_dur == "na") 
# Remove studies missing age at baseline
FMA<-subset(FMA,!age_base == "na") 
# Remove studies missing days post-stroke
FMA<-subset(FMA,!days_ps == "na") 
# Remove studies missing terminal effect-sizes
FMA<-subset(FMA,!term_d == "na") 
length(FMA$author)
# Length for the FMA data should be 79
sum(FMA$base_n)
# Sum of baseline n should be 1,741

# Descriptive Statistics for the FMA data
summary(FMA$age_base)
sd(FMA$age_base, na.rm=TRUE)

summary(FMA$days_ps)
sd(FMA$days_ps, na.rm=TRUE)

summary(FMA$exp_dur)
sd(FMA$exp_dur, na.rm=TRUE)

summary(FMA$time_50)
sd(FMA$time_50, na.rm=TRUE)

summary(FMA$group)

summary(FMA$base_n)
sd(FMA$base_n, na.rm=TRUE)

summary(FMA$term_g)
sd(FMA$term_g, na.rm=TRUE)

#Regression models for FMA
f0<-rma(term_g, term_Vg, data=FMA, method="ML")
f0

# Moderator Analysis of Group
summary(FMA$group)
FMA$exp.c<-as.numeric(FMA$group)-1.5
head(FMA)
f1<-rma(term_g~exp.c, term_Vg, data=FMA, method="ML")
f1

# Moderator Analysis of Groups and DPS
summary(FMA$days_ps)
FMA$DPS.c<-FMA$days_ps-mean(FMA$days_ps)
hist(FMA$DPS.c)

FMA$rtDPS<-sqrt(FMA$days_ps)
summary(FMA$rtDPS)
FMA$rtDPS.c<-FMA$rtDPS-mean(FMA$rtDPS)
hist(FMA$rtDPS.c)
mean(FMA$rtDPS.c)

f2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=FMA, method="ML")
f2

# Moderator Analysis of Group and Estimated Time Scheduled
summary(FMA$time_50)
hist(FMA$time_50)
FMA$time50.c<-FMA$time_50-mean(FMA$time_50)

FMA$rtTIME<-sqrt(FMA$time_50)
summary(FMA$rtTIME)
FMA$rtTIME.c<-FMA$rtTIME-mean(FMA$rtTIME)
hist(FMA$rtTIME)

f3<-rma(term_g~exp.c+rtDPS.c+rtTIME.c, term_Vg, data=FMA, method="ML")
f3
plot(f3)

# Adding the interaction of DPS and TIME
f4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=FMA, method="ML")
f4
confint(f4)
plot(f4)

# Moderator Analysis of Group, DPS, and Age
summary(FMA$age_base)
hist(FMA$age_base)

FMA$age.c<-FMA$age_base-mean(FMA$age_base)

f5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=FMA, method="ML")
f5
confint(f5)

# Controlling for the effect of therapy duration
summary(FMA$exp_dur)
hist(FMA$exp_dur)
FMA$DUR.c<-FMA$exp_dur-mean(FMA$exp_dur)

FMA$rtDUR<-sqrt(FMA$exp_dur)
summary(FMA$rtDUR)
hist(FMA$rtDUR)
FMA$rtDUR.c<-FMA$rtDUR-mean(FMA$rtDUR)
## Checking for colinearity between TIME and DURATION
plot(FMA$rtDUR.c,FMA$rtTIME.c)
cor.test(FMA$rtDUR.c,FMA$rtTIME.c)

f6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=FMA, method="ML")
f6

# Comparisons between the different FMA models
summary(f0)
# logLik  deviance       AIC       BIC      AICc  
# -26.1594  129.2677   56.3187   61.0576   56.4766 
summary(f1)
# logLik  deviance       AIC       BIC      AICc  
# -22.5504  122.0497   51.1007   58.2091   51.4207
summary(f2)
# logLik  deviance       AIC       BIC      AICc  
# -13.7268  104.4027   35.4537   44.9315   35.9942  
summary(f3)
# logLik  deviance       AIC       BIC      AICc  
# -13.5434  104.0359   37.0869   48.9341   37.9088
summary(f4)
# logLik  deviance       AIC       BIC      AICc  
# -13.5408  104.0307   39.0817   53.2984   40.2484
summary(f5)
# logLik  deviance       AIC       BIC      AICc  
# -10.2347   97.4185   34.4695   51.0556   36.0469 
summary(f6)
# logLik  deviance       AIC       BIC      AICc  
# -9.5106   95.9703   35.0213   53.9769   37.0784

# So f5 is our best fitting model based on the AIC and AICc
# but lets double-check our statistical assumptions/diagnostics
# The distribution of residuals for this model is approximately normal
plot(density(resid(f5))) 
# The distribution of residuals for this model show some heteroscedasticity
plot(resid(f5)~fitted(f5)) 
# but not a significant relationship between fitted values and residuals.
cor.test(resid(f5),fitted(f5))

# Exporting the FMA data to a .csv file
write.csv(FMA, file="data SCOAR FMA ANALYSIS.csv")

# Plots of the FMA data
palette(c("purple","chartreuse","cyan"))

plot(FMA$term_g~FMA$rtTIME, bty='n', type='p', pch=21, bg=FMA$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex.axis=1.25, cex.lab=1.25,
     cex=sqrt(FMA$base_n)/2, 
     ylab="Terminal FMA (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")

#Plotting relationships between the baseline and terminal means 
g1<-ggplot(FMA, aes(x = base_m, y = term_m, color = as.factor(outcome_name))) +
    geom_point(size = 4, alpha=0.5) + theme(text = element_text(size=20)) +
    #geom_smooth(method=lm)+
    #scale_color_manual(values=c("#007360","#9b1d36","#369cba")) +
    theme(panel.background=element_rect(fill="white", color="gray"), 
          panel.grid.major=element_line(color="gray"),
          panel.grid.minor=element_line(color="gray")) +
    labs(x = "FMA at Baseline", y = "FMA at Terminal Assessment")
print(g1)



################################################
### Analysis of the gait speed measures only ###
################################################
# only one outcome per study
SPEED<-subset(DATA, speed_outcome == "speed")
# Remove studies missing time for therapy
SPEED<-subset(SPEED,!time_50 == "na") 
# Remove studies missing duration of therapy
SPEED<-subset(SPEED,!exp_dur == "na") 
# Remove studies missing age at baseline
SPEED<-subset(SPEED,!age_base == "na")
# Remove studies missing days post-stroke
SPEED<-subset(SPEED,!days_ps == "na")
# Remove studies missing terminal effect-sizes
SPEED<-subset(SPEED,!term_d == "na") 
length(SPEED$author)
# Length for the SPEED data should be 116
sum(SPEED$base_n)
# Sum of baseline n for SPEED data should be 2,609

#Descriptive Statistics for the SPEED data
summary(SPEED$age_base)
sd(SPEED$age_base, na.rm=TRUE)

summary(SPEED$days_ps)
sd(SPEED$days_ps, na.rm=TRUE)

summary(SPEED$exp_dur)
sd(SPEED$exp_dur, na.rm=TRUE)

summary(SPEED$time_50)
sd(SPEED$time_50, na.rm=TRUE)

summary(SPEED$group)

summary(SPEED$base_n)
sd(SPEED$base_n, na.rm=TRUE)

summary(SPEED$term_g)
sd(SPEED$term_g, na.rm=TRUE)

# Statistical models for the speed data
s0<-rma(term_g, term_Vg, data=SPEED, method="ML")
s0

# Moderator Analysis of Group
summary(SPEED$group)
SPEED$exp.c<-as.numeric(SPEED$group)-1.5
head(SPEED)
s1<-rma(term_g~exp.c, term_Vg, data=SPEED, method="ML")
s1

# Moderator Analysis of Groups and DPS
summary(SPEED$days_ps)
SPEED$DPS.c<-SPEED$days_ps-mean(SPEED$days_ps)
hist(SPEED$DPS.c)

SPEED$rtDPS<-sqrt(SPEED$days_ps)
summary(SPEED$rtDPS)
SPEED$rtDPS.c<-SPEED$rtDPS-mean(SPEED$rtDPS)
hist(SPEED$rtDPS.c)

s2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=SPEED, method="ML")
s2

# Moderator Analysis of Group and Estimated Time Scheduled
summary(SPEED$time_50)
hist(SPEED$time_50)
SPEED$time50.c<-SPEED$time_50-mean(SPEED$time_50)

SPEED$rtTIME<-sqrt(SPEED$time_50)
summary(SPEED$rtTIME)
SPEED$rtTIME.c<-SPEED$rtTIME-mean(SPEED$rtTIME)
hist(SPEED$rtTIME.c)

s3<-rma(term_g~exp.c+rtDPS.c+rtTIME.c, term_Vg, data=SPEED, method="ML")
s3
plot(s3)

# Adding the interaction of DPS and TIME
s4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=SPEED, method="ML")
s4
plot(s4)

# Moderator Analysis of Group, DPS, and Age
summary(SPEED$age_base)
hist(SPEED$age_base)

SPEED$age.c<-SPEED$age_base-mean(SPEED$age_base)

s5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=SPEED, method="ML")
s5
confint(s5)

# Controlling for the effect of therapy duration
summary(SPEED$exp_dur)
hist(SPEED$exp_dur)
SPEED$DUR.c<-SPEED$exp_dur-mean(SPEED$exp_dur)

SPEED$rtDUR<-sqrt(SPEED$exp_dur)
summary(SPEED$rtDUR)
hist(SPEED$rtDUR)
SPEED$rtDUR.c<-SPEED$rtDUR-mean(SPEED$rtDUR)
# Checking for colinearity between TIME and DURATION
plot(SPEED$rtDUR.c,SPEED$rtTIME.c)
cor.test(SPEED$rtDUR.c,SPEED$rtTIME.c)

s6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=SPEED, method="ML")
s6

# Comparisons between the different SPEED models
summary(s0)
# logLik  deviance       AIC       BIC      AICc  
# -80.4021  265.5184  164.8041  170.3113  164.9103  
summary(s1)
# logLik  deviance       AIC       BIC      AICc  
# -76.6772  258.0687  159.3544  167.6151  159.5686
summary(s2)
# logLik  deviance       AIC       BIC      AICc  
# -52.6510  210.0163  113.3019  124.3163  113.6623
summary(s3)
# logLik  deviance       AIC       BIC      AICc  
# -49.8661  204.4464  109.7321  123.5001  110.2776
summary(s4)
# logLik  deviance       AIC       BIC      AICc  
# -48.0225  200.7593  108.0450  124.5665  108.8156 
summary(s5)
# logLik  deviance       AIC       BIC      AICc  
# -44.9879  194.6902  103.9758  123.2510  105.0129
summary(s6)
# logLik  deviance       AIC       BIC      AICc  
# -43.8895  192.4934  103.7791  125.8078  105.1249 

# So s5 is our best fitting model based on the AIC and AICc
# but lets double-check our statistical assumptions/diagnostics
# The distribution of residuals for this model is approximately normal
plot(density(resid(s5))) 
# The distribution of residuals for this model show some heteroscedasticity
plot(resid(s5)~fitted(s5)) 
# but not a significant relationship between fitted values and residuals.
cor.test(resid(s5),fitted(s5))


#Exporting the SPEED data to a .csv file
#write.csv(SPEED, file="data SCOAR SPEED ANALYSIS.csv")

# Plots of the SPEED data
palette(c("purple","chartreuse","cyan"))

plot(SPEED$term_g~SPEED$rtTIME, bty='n', type='p', pch=21, bg=SPEED$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex.axis=1.25, cex.lab=1.25,
     cex=1/sqrt(SPEED$term_Vg)/(1/sqrt(max(SPEED$term_Vg))), 
     ylab="Terminal Gait Speed (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")

# There is an issue with plotting baseline means against terminal means
# First the Fuzaro study need to be removed because it used %Time and not m/s
# Second all of the studies that used cm/s need to be converted to m/s
# This has already been done in the SCOAR SPEED ANALYSIS.csv file on github
# ... so download it from github or recreate the above steps before graphing!
SPEED_Graph<-read.table("data SCOAR SPEED ANALYSIS.csv", header = TRUE, sep=",") 
head(SPEED_Graph)

#Plotting relationships 
g1<-ggplot(SPEED_Graph, aes(x = base_mCON, y = term_mCON, color = as.factor(outcome_name))) +
    geom_point(size = 4, alpha=0.5) + theme(text = element_text(size=20)) +
    xlim(0,2)+ylim(0,2)+
    #geom_smooth(method=lm)+
    #scale_color_manual(values=c("#007360","#9b1d36","#369cba")) +
    theme(panel.background=element_rect(fill="white", color="gray"), 
          panel.grid.major=element_line(color="gray"),
          panel.grid.minor=element_line(color="gray")) +
    labs(x = "Speed at Baseline (m/s)", y = "Speed at Terminal Assessment (m/s)")
print(g1)

