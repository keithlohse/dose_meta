library("metafor"); library("car"); library("dplyr"); library("ggplot2")

#Read in the full data set
DATA<-read.table("data SCOAR TEXT OUTLIERS REMOVED.txt", header = TRUE, sep="\t") 
head(DATA)

#Alternatively, you can read the text file into R directly from my GitHub repo:
#DATA<-read.table("https://raw.github.com/keithlohse/dose_meta/master/data%20SCOAR%20TEXT%20OUTLIERS%20REMOVED.txt",
#    header = TRUE, sep="\t")
#head(DATA)

#Note that all effect sizes greater than d = 3.0 have been removed

########################################
###Analysis of the primary extraction###

#In this first analysis, we took only the primary outcome or
#the first usable secondary outcome if no primary was listed.
#Below, we will conduct separate analyses for the FMA and gait
#speed (regardless of if those measures were primary outcomes).

#Create different sets of data
#1. DATA <- contains all data for experimental and control groups
#2. LOHSE <- contains only primary outcomes and excludes missing cases
#3. BIGN <- contains only primary outcomes for groups with base n > 30
#3. CTRLS <- contains only the control group data from LOHSE
#4. EXPS <- contains on the experimental group data from LOHSE
#5. FMA <- contains data for Fugl-Meyer Assessment outcomes excludes missing cases
#6. SPEED <- contains data for gait speed outcomes excludes missing cases

LOHSE<-subset(DATA, SCOAR_outcome == "lohse")
length(LOHSE$SCOAR_outcome)
sum(LOHSE$base_n)

LOHSE<-subset(LOHSE,!time_50 == "na") #Remove studies missing time for therapy
LOHSE<-subset(LOHSE,!exp_dur == "na") #Remove studies missing duration of therapy
LOHSE<-subset(LOHSE,!age_base == "na") #Remove studies missing age at baseline
LOHSE<-subset(LOHSE,!days_ps == "na") #Remove studies missing days post-stroke
LOHSE<-subset(LOHSE,!term_d == "na") #Remove studies missing terminal effect-sizes
length(LOHSE$SCOAR_outcome)
sum(LOHSE$base_n)

#Creating a subset of only those studies with large sample sizes
BIGN<-subset(LOHSE, base_n >= 20)
length(BIGN$SCOAR_outcome)
sum(BIGN$base_n)

#Creating a subset of control groups only
CTRLS<-subset(LOHSE, group == "ctrl")
summary(CTRLS$group)

#Creating a subset of experimental groups only
EXPS<-subset(LOHSE, group == "exp")
summary(EXPS$group)

#FMA and SPEED datasets are defined below

###############################################
#Descriptive Statistics of the included Groups#
###############################################
head(LOHSE)
summary(LOHSE$age_base) #average age of patients at baseline
sd(LOHSE$age_base)

summary(LOHSE$days_ps) #average days from stroke to the beginning of therapy
sd(LOHSE$days_ps)

summary(LOHSE$exp_dur) #duration of the intervention
sd(LOHSE$exp_dur)

summary(LOHSE$time_MAX) #time scheduled for therapy based on the Max Time calculation
sd(LOHSE$time_MAX)

summary(LOHSE$time_50) #time scheduled for therapy based on the 50% time calculation
sd(LOHSE$time_50)

summary(LOHSE$time_MIN) #time scheduled for therapy based on the Min time calculation
sd(LOHSE$time_MIN)

summary(as.factor(LOHSE$detailed_time)) #How was time in therapy quantified? 1 = time scheduled; 4 = repetitions
summary(as.factor(LOHSE$itt_analysis)) #Did the authors specify an intention to treat analysis

summary(LOHSE$base_n) #Number of subjects contributing to baseline means
sd(LOHSE$base_n)

summary(LOHSE$term_n) #Number of subjects contributing to the terminal means
sd(LOHSE$term_n)

summary(LOHSE$fu_n) #Number of subjects contributing for follow up means
sd(LOHSE$fu_n, na.rm=TRUE) #Note that we did not filter the data by follow up outcomes,
#so there are missing cases that we need to exclude in this calculation

summary(LOHSE$fu_dur) #time from the baseline assessment to the follow-up assessment
sd(LOHSE$fu_dur, na.rm=TRUE)

summary(LOHSE$group) #experimental versus control groups

summary(LOHSE$term_g) #Number of subjects contributing to the terminal means
sd(LOHSE$term_g)


#############################
#Analysis of All Group Data #
#############################

#Overall analysis of all trials
m0<-rma(term_g, term_Vg, data=LOHSE, method="ML")
m0
confint(m0)

m0A<-rma(term_g, term_Vg, data=CTRLS, method="ML") #intecept only model for control groups
m0A

m0B<-rma(term_g, term_Vg, data=EXPS, method="ML") #intercept only model for experimental groups
m0B


#Creating a forest plot to show the RE model of all of the data
forest(m0, slab=paste(LOHSE$author, LOHSE$year, sep=", "), cex=1.5)

#Creating a funnel plot to show potential bias in the full dataset
palette(c("orange","dodgerblue"))
funnel(m0, pch=21, bg=LOHSE$group, cex = 1.5, cex.axis=1.25, cex.lab=1.25, xlab="Terminal Outcome (g)")
#Statistical test of symmetry
regtest(m0, model = "lm") #Applies the test of asymmetry from Eggers (1997)
#Significant result indicates bias in the distribution of results

#Moderator Analysis of Treatment Versus Control Groups
summary(as.numeric(LOHSE$group))
LOHSE$exp.c<-as.numeric(LOHSE$group)-1.5
head(LOHSE)
m1<-rma(term_g~exp.c, term_Vg, data=LOHSE, method="ML")
m1
confint(m1)

#Moderator Analysis of Groups and DPS
summary(LOHSE$days_ps)
LOHSE$DPS.c<-LOHSE$days_ps-mean(LOHSE$days_ps)
hist(LOHSE$days_ps)

LOHSE$rtDPS<-sqrt(LOHSE$days_ps)
summary(LOHSE$rtDPS)
LOHSE$rtDPS.c<-LOHSE$rtDPS-mean(LOHSE$rtDPS)
hist(LOHSE$rtDPS.c)
mean(LOHSE$rtDPS.c)

m2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=LOHSE, method="ML")
m2
confint(m2)

#Moderator Analysis of Group and Estimated Time Scheduled
summary(LOHSE$time_50)
hist(LOHSE$time_50)
LOHSE$time50.c<-LOHSE$time_50-mean(LOHSE$time_50)

LOHSE$rtTIME<-sqrt(LOHSE$time_50)
summary(LOHSE$rtTIME)
LOHSE$rtTIME.c<-LOHSE$rtTIME-mean(LOHSE$rtTIME)
hist(LOHSE$rtTIME)

m3<-rma(term_g~exp.c+rtDPS.c+age.c+rtTIME.c, term_Vg, data=LOHSE, method="ML")
m3
confint(m3)
plot(m3)

#INTERACTION OF DOSE WITH TIME
m4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=LOHSE, method="ML")
m4
confint(m4)
plot(fitted(m4),resid(m4))

#Moderator Analysis of Group, DPS, TIME and Age
summary(LOHSE$age_base)
hist(LOHSE$age_base)

LOHSE$age.c<-LOHSE$age_base-mean(LOHSE$age_base)

m5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=LOHSE, method="ML")
m5
confint(m5)


#Controlling for the effect of therapy duration
summary(LOHSE$exp_dur)
hist(LOHSE$exp_dur)
LOHSE$DUR.c<-LOHSE$exp_dur-mean(LOHSE$exp_dur)

LOHSE$rtDUR<-sqrt(LOHSE$exp_dur)
summary(LOHSE$rtDUR)
hist(LOHSE$rtDUR)
LOHSE$rtDUR.c<-LOHSE$rtDUR-mean(LOHSE$rtDUR)

plot(LOHSE$rtDUR.c,LOHSE$rtTIME.c)
cor.test(LOHSE$rtDUR.c,LOHSE$rtTIME.c)

###
#Adding controlling for the effects of duration
m6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=LOHSE, method="ML")
m6
confint(m6)

summary(LOHSE$rtDPS.c)
summary(LOHSE$rtTIME.c)
quantile(LOHSE$rtTIME.c, c(.10,.25,.50,.75,.90))
summary(LOHSE$age.c)

write.csv(LOHSE, file="data SCOAR R MAIN ANALYSIS.csv")

#Comparison between the different models
summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)


plot(density(resid(m5))) #The distribution of residuals for this model is approximately normal
plot(resid(m5)~fitted(m5)) #The distribution of residuals for this model show some heteroscedasticity
cor.test(resid(m5),fitted(m5)) #but not a significant relationship between fitted values and residuals.


#######################################
#Plots for Multivariable Relationships#
#######################################
palette(c("purple","chartreuse","cyan"))

#Plot of outcomes by DOSE and DPS
plot(LOHSE$term_g~LOHSE$rtTIME, bty='n', type='p', pch=21, bg=LOHSE$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex=sqrt(LOHSE$base_n)/2, cex.axis=1.25, cex.lab=1.25,
     ylab="Terminal Outcome (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")
##Cyan = <90 days post stroke, green = <365 days post stroke, purple = > 365 days post stroke


#Plot of outcomes by AGE
LOHSE$AGEq <- ntile(LOHSE$age_base, 4) #the ntile function chops the variable into X equally sized groups. 
palette(topo.colors(4)) #Topocolors is just a color pallette that I like

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



#############################################
### Analysis of the "Big" groups only, n > 20

#Descriptive Statistics for the BIGN data
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

#Overall analysis of all trials
b0<-rma(term_g, term_Vg, data=BIGN, method="ML")
b0
confint(b0)

#Creating a forest plot to show the RE model of all of the data
forest(b0, slab=paste(BIGN$author, BIGN$year, sep=", "), cex=1.5)

#Creating a funnel plot to show potential bias in the full dataset
palette(c("orange","dodgerblue"))
funnel(b0, pch=21, bg=BIGN$group, cex = 1.5, cex.axis=1.25, cex.lab=1.25, xlab="Terminal Outcome (g)")
#Statistical test of symmetry
regtest(m0, model = "lm") #Applies the test of asymmetry from Eggers (1997)
#Significant result indicates bias in the distribution of results

#Moderator Analysis of Treatment Versus Control Groups
summary(as.numeric(BIGN$group))
BIGN$exp.c<-as.numeric(BIGN$group)-1.5
head(BIGN)
b1<-rma(term_g~exp.c, term_Vg, data=BIGN, method="ML")
b1
confint(b1)

#Moderator Analysis of Groups and DPS
summary(BIGN$days_ps)
BIGN$DPS.c<-BIGN$days_ps-mean(BIGN$days_ps)
hist(BIGN$days_ps)

BIGN$rtDPS<-sqrt(BIGN$days_ps)
summary(BIGN$rtDPS)
BIGN$rtDPS.c<-BIGN$rtDPS-mean(BIGN$rtDPS)
hist(BIGN$rtDPS.c)
mean(BIGN$rtDPS.c)

b2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=BIGN, method="ML")
b2
confint(b2)

#Moderator Analysis of Group and Estimated Time Scheduled
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

#Adding the interaction of DPS and TIME
#Adding interactions with therapy type
b4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=BIGN, method="ML")
b4
confint(b4)
plot(fitted(b4),resid(b4))
plot(density(resid(b4)))

#Moderator Analysis of Group, DPS, TIME and Age
summary(BIGN$age_base)
hist(BIGN$age_base)

BIGN$age.c<-BIGN$age_base-mean(BIGN$age_base)

b5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=BIGN, method="ML")
b5
confint(b5)


#Controlling for the effect of therapy duration
summary(BIGN$exp_dur)
hist(BIGN$exp_dur)
BIGN$DUR.c<-BIGN$exp_dur-mean(BIGN$exp_dur)

BIGN$rtDUR<-sqrt(BIGN$exp_dur)
summary(BIGN$rtDUR)
BIGN$rtDUR.c<-BIGN$rtDUR-mean(BIGN$rtDUR)
hist(BIGN$rtDUR.c)

plot(BIGN$rtDUR.c,BIGN$rtTIME.c)
cor.test(BIGN$rtDUR.c,BIGN$rtTIME.c)

###
#Controlling for the duration of therapy
b6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=BIGN, method="ML")
b6
confint(b6)

summary(b0)
summary(b1)
summary(b2)
summary(b3)
summary(b4)
summary(b5)
summary(b6)

write.csv(BIGN, file="data SCOAR BIG N ANALYSIS.csv")

#Plot of outcomes by AGE
BIGN$AGEq <- ntile(BIGN$age_base, 4) #the ntile function chops the variable into X equally sized groups. 
palette(topo.colors(4)) #Topocolors is just a color pallette that I like

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




#########################################
###Analysis of the FMA for the UE only###

#All outcomes are fma, fma-ue, fma-se, or fma-wh
#fma-le had been excluded
#only one outcome per study
FMA<-subset(DATA, fma_outcome == "fma")
FMA<-subset(FMA,!time_50 == "na") #Remove studies missing time for therapy
FMA<-subset(FMA,!exp_dur == "na") #Remove studies missing duration of therapy
FMA<-subset(FMA,!age_base == "na") #Remove studies missing age at baseline
FMA<-subset(FMA,!days_ps == "na") #Remove studies missing days post-stroke
FMA<-subset(FMA,!term_d == "na") #Remove studies missing terminal effect-sizes
length(FMA$author)
summary(FMA$term_d)
sum(FMA$base_n)

#Descriptive Statistics for the FMA data
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

summary(SPEED$term_g)
sd(SPEED$term_g, na.rm=TRUE)

#Regression models for FMA
f0<-rma(term_g, term_Vg, data=FMA, method="ML")
f0

#Moderator Analysis of Group
summary(FMA$group)
FMA$exp.c<-as.numeric(FMA$group)-1.5
head(FMA)
f1<-rma(term_g~exp.c, term_Vg, data=FMA, method="ML")
f1

#Moderator Analysis of Groups and DPS
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

#Moderator Analysis of Group and Estimated Time Scheduled
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

#Adding the interaction of DPS and TIME
f4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=FMA, method="ML")
f4

plot(fitted(f4),resid(f4))

#Moderator Analysis of Group, DPS, and Age
summary(FMA$age_base)
hist(FMA$age_base)

FMA$age.c<-FMA$age_base-mean(FMA$age_base)

f5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=FMA, method="ML")
f5
confint(f5)

#Controlling for the effect of therapy duration
summary(FMA$exp_dur)
hist(FMA$exp_dur)
FMA$DUR.c<-FMA$exp_dur-mean(FMA$exp_dur)

FMA$rtDUR<-sqrt(FMA$exp_dur)
summary(FMA$rtDUR)
hist(FMA$rtDUR)
FMA$rtDUR.c<-FMA$rtDUR-mean(FMA$rtDUR)

plot(FMA$rtDUR.c,FMA$rtTIME.c)
cor.test(FMA$rtDUR.c,FMA$rtTIME.c)

f6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=FMA, method="ML")
f6

#Comparisons between the different FMA models
summary(f0)
summary(f1)
summary(f2)
summary(f3)
summary(f4)
summary(f5)
summary(f6)

#Exporting the FMA data to a .csv file
write.csv(FMA, file="data SCOAR FMA ANALYSIS.csv")

###Plots of the FMA data###
palette(c("purple","chartreuse","cyan"))

plot(FMA$term_g~FMA$rtTIME, bty='n', type='p', pch=21, bg=FMA$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex.axis=1.25, cex.lab=1.25,
     cex=sqrt(FMA$base_n)/2, 
     ylab="Terminal FMA (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")

#Plotting relationships 
g1<-ggplot(FMA, aes(x = base_m, y = term_m, color = as.factor(outcome_name))) +
    geom_point(size = 4, alpha=0.5) + theme(text = element_text(size=20)) +
    #geom_smooth(method=lm)+
    #scale_color_manual(values=c("#007360","#9b1d36","#369cba")) +
    theme(panel.background=element_rect(fill="white", color="gray"), 
          panel.grid.major=element_line(color="gray"),
          panel.grid.minor=element_line(color="gray")) +
    labs(x = "FMA at Baseline", y = "FMA at Terminal Assessment")
print(g1)



#################################################
#Analysis of the gait speed measures only
#All outcomes are fma, fma-ue, fma-se, or fma-wh
#fma-le had been excluded
#only one outcome per study
SPEED<-subset(DATA, speed_outcome == "speed")
SPEED<-subset(SPEED,!time_50 == "na") #Remove studies missing time for therapy
SPEED<-subset(SPEED,!exp_dur == "na") #Remove studies missing duration of therapy
SPEED<-subset(SPEED,!age_base == "na") #Remove studies missing age at baseline
SPEED<-subset(SPEED,!days_ps == "na") #Remove studies missing days post-stroke
SPEED<-subset(SPEED,!term_d == "na") #Remove studies missing terminal effect-sizes
summary(SPEED$term_d)
length(SPEED$author)
sum(SPEED$base_n)

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

#Statistical models for the speed data
s0<-rma(term_g, term_Vg, data=SPEED, method="ML")
s0

#Moderator Analysis of Group
summary(SPEED$group)
SPEED$exp.c<-as.numeric(SPEED$group)-1.5
head(SPEED)
s1<-rma(term_g~exp.c, term_Vg, data=SPEED, method="ML")
s1

#Moderator Analysis of Groups and DPS
summary(SPEED$days_ps)
SPEED$DPS.c<-SPEED$days_ps-mean(SPEED$days_ps)
hist(SPEED$DPS.c)

SPEED$rtDPS<-sqrt(SPEED$days_ps)
summary(SPEED$rtDPS)
SPEED$rtDPS.c<-SPEED$rtDPS-mean(SPEED$rtDPS)
hist(SPEED$rtDPS.c)
mean(FMA$rtDPS.c)

s2<-rma(term_g~exp.c+rtDPS.c, term_Vg, data=SPEED, method="ML")
s2

#Moderator Analysis of Group and Estimated Time Scheduled
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

#Adding the interaction of DPS and TIME
s4<-rma(term_g~exp.c+rtDPS.c*rtTIME.c, term_Vg, data=SPEED, method="ML")
s4
plot(fitted(s4),resid(s4))

#Moderator Analysis of Group, DPS, and Age
summary(SPEED$age_base)
hist(SPEED$age_base)

SPEED$age.c<-SPEED$age_base-mean(SPEED$age_base)

s5<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c, term_Vg, data=SPEED, method="ML")
s5
confint(s5)

#Controlling for the effect of therapy duration
summary(SPEED$exp_dur)
hist(SPEED$exp_dur)
SPEED$DUR.c<-SPEED$exp_dur-mean(SPEED$exp_dur)

SPEED$rtDUR<-sqrt(SPEED$exp_dur)
summary(SPEED$rtDUR)
hist(SPEED$rtDUR)
SPEED$rtDUR.c<-SPEED$rtDUR-mean(SPEED$rtDUR)

plot(SPEED$rtDUR.c,SPEED$rtTIME.c)
cor.test(SPEED$rtDUR.c,SPEED$rtTIME.c)

s6<-rma(term_g~exp.c+rtDPS.c*rtTIME.c+age.c+rtDUR.c, term_Vg, data=SPEED, method="ML")
s6

#Comparisons between the different SPEED models
summary(s0)
summary(s1)
summary(s2)
summary(s3)
summary(s4)
summary(s5)
summary(s6)

#Exporting the SPEED data to a .csv file
#write.csv(SPEED, file="data SCOAR SPEED ANALYSIS.csv")

###Plots of the SPEED data###
palette(c("purple","chartreuse","cyan"))

plot(SPEED$term_g~SPEED$rtTIME, bty='n', type='p', pch=21, bg=SPEED$days_cat, col="black",  
     lwd=1.0, xlim=c(-1,15), ylim=c(-0.5,3.5), cex.axis=1.25, cex.lab=1.25,
     cex=1/sqrt(SPEED$term_Vg)/(1/sqrt(max(SPEED$term_Vg))), 
     ylab="Terminal Gait Speed (g)", xlab="Time Scheduled for Therapy (sqrt(hrs))")

#There is an issue with plotting baseline means against terminal means
#First the Fuzaro study need to be removed because it used %Time and not m/s
#Second all of the studies that used cm/s need to be converted to m/s
#This has already been done in the SCOAR SPEED ANALYSIS.csv file
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



