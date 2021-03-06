# Example for spline-based AFT model using colon dataset, which is available in 'survival' package #

source("SplineAFT.R")
library(survival)
library(splines)

colon <- colon[complete.cases(colon) & colon$etype==2, ]
colon$time.yr <- colon$time/365.25

colon$rx_Lev <- ifelse(colon$rx=="Lev",1,0)
colon$rx_Lev5FU <- ifelse(colon$rx=="Lev+5FU",1,0)

# To avoid potential numerical problems, it is recommended to standardze/normalize continuous variables using Z-score transformation:
# X* = (X - mean)/SD
# (If preferable, "mean" can be replaced by e.g, Median or other value within the range of the data)
colon$age_c <- scale(colon$age)

colon$differ_mod <- ifelse(colon$differ==2,1,0)
colon$differ_poor <- ifelse(colon$differ==3,1,0)

colon$extent_muscle <- ifelse(colon$extent==2,1,0)
colon$extent_serosa <- ifelse(colon$extent==3,1,0)
colon$extent_cs <- ifelse(colon$extent==4,1,0)

colon_flex <- cbind(colon$id,colon$time.yr,colon$rx_Lev,colon$rx_Lev5FU,
                    colon$sex,colon$age_c,colon$obstruct,colon$perfor,colon$adhere,
                    colon$differ_mod,colon$differ_poor,
                    colon$extent_muscle,colon$extent_serosa,colon$extent_cs,
                    colon$surg,colon$node4,colon$status,colon$age)

colon_flex <- data.frame(colon_flex)
covname <- c("id","time.yr","rx_Lev","rx_Lev5FU",
             "sex","age_c","obstruct","perfor","adhere",
             "differ_mod","differ_poor",
             "extent_muscle","extent_serosa","extent_cs",
             "surg","node4","status","age")

colnames(colon_flex) <- covname

##### Fit the spline-based AFT model #####

## Input for SplineAFT function:
# Data: dataset (no missing data allowed; categorial variables have to be coded as dummy variables)
# Var: names of the independent variables to be included in the spline-based AFT model
# Time.Obs: name of the observed time variable 
# (NOTE: convert the time variable to years or months, instead of days, can help reduce potentially long computation time)
# Delta: name of the event indicator
# degree.bh: degree of the splines for modeling the baseline hazard
# nknot.bh: number of interior knots of the splines for modeling the baseline hazard
# tol: convergence criterion evaluating two consecutive log-likelihoods; default is 1^-5
# ndivision: the number of intervals defining the granularity of the numerical computation of integrals;
# the range of the integral is divided by ndivision small intervals.
# The larger the number is, the longer the estimation takes.
# The default value is 100*floor(max(Time.Obs)) and can be requested by ndivision="maxtime".  
# When deciding a customized value for ndivision for the trade-off between precision and computation time,
# one may take into account the maximum of Time.Obs and how precise the event is being recorded.
# For example, if the maximum of Time.Obs is 3 years and the precision of event times is recorded in months,
# then one may choose ndivision=36, so that the numerical integration uses one month for each interval.
# If the maximum of Time.Obs is 3 years and the precision of event times is recorded in days, 
# then we need a better granularity for the numerical computation.
# One may choose ndivision=300, so that the numerical integration uses 3.65 days for each interval ((3*365)/300).

sink("colon_cancer_splineaft.txt")
flex_colon <- SplineAFT(Data=colon_flex,Var=covname[3:16],Time.Obs="time.yr",Delta="status",degree.bh=3,nknot.bh=2,tol=1e-5,ndivision="maxtime")
sink()
save(flex_colon,file="colon_splineaft.RData")

##### Estimated hazard function for a specific covariate pattern #####

cov.val <- rep(NA,14)
for (i in 1:14){
  cov.val[i] <- 0
}

## Input for HazardEst function:
# fit: object obtained after fitting the spline-based AFT model using the SplineAFT function
# time: a vector of time values at which the hazard function is estimated
# cov: the specific covariate pattern for which hazard function is estimated
# Data: dataset that was used to fit the spline-based AFT model

par(mfrow=c(1,2))
haz.est <- HazardEst(fit=flex_colon,time=seq(0,6,0.1),cov=cov.val,Data=colon_flex,plot=FALSE)
plot(haz.est$time,haz.est$hazard,type="l",col="blue",xlab="Time (Year)",ylab="Hazard")

##### Estimated survival function for a specific covariate pattern #####

## Input for SurvEst function:
# fit: object obtained after fitting the spline-based AFT model using the SplineAFT function
# time: a vector of time values at which the survival curve is estimated
# cov: the specific covariate pattern for which survival curve is estimated
# Data: dataset that was used to fit the spline-based AFT model

surv.est <- SurvEst(fit=flex_colon,time=seq(0,6,0.1),cov=cov.val,Data=colon_flex,plot=FALSE)
plot(surv.est$time,surv.est$survival,type="l",col="blue",xlab="Time (Year)",ylab="Survival")


##### Smoothed error AFT model developed by Komárek et al.#####
# Komárek A, Lesaffre E, Hilton JF. Journal of Computational and Graphical Statistics. 2005;14(3):726-45.

library(smoothSurv)
smooth_colon <- smoothSurvReg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                              +sex+age_c+obstruct+perfor+adhere
                              +differ_mod+differ_poor
                              +extent_muscle+extent_serosa+extent_cs
                              +surg+node4,data=colon_flex,info=FALSE)
save(smooth_colon,file=paste0("colon_smooth.RData"))

