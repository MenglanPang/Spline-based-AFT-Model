#Example for spline-based AFT model using colon dataset, which is available in 'survival' package

setwd("~/Desktop/OneDrive - McGill University/Flexible AFT Model R Program/Spline-Based AFT")
source("SplineAFT.R")
library(survival)
library(dummies)
library(splines)

colon<-colon[complete.cases(colon) & colon$etype==2, ]
colon$time.yr<-colon$time/365.25

colon$rx_Lev<-dummy(colon$rx)[,2]
colon$rx_Lev5FU<-dummy(colon$rx)[,3]

colon$age_c<-scale(colon$age)

colon$differ_mod<-ifelse(colon$differ==2,1,0)
colon$differ_poor<-ifelse(colon$differ==3,1,0)

colon$extent_muscle<-ifelse(colon$extent==2,1,0)
colon$extent_serosa<-ifelse(colon$extent==3,1,0)
colon$extent_cs<-ifelse(colon$extent==4,1,0)

colon_flex<-cbind(colon$id,colon$time.yr,colon$rx_Lev,colon$rx_Lev5FU,
                  colon$sex,colon$age_c,colon$obstruct,colon$perfor,colon$adhere,
                  colon$differ_mod,colon$differ_poor,
                  colon$extent_muscle,colon$extent_serosa,colon$extent_cs,
                  colon$surg,colon$node4,colon$status,colon$age)

colon_flex<-data.frame(colon_flex)
covname<-c("id","time.yr","rx_Lev","rx_Lev5FU",
           "sex","age_c","obstruct","perfor","adhere",
           "differ_mod","differ_poor",
           "extent_muscle","extent_serosa","extent_cs",
           "surg","node4","status","age")

colnames(colon_flex)<-covname

##spline-based AFT model
##Input:
#Data: Dataset (no missing data;categorial variables have to be coded as dummy variables)
#Var: names of the variable to be included in the spline-based AFT model
#Time.Obs: name of the observed time variable
#Delta: name of the event indicator
#degree.bh: degree of the splines for modeling the baseline hazard
#nknot.bh: number of the interior knots of the splines for modeling the baseline hazard
#tol: convergece criter evaluating two consecutive log-likelihood;default is 1^-5
sink("colon_cancer_splineaft.txt")
flex_colon<-SplineAFT(Data=colon_flex,Var=covname[3:16],Time.Obs="time.yr",Delta="status",degree.bh=3,nknot.bh=2,tol=1e-5)
sink()
save(flex_colon,file="colon_splineaft.rda")

###estimated hazard function for a specific covariate pattern
cov.val<-rep(NA,14)
for (i in 1:14){
  cov.val[i]<-0
}

##Input
#fit: object after fitting the spline-based AFT model using the SplineAFT function
#time: a vector of time at which the hazard function is estimated
#cov: the specific covariate pattern for which hazard function is estimated
#Data: Dataset that used to fit the spline-based AFT model

par(mfrow=c(1,2))
haz.est<-HazardEst(fit=flex_colon,time=seq(0,6,0.1),cov=cov.val,Data=colon_flex,plot=FALSE)
plot(haz.est$time,haz.est$hazard,type="l",col="blue",xlab="Time (Year)",ylab="Hazard")

###estimated survival function for a specific covariate pattern
surv.est<-SurvEst(fit=flex_colon,time=seq(0,6,0.1),cov=cov.val,Data=colon_flex,plot=FALSE)
plot(surv.est$time,surv.est$survival,type="l",col="blue",xlab="Time (Year)",ylab="Survival")


##smoothed error AFT model
library(smoothSurv)
smooth_colon<-smoothSurvReg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                            +sex+age_c+obstruct+perfor+adhere
                            +differ_mod+differ_poor
                            +extent_muscle+extent_serosa+extent_cs
                            +surg+node4,data=colon_flex,info=FALSE)
save(smooth_colon,file=paste0("colon_smooth.rda"))

