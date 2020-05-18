##Diagnostic tools for checking the AFT and PH asssumptions
library(smoothSurv)
####Cox-Snell Residual
cov<-covname[3:16]
##Weibull AFT
weibull_fit<-survreg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                     +sex+age_c+obstruct+perfor+adhere
                     +differ_mod+differ_poor
                     +extent_muscle+extent_serosa+extent_cs
                     +surg+node4,dist="weibull",data=colon_flex)

shape.weibull<-1/weibull_fit$scale
scale0.weibull<-exp(weibull_fit$coef[1])
scale.weibull<-scale0.weibull/exp(as.matrix(colon_flex[,cov])%*%(-weibull_fit$coef[cov]))

weibull.res<-(colon_flex$time.yr/c(scale.weibull))^shape.weibull

#Exponential AFT
exp_fit<-survreg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                 +sex+age_c+obstruct+perfor+adhere
                 +differ_mod+differ_poor
                 +extent_muscle+extent_serosa+extent_cs
                 +surg+node4,dist="exponential",data=colon_flex)
exp.rate0<-exp(-exp_fit$coef[1])
exp.rate1<-exp.rate0*exp(as.matrix(colon_flex[,cov])%*%(-exp_fit$coef[cov])) 
exp.res<-c(exp.rate1)*colon_flex$time.yr


#Log-logistic AFT
llogistic_fit<-survreg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                       +sex+age_c+obstruct+perfor+adhere
                       +differ_mod+differ_poor
                       +extent_muscle+extent_serosa+extent_cs
                       +surg+node4,dist="loglogistic",data=colon_flex)
llogistic.shape<-1/llogistic_fit$scale
llogistic.scale0<-exp(llogistic_fit$coef[1])
llogistic.scale<-llogistic.scale0/exp(as.matrix(colon_flex[,cov])%*%(-llogistic_fit$coef[cov]))
llogistic.S<-1/(1+(colon_flex$time.yr/c(llogistic.scale))^llogistic.shape)
llogistic.res<--log(llogistic.S)

#Log-normal AFT
lnormal_fit<-survreg(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
                     +sex+age_c+obstruct+perfor+adhere
                     +differ_mod+differ_poor
                     +extent_muscle+extent_serosa+extent_cs
                     +surg+node4,dist="lognormal",data=colon_flex)

lnormal.meanlog0<-lnormal_fit$coef[1]
lnormal.sdlog<-lnormal_fit$scale
lnormal.meanlog<-lnormal.meanlog0-as.matrix(colon_flex[,cov])%*%(-lnormal_fit$coef[cov])
lnormal.S<-pnorm(log(colon_flex$time.yr),mean=lnormal.meanlog,sd=lnormal.sdlog,lower.tail=FALSE)
lnormal.res<--log(lnormal.S)

##Spline-based AFT
load("colon_splineaft.rda")
source("SplineAFT.R")
flex.S<-rep(NA,dim(colon_flex)[1])

for (i in 1:dim(colon_flex)[1]){
  cov.val<-as.numeric(colon_flex[i,cov])
  timeobs<-colon_flex[i,"time.yr"]
  flex.S[i]<-SurvEst(fit=flex_colon,time=timeobs,cov=cov.val,Data=colon_flex,plot=FALSE)$survival
  
}
flex.res<--log(flex.S)



###
load("colon_smooth.rda")
smooth.S<-rep(NA,dim(colon_flex)[1])
for (i in 1:dim(colon_flex)[1]){
  cov.val<-as.numeric(colon_flex[i,cov])
  timeobs<-colon_flex[i,"time.yr"]
  sm_t<-survfit(smooth_colon,cov.val,by=0.001,plot=FALSE)$x
  sm_Sest<-survfit(smooth_colon,cov.val,by=0.001,plot=FALSE)$y1
  smooth.S[i]<-sm_Sest[which(abs(sm_t-timeobs)<=0.001)[1]]
}
smooth.res<--log(smooth.S)

cox_fit<-coxph(Surv(time.yr,status)~rx_Lev+rx_Lev5FU
               +sex+age_c+obstruct+perfor+adhere
               +differ_mod+differ_poor
               +extent_muscle+extent_serosa+extent_cs
               +surg+node4,data=colon_flex)

cox.res <- colon_flex$status-resid(cox_fit, type="martingale")

jpeg("coxsnell.jpeg",width = 3000, height = 3000,res=250)

par(mfrow=c(3,3))
par(oma=c(0,0,2,0))
plot(1, type="n", axes=F, xlab="", ylab="")


km.res<-summary(survfit(Surv(flex.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Spline-based AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)
plot(1, type="n", axes=F, xlab="", ylab="")

km.res<-summary(survfit(Surv(smooth.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Smooth Error AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)

km.res<-summary(survfit(Surv(weibull.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Weibull AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)

km.res<-summary(survfit(Surv(exp.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Exponential AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)

km.res<-summary(survfit(Surv(lnormal.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Log-normal AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)

km.res<-summary(survfit(Surv(llogistic.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Log-logistic AFT Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)

km.res<-summary(survfit(Surv(cox.res,colon_flex$status)~1))
plot(km.res$time,km.res$cumhaz,ylab="Cumulative Hazard",xlab="Cox-Snell Residual",main="Cox PH Model",xlim=c(0,2.5),ylim=c(0,2.5))
lines(c(0,2.5),c(0,2.5),lty=2,col="grey",lwd=2)
dev.off()

##Q-Q plot for AFT model
title<-c("Levamisole",
         "Levamisole +5FU",
         "Male",
         "Age",
         "Obstruction",
         "Perforation",
         "Adhere",
         "Moderate differentiation",
         "Poor differentiation",
         "Invasion to muscle",
         "Invasion to serosa",
         "Invasion to contiguous structures",
         "Time since enrollment",
         "Node>4")
jpeg("qqplot.jpeg",width = 3000, height = 3000,res=250)

par(mfrow=c(4,4),mar=c(4,4,2,2))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-summary(survfit(formula, data = colon_flex))
  
  km.fit.surv1<-c(1,km.fit$surv[km.fit$strata==paste0(cov[k],"=1")])
  km.fit.surv0<-c(1,km.fit$surv[km.fit$strata==paste0(cov[k],"=0")])
  km.fit.time1<-c(0,km.fit$time[km.fit$strata==paste0(cov[k],"=1")])
  km.fit.time0<-c(0,km.fit$time[km.fit$strata==paste0(cov[k],"=0")])
  ss<-seq(0.1,0.99,0.1)
  qt1<-rep(NA,length(ss))
  qt0<-rep(NA,length(ss))
  for (i in 1:length(ss)){
    
    qt1[i]<-tail(km.fit.time1[km.fit.surv1>=ss[i]],1)
    qt0[i]<-tail(km.fit.time0[km.fit.surv0>=ss[i]],1)
  }
  
  plot(qt0,qt1,type="l",xlab="T0",ylab="T1",main=title[k],xlim=c(1,7),ylim=c(1,7),cex.main=0.75)
  
}
dev.off()

jpeg("coxph.jpeg",width = 3000, height = 3000,res=250)
##log-log plot for Cox model
par(mfrow=c(4,4),mar=c(4,4,2,2))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-survfit(formula, data = colon_flex)
  
  plot(km.fit, fun = "cloglog", xlab = "Time",
       ylab = "log(-log(S))", main =title[k],cex.main=0.75,xlim=c(0.01,10),lty=2:1)
}
dev.off()

####schoenfield residual
cox_log<-cox.zph(cox_fit, transform ='log')
cox_identity<-cox.zph(cox_fit, transform ='identity')
cox_km<-cox.zph(cox_fit, transform ='km')
cox_rank<-cox.zph(cox_fit, transform ='rank')
jpeg("Schoenfeldres_log.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in 1:14){
  plot(cox_log[k],ann=F,se=F)
  title(xlab="Log Time",ylab="Schoenfeld Residuals",main=title[k])
  abline(h=0, col="red")
}
mtext("Scatterplot of the scaled Schoenfeld residuals based on log time",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

jpeg("Schoenfeldres_identity.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in 1:14){
  plot(cox_identity[k],ann=F,se=F)
  title(xlab="Time",ylab="Schoenfeld Residuals",main=title[k])
  abline(h=0, col="red")
}
mtext("Scatterplot of the scaled Schoenfeld residuals based on time",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

jpeg("Schoenfeldres_km.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in 1:14){
  plot(cox_km[k],ann=F,se=F)
  title(xlab="Kaplan−Meier Estimate",ylab="Schoenfeld Residuals",main=title[k])
  abline(h=0, col="red")
}
mtext("Scatterplot of the scaled Schoenfeld residuals based on KM estimate",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

jpeg("Schoenfeldres_rank.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in 1:14){
  plot(cox_rank[k],ann=F,se=F)
  title(xlab="Rank of Time",ylab="Schoenfeld Residuals",main=title[k])
  abline(h=0, col="red")
}
mtext("Scatterplot of the scaled Schoenfeld residuals based on rank time",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

#exponential  H(t) vs t
jpeg("exp_diag.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-summary(survfit(formula, data = colon_flex))
  
  plot(km.fit$time[km.fit$strata==paste0(cov[k],"=0")],km.fit$cumhaz[km.fit$strata==paste0(cov[k],"=0")],
       xlim=c(0,7),ylim=c(0,1.2),type="s",lty=2,xlab="Time",ylab="Cumulative Hazard",
       main =title[k],cex.main=0.75)
  lines(km.fit$time[km.fit$strata==paste0(cov[k],"=1")],km.fit$cumhaz[km.fit$strata==paste0(cov[k],"=1")],type="s")
}
mtext("Diagnostic plot for exponential distribution",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

##log-log plot (Weibull and Cox)
jpeg("weibull_diag.jpeg",width = 3000, height = 3000,res=250)
logneglog<-function(x) log(-log(x))
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-summary(survfit(formula, data = colon_flex))
  
  plot(log(km.fit$time[km.fit$strata==paste0(cov[k],"=0")]),logneglog(km.fit$surv[km.fit$strata==paste0(cov[k],"=0")]),
       xlim=c(-4,2),ylim=c(-6,1),type="s",lty=2,xlab="log(Time)",ylab="log(-log(S))",
       main=title[k],cex.main=0.75)
  lines(log(km.fit$time[km.fit$strata==paste0(cov[k],"=1")]),logneglog(km.fit$surv[km.fit$strata==paste0(cov[k],"=1")]),type="s")
}
mtext("Diagnostic plot for Weibull distribution",
      side = 3, line = 0.7, outer = TRUE)
dev.off()

#log-logistic  logit(S(t))  vs log(t)
logit<-function(x) log(x/(1-x))
jpeg("llogistic_diag.jpeg",width = 3000, height = 3000,res=250)
par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-summary(survfit(formula, data = colon_flex))
  
  plot(log(km.fit$time[km.fit$strata==paste0(cov[k],"=0")]),logit(km.fit$surv[km.fit$strata==paste0(cov[k],"=0")]),
       xlim=c(-4,2),ylim=c(-2,6),type="s",lty=2,xlab="log(Time)",ylab="logit(S)",
       main=title[k],cex.main=0.75)
  lines(log(km.fit$time[km.fit$strata==paste0(cov[k],"=1")]),logit(km.fit$surv[km.fit$strata==paste0(cov[k],"=1")]),type="s")
}
mtext("Diagnostic plot for log-logistic distribution",
      side = 3, line = 0.7, outer = TRUE)
dev.off()
#log-normal inverse cdf (1-S) vs log(T)

inversecdf<-function(x) qnorm(1-x)
jpeg("lnormal_diag.jpeg",width = 3000, height = 3000,res=250)

par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(1,1,3,1))
for (k in c(1:3,5:14)){
  formula<-as.formula(paste0("Surv(time.yr,status)~",cov[k]))
  km.fit<-summary(survfit(formula, data = colon_flex))
  
  plot(log(km.fit$time[km.fit$strata==paste0(cov[k],"=0")]),inversecdf(km.fit$surv[km.fit$strata==paste0(cov[k],"=0")]),
       xlim=c(-4,2),ylim=c(-6,2),type="s",lty=2,xlab="log(Time)",ylab="logit(S)",
       main=title[k],cex.main=0.75)
  lines(log(km.fit$time[km.fit$strata==paste0(cov[k],"=1")]),inversecdf(km.fit$surv[km.fit$strata==paste0(cov[k],"=1")]),type="s")
}
mtext("Diagnostic plot for log-normal distribution",
      side = 3, line = 0.7, outer = TRUE)
dev.off()