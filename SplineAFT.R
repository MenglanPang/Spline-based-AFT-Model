###
### Program for spline-based AFT model in time-to-event analysis
###

# Author: Menglan Pang
# Last update: May 4, 2018 

# Reference:(under revision at Stat in Med)
# Pang M, Platt RW, Schuster T, Abrahamowicz M. Spline-based Accelerated Failure Time model 
# (2019)

# Functions in this program:
# SplineAFT
# AFT.logLik.gamma.der
# AFT.logLik.beta.der
# AFT.logLik
# integrand
# integrat
# HazardEst
# S.integrand
# SurvEst


library(rootSolve)

get_args_for <- function(fun, env = parent.frame(), inherits = FALSE, ..., dots) {
  potential <- names(formals(fun))
  
  if ("..." %in% potential) {
    if (missing(dots)) {
      # return everything from parent frame
      return(as.list(env))
    }
    else if (!is.list(dots)) {
      stop("If provided, 'dots' should be a list.")
    }
    
    potential <- setdiff(potential, "...")
  }
  
  # get all formal arguments that can be found in parent frame
  args <- mget(potential, env, ..., ifnotfound = list(NULL), inherits = inherits)
  # remove not found
  args <- args[sapply(args, Negate(is.null))]
  # return found args and dots
  c(args, dots)
}

SplineAFT<-function(Data,Var,Time.Obs,Delta,degree.bh,nknot.bh,tol=1e-5,ndivision="maxtime"){
  
  #Covariate matrix from the data
  Cov<-as.matrix(Data[,Var])
  #Observed event time vector
  Tt<-as.matrix(Data[,Time.Obs])
  #event indicator
  delta<-Data[,Delta]
  #number of covariates
  nvar<-length(Var)
  
  #number of events
  num.event<-sum(delta)
  
  # initial values for beta (regression coefficients)
  df.beta<-nvar
  betaVec<-matrix(rep(0,df.beta),ncol=1)
  # initial values for gamma (regression coefficients for the spline basis for approximating the baselind hazard)
  df.gamma<-degree.bh+nknot.bh+1
  gammaVec<-matrix(rep(1,df.gamma),ncol=1)
  #knots position for the splines
  bh.quantile<-seq(1,nknot.bh)/(nknot.bh+1)
  bh.knots<-quantile(Tt,probs=bh.quantile)
  
  loglikelihood.new<-9999
  diff<-1  
  j<-0
  starttime<-proc.time()[3] 
  while(diff>tol){
    j<-j+1
    cat("Iteration=",j,"\n")
    loglikelihood.old<-loglikelihood.new
    
    starttime_gamma<-proc.time()[3] 
    #Estimate gamma conditional on beta
    environment(AFT.logLik.gamma.der) <- environment()
    
    f<-AFT.logLik.gamma.der
    start<-gammaVec
    maxiter = 100;
    rtol = 1e-6; atol = 1e-8; ctol = 1e-8;
    useFortran = TRUE; positive = FALSE;
    jacfunc = NULL; jactype = "fullint";
    verbose = FALSE; bandup = 1; banddown = 1;
    arg_list <- get_args_for(multiroot, dots = list(NULL))
    update_gamma<- do.call(multiroot, arg_list)
    
    #update_gamma<-multiroot(f=AFT.logLik.gamma.der,start=gammaVec, 
    #                        parms=list(betaVec=betaVec,Cov=Cov,
    #                                   Tt=Tt,delta=delta,degree.bh=degree.bh,
    #                                   bh.knots=bh.knots,
    #                                   df.beta=df.beta,
    #                                   df.gamma=df.gamma))
    cat("gamma=",update_gamma$root,"\n")
    cat("time_gamma",proc.time()[3]-starttime_gamma,"\n")
    
    gamma.new<-update_gamma$root
    
    
    gamma.new<-matrix(gamma.new,ncol=1)
    
    starttime_beta<-proc.time()[3] 
    #estimate beta conditional on gamma
    environment(AFT.logLik.beta.der) <- environment()
    f<-AFT.logLik.beta.der
    start<-betaVec
    arg_list <- get_args_for(multiroot, dots = list(NULL))
    update_beta<- do.call(multiroot, arg_list)
    
    #update_beta<-multiroot(f=AFT.logLik.beta.der,start=betaVec, 
    #                       parms=list(gammaVec=gamma.new,Cov=Cov,
    #                                  Tt=Tt,delta=delta,degree.bh=degree.bh,
    #                                  bh.knots=bh.knots,
    #                                  df.beta=df.beta,
    #                                  df.gamma=df.gamma))
    cat("beta=",update_beta$root,"\n")
    cat("time_beta",proc.time()[3]-starttime_beta,"\n")
    
    beta.new<-update_beta$root
    
    
    
    gammaVec<-matrix(gamma.new,ncol=1)
    betaVec<-matrix(beta.new,ncol=1)
    
    environment(AFT.logLik) <- environment()
    #calculate the log-likelihood using the most updated beta and gamma
    loglikelihood.new<-unlist(AFT.logLik(betaVec=betaVec,gammaVec=gamma.new
                                         #,
                                         #Cov=Cov,Tt=Tt,delta=delta,
                                         #degree.bh=degree.bh,bh.knots=bh.knots,
                                         #df.beta=df.beta,df.gamma=df.gamma
                                         ))
    
    cat("loglikelihood=",loglikelihood.new,"\n")
    diff<-abs(loglikelihood.new-loglikelihood.old)
    
  }
  
  elapsedtime<-unname(proc.time()[3]-starttime)
  rownames(betaVec)<-Var
  
  names(bh.knots)<-bh.quantile
  
  df<-df.beta+df.gamma
  
  loglikelihood.new<-unname(loglikelihood.new)
  
  return(list("Var"=Var,"Time.Obs"=Time.Obs,"Delta"=Delta,
              "coefficients"=betaVec,"spline_coef_bh"=gammaVec,
              "degree.bh"=degree.bh,"nknot.bh"=nknot.bh,"bh.knots"=bh.knots,
              "num.events"=num.event,"df"=df,"logLikelihood"=loglikelihood.new,
              "ndivision"=ndivision,
              "runtime"=elapsedtime))
}

#first derivative of the loglikelihood with respect to gamma
AFT.logLik.gamma.der<-function(gammaVec
  #parms
  ){
  #gamma<-gammaVec
  #beta<-parms$betaVec
  #Cov<-parms$Cov
  #Tt<-parms$Tt
  #delta<-parms$delta
  #degree.bh<-parms$degree.bh
  #bh.knots<-parms$bh.knots
  #df.beta<-parms$df.beta
  #df.gamma<-parms$df.gamma
  gamma<-gammaVec
  beta<-betaVec

  W<-exp(Cov%*%beta)*Tt
  W.range<-c(0,W)
  BsW<-bs(W,degree=degree.bh,intercept=TRUE,knots=bh.knots,Boundary.knots=range(W.range))
  environment(integrat) <- environment()
  Interg<-integrat("gamma"
                   #,Cov,beta,gamma,df.beta,df.gamma,
                   #Tt,degree.bh,bh.knots
                   )
  Interg.matrix<-matrix(unlist(Interg),nrow=length(Tt))
  loglik.der<-drop(delta%*%BsW-apply(Interg.matrix,2,sum))
  
  return(LogLik.der=loglik.der)
}

#first derivative of the loglikelihood with respect to beta
AFT.logLik.beta.der<-function(betaVec
                              #,parms
  ){
  beta<-betaVec
  gamma<-gamma.new
  #gamma<-parms$gammaVec
  #Cov<-parms$Cov
  #Tt<-parms$Tt
  #delta<-parms$delta
  #degree.bh<-parms$degree.bh
  #bh.knots<-parms$bh.knots
  #df.beta<-parms$df.beta
  #df.gamma<-parms$df.gamma

  environment(integrat) <- environment()
  Interg<-integrat("beta"
                   #,Cov,beta,gamma,df.beta,df.gamma,
                   #Tt,degree.bh,bh.knots
                   )
  Interg.matrix<-matrix(unlist(Interg),nrow=length(Tt))
  
  loglik.der<-rep(NA,length=df.beta)
  
  for (j in 1:df.beta){
    loglik.der[j]<-sum(delta*Cov[,j]-Interg.matrix[,j])
  }
  
  return(LogLik.der=loglik.der)
}

#loglikelihood calculation
AFT.logLik<-function(betaVec,gammaVec
                     #,Cov,Tt,delta,
                     #degree.bh,bh.knots,
                     #df.beta,df.gamma
                     ){
  beta<-betaVec  
  gamma<-gammaVec

  ### Compute the spline basis for the full dataset, which is need for log likelihood calculation
  #compute W to determine the boundary knots for the splines
  W<-exp(Cov%*%beta)*Tt
  W.range<-c(0,W)
  BsW<-bs(W,degree=degree.bh,intercept=TRUE,knots=bh.knots,Boundary.knots=range(W.range))
  logL1<-delta%*%(Cov%*%beta+BsW%*%gamma)
  
  environment(integrat) <- environment()
  logL2<-integrat("none"
                  #,Cov,beta,gamma,df.beta,df.gamma,
                  #Tt,degree.bh,bh.knots
                  )
  LogLik<-logL1-sum(unlist(logL2))
  
  return(list(LogLik=LogLik))
}

integrand<-function(u,wrt
                    #,Cov,beta,gamma,df.beta,df.gamma,
                    #Tt,degree.bh,bh.knots
                    ){
  n.sam<-length(Tt)
  w<-matrix(NA,nrow=n.sam,ncol=ncol(u))
  exp1<-exp(Cov%*%beta)
  exp2<-matrix(NA,nrow=n.sam,ncol=ncol(u))
  
  W.matrix<-matrix(NA,nrow=n.sam,ncol=ncol(u))
  
  W.matrix<-drop(exp(Cov%*%beta))*u
  if (wrt=="beta") {
    Return.ls<-vector("list",df.beta)
    for (i in 1:df.beta) {Return.ls[[i]]<-matrix(NA,nrow=n.sam,ncol=ncol(u))}
  }
  
  if (wrt=="gamma") {
    Return.ls<-vector("list",df.gamma)
    for (i in 1:df.gamma) {Return.ls[[i]]<-matrix(NA,nrow=n.sam,ncol=ncol(u))}
  }
  
  if (wrt=="none"){
    Return.ls<-vector("list",1)
  }
  
  for (i in 1:ncol(u)){
    Bs.bh<-bs(W.matrix[,i],degree=degree.bh,intercept=TRUE,knots=bh.knots,Boundary.knots=c(0,max(W.matrix)))
    exp2[,i]<-exp(Bs.bh%*%gamma)
    
    if (wrt=="beta"){
      for (k in 1:df.beta){
        Return.ls[[k]][,i]<-exp1*exp2[,i]*Cov[,k]
      }
    } else if (wrt=="gamma"){
      for (k in 1:df.gamma){
        Return.ls[[k]][,i]<-exp1*exp2[,i]*Bs.bh[,k]}
    } 
  }
  if (wrt=="none"){
    Return.ls[[1]]<-drop(exp1)*exp2
  }
  
  return(Return.ls)
}

integrat<-function(wrt
                   #,Cov,beta,gamma,df.beta,df.gamma,
                   #Tt,degree.bh,bh.knots
                   ){
  bound<-cbind(0,Tt)
  #divide the entire range to many small intervals;
  #the number of intervals is to be specified by user (ndivision argument)
  #it defines the granularity of the calculation
  #the larger the number is, the longer the estimation takes
  #we may take into account the the maximum of the event time and data granularity,i.e., how precise the event is being recorded
  #by default, the number of interval=100*floor(maximum of the observed time)
  
  if (ndivision=="maxtime"){
    max_obsT<-floor(max(Tt))  
    num_divide<-max_obsT*100
  } else{
    num_divide<-ndivision
  }
  
  xmatrix<-t(apply(bound,1,function(x) {seq(x[1],x[2],length=num_divide)}))
  step<-apply(xmatrix,1,function(x) (x[2]-x[1]))
  
  xmatrix<-(xmatrix+step/2)[,-ncol(xmatrix)]
  environment(integrand) <- environment()
  yvalue<-integrand(xmatrix,wrt
                    #,Cov,beta,gamma,df.beta,df.gamma,
                    #Tt,degree.bh,bh.knots
                    )
  
  value<-lapply(yvalue,function(x) apply(x*step,1,sum))
  
  return(value)
}




# Estimate of hazard function by spline-based AFT
HazardEst<-function(fit,time,cov,Data,plot=TRUE){
  
  Var<-fit$Var
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  nknot.bh<-fit$nknot.bh
  degree.bh<-fit$degree.bh
  
  #Observed event time vector
  Tt<-as.matrix(Data[,Time.Obs])
  delta<-as.matrix(Data[,Delta])
  Cov<-as.matrix(Data[,Var])
  knots<-quantile(Tt,probs=seq(0,1,by=1/(nknot.bh+1))[-c(1,nknot.bh+2)])
  beta<-fit$coefficients
  gamma<-fit$spline_coef_bh

  
  #Define the boundary knots for the spline basis
  W<-c(0,c(exp(Cov%*%beta))*Tt)
  BsW<-bs(c(exp(cov%*%beta))*time,degree=degree.bh,intercept=TRUE,knots=knots,Boundary.knots=range(W))
  flex.hazard<-c(exp(cov%*%beta))*exp(BsW%*%gamma)
  if (plot==TRUE){
    plot(time,flex.hazard,col="grey",type="l",main="Spline-based AFT",xlim=c(0,ceiling(max(time))),ylab="Hazard",xlab="Time")
    Time.event<-Tt[delta==1 & Tt<=max(time)]
    Time.event.freq<-table(Time.event)
    Time.event.uniq<-as.numeric(names(Time.event.freq))
    for (i in 1:length(Time.event.uniq)){
      rug(Time.event.uniq[i],ticksize=0.03*Time.event.freq[i])
    } 
  }
  return(list("time"=time,"hazard"=flex.hazard))
}

#Estimation of survival function by spline-based AFT
S.integrand<-function(u,cov=cov,degree.bh=degree.bh,knots=knots,W=W,beta=beta,gamma=gamma)
{
  BsW<-bs(exp(cov%*%beta)*u,degree=degree.bh,intercept=TRUE,knots=knots,Boundary.knots=range(W))
  return(exp(cov%*%beta)*exp(BsW%*%gamma))
}

SurvEst<-function(fit,time,cov,Data,plot=TRUE){
  Var<-fit$Var
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  nknot.bh<-fit$nknot.bh
  degree.bh<-fit$degree.bh
  
  #Observed event time vector
  Tt<-as.matrix(Data[,Time.Obs])
  delta<-as.matrix(Data[,Delta])
  Cov<-as.matrix(Data[,Var])
  knots<-quantile(Tt,probs=seq(0,1,by=1/(nknot.bh+1))[-c(1,nknot.bh+2)])
  beta<-fit$coefficients
  gamma<-fit$spline_coef_bh

  
  W<-c(0,c(exp(Cov%*%beta))*Tt)
  
  flex.S<-exp(-1*sapply(time,function(x) integrate(Vectorize(S.integrand,"u"),
                                                   lower=0,upper=x,cov=cov,degree.bh=degree.bh,knots=knots,W=W,beta=beta,gamma=gamma)$value))
  if (plot==TRUE){
    plot(time,flex.S,col="grey",type="l",xlim=c(0,max(time)),ylim=c(0,1),ylab="Survival",xlab="Time",main="Spline-based AFT")
    Time.event<-Tt[delta==1 & Tt<=max(time)]
    Time.event.freq<-table(Time.event)
    Time.event.uniq<-as.numeric(names(Time.event.freq))
    for (i in 1:length(Time.event.uniq)){
      rug(Time.event.uniq[i],ticksize=0.03*Time.event.freq[i])
    }
  }
  
  return(list("time"=time,"survival"=flex.S))
}