
make.calfun<-function(Fm1,dF, name){
  if (!identical(names(formals(Fm1)), c("u","bounds")))
    stop("wrong argument names for Fm1")
  if(!identical(names(formals(dF)), c("u","bounds")))
    stop("wrong argument names for dF")
  rval<-list(Fm1=Fm1, dF=dF, name=name)
  class(rval)<-"calfun"
  rval
}

print.calfun<-function(x,...) cat("calibration metric: ",x$name,"\n")

calibrate<-function(design, ...) UseMethod("calibrate")

calibrate.survey.design2<-function(design, formula, population,
                                    aggregate.stage=NULL, stage=0, variance=NULL,
                                    bounds=c(-Inf,Inf), calfun=c("linear","raking","logit"),
                                    maxit=50, epsilon=1e-7, verbose=FALSE, force=FALSE,
                                    ...){
  if (is.character(calfun)) calfun<-match.arg(calfun)
  if (is.character(calfun) && calfun=="linear" && (bounds==c(-Inf,Inf))){
    ## old code is better for ill-conditioned linear calibration
    rval<-regcalibrate(design,formula,population,
                       aggregate.stage=aggregate.stage, stage=stage,
                       lambda=variance,...)
    rval$call<-sys.call(-1)
    return(rval)
  }

  if(is.character(calfun))
    calfun<-switch(calfun,linear=cal.linear, raking=cal.raking, logit=cal.logit)
  else
    if(!inherits(calfun,"calfun"))
      stop("'calfun' must be a string or of class 'calfun'.")
  
  if (!is.null(aggregate.stage)){
    aggindex<-design$cluster[[aggregate.stage]]
  }

  expit<-function(x) 1-1/(1+exp(x))
  
  ## calibration to population totals
  mm<-model.matrix(formula, model.frame(formula, model.frame(design)))
  ww<-weights(design)
  
  if (!is.null(aggregate.stage)){
    mm<-apply(mm,2,function(mx) ave(mx,aggindex))
    ww<-ave(ww,aggindex)
  }
  whalf<-sqrt(ww)
  sample.total<-colSums(mm*ww)
  
  if(any(sample.total==0)){
    ## drop columsn where all sample and population are zero
    zz<-(population==0) & (apply(mm,2,function(x) all(x==0)))
    mm<-mm[,!zz]
    population<-population[!zz]
    sample.total<-sample.total[!zz]
  }

    
  if (length(sample.total)!=length(population))
    stop("Population and sample totals are not the same length.")

  if(!is.null(names(population))){
    if (!all(names(sample.total) %in% names(population)))
      warning("Sampling and population totals have different names.")
    else if (!all(names(sample.total) == names(population))){
      warning("Sample and population totals reordered to make names agree: check results.")
      population <- population[match(names(sample.total), names(population))]
    }
  }
  
  tqr<-qr(mm*whalf)
  if (!all(abs(qr.resid(tqr,whalf))<1e-10))
    warning("G-calibration models must have an intercept")

  g<-grake(mm,ww,calfun, bounds=bounds,population=population,
           verbose=verbose,epsilon=epsilon,maxit=maxit)

  if (!force && !is.null(attr(g,"failed"))) stop("Calibration failed")
  design$prob<-design$prob/g
  
  caldata <- list(qr=tqr, w=g*whalf, stage=0, index=NULL)
  
  class(caldata) <- c("greg_calibration","gen_raking")
  
  design$postStrata <- c(design$postStrata, list(caldata))
  design$call <- sys.call(-1)
  
  design
}


grake<-function(mm,ww,calfun,eta=rep(0,NCOL(mm)),bounds,population,epsilon, verbose,maxit){

  sample.total<-colSums(mm*ww)
  ## No longer needed: ReGenesees now IMPORTS MASS
  # require(MASS) ##ginv
  if(!inherits(calfun,"calfun")) stop("'calfun' must be of class 'calfun'")
  
  Fm1<-calfun$Fm1
  dF<-calfun$dF

  xeta<-drop(mm%*%eta)
  g<-1+Fm1(xeta, bounds)

  iter<-1

  repeat({
    Tmat<-crossprod(mm*ww*dF(xeta, bounds), mm)

    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta, bounds)))
    deta<-ginv(Tmat, tol=256*.Machine$double.eps)%*%misfit
    eta<-eta+deta

    xeta<- drop(mm%*%eta)
    g<-1+Fm1(xeta, bounds)
    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta, bounds)))
    
    if (verbose)
      print(misfit)

    if (all(abs(misfit)/(1+abs(population))<epsilon)) break

    iter <- iter+1
    if (iter>maxit) {
       achieved<-max((abs(misfit)/(1+abs(population))))
       warning("Failed to converge: eps=",achieved," in ",iter," iterations")
       attr(g,"failed")<-achieved
       break;
     }
  })

  attr(g,"eta")<-eta
  g
}
