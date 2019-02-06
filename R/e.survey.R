# |> DEFUNCT <| #
make.formula<-function(names) {
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
formula(paste("~",paste(names,collapse="+")))}

# |> DEFUNCT <| #
dimnames.survey.design<-function(x) {
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
dimnames(x$variables)}

# |> DEFUNCT <| #
oldsvydesign<-function(ids,probs=NULL,strata=NULL,variables=NULL, fpc=NULL,
                    data=NULL, nest=FALSE, check.strata=!nest,weights=NULL){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

    .Deprecated("svydesign")
  
    ## less memory-hungry version for sparse tables
    interaction<-function (..., drop = TRUE) {
        args <- list(...)
        narg <- length(args)
        if (narg == 1 && is.list(args[[1]])) {
            args <- args[[1]]
            narg <- length(args)
        }
        
        ls<-sapply(args,function(a) length(levels(a)))
        ans<-do.call("paste",c(lapply(args,as.character),sep="."))
        ans<-factor(ans)
        return(ans)
        
    }


    na.failsafe<-function(object,...){
      if (NCOL(object)==0)
        object
      else na.fail(object)
    }
    
     if(inherits(ids,"formula")) {
     mf<-substitute(model.frame(ids,data=data,na.action=na.failsafe))   
     ids<-eval.parent(mf)
    } else if (!is.null(ids))
            ids<-na.fail(data.frame(ids))

     if(inherits(probs,"formula")){
    mf<-substitute(model.frame(probs,data=data,na.action=na.failsafe))
    probs<-eval.parent(mf)
    }
     
     if(inherits(weights,"formula")){
       mf<-substitute(model.frame(weights,data=data,na.action=na.failsafe))
       weights<-eval.parent(mf)
     } else if (!is.null(weights))
         weights<-na.fail(data.frame(weights))
     
     if(!is.null(weights)){
       if (!is.null(probs))
         stop("Can't specify both sampling weights and probabilities")
       else
         probs<-1/weights
     }

      

    if (!is.null(strata)){
      if(inherits(strata,"formula")){
        mf<-substitute(model.frame(strata,data=data, na.action=na.failsafe))
        strata<-eval.parent(mf)
      }
      if(is.list(strata))
        strata<-na.fail(do.call("interaction", strata))
      if (!is.factor(strata))
        strata<-factor(strata)
      has.strata<-TRUE
    } else {
      strata<-factor(rep(1,NROW(ids)))
      has.strata <-FALSE
    }
    
    if (inherits(variables,"formula")){
        mf<-substitute(model.frame(variables,data=data,na.action=na.pass))
        variables <- eval.parent(mf)
    } else if (is.null(variables)){
        variables<-data
    } else
        variables<-data.frame(variables)

    
     if (inherits(fpc,"formula")){
       mf<-substitute(model.frame(fpc,data=data,na.action=na.failsafe))
       fpc<-eval.parent(mf)
       if (length(fpc))
         fpc<-fpc[,1]
     }
     
    if (is.null(ids) || NCOL(ids)==0)
    ids<-data.frame(.id=seq(length=NROW(variables)))

     ## force subclusters nested in clusters
     if (nest && NCOL(ids)>1){
      N<-ncol(ids)
      for(i in 2:(N)){
          ids[,i]<-do.call("interaction", ids[,1:i,drop=TRUE])
      }
    }
     ## force clusters nested in strata
     if (nest && has.strata && NCOL(ids)){
       N<-NCOL(ids)
       for(i in 1:N)
         ids[,i]<-do.call("interaction", list(strata, ids[,i]))
     }

    ## check if clusters nested in strata 
     if (check.strata && nest)
      warning("No point in check.strata=TRUE if nest=TRUE")
    if(check.strata && !is.null(strata) && NCOL(ids)){
       sc<-rowSums(table(ids[,1],strata)>0)
       if(any(sc>1)) stop("Clusters not nested in strata")
    }

    ## Put degrees of freedom (# of PSUs in each stratum) in object, to 
    ## allow subpopulations
    if (NCOL(ids)){
        nPSU<-table(strata[!duplicated(ids[,1])])
    }


    if (!is.null(fpc)){

       if (NCOL(ids)>1){
         if (all(fpc<1))
           warning("FPC is not currently supported for multi-stage sampling")
         else
           stop("Can't compute FPC from population size for multi-stage sampling")
       }
       
       ## Finite population correction: specified per observation
       if (is.numeric(fpc) && length(fpc)==NROW(variables)){
         tbl<-by(fpc,list(strata),unique)
         if (any(sapply(tbl,length)!=1))
           stop("fpc not constant within strata")
         fpc<-data.frame(strata=factor(rownames(tbl),levels=levels(strata)),
                         N=as.vector(tbl))
       }
       ## Now reduced to fpc per stratum
       nstr<-table(strata[!duplicated(ids[[1]])])
       
       if (all(fpc[,2]<=1)){
         fpc[,2]<- nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
       } else if (any(fpc[,2]<nstr[match(as.character(fpc[,1]), names(nstr))]))
         stop("Over 100% sampling in some strata")
       
     }

    ## if FPC specified, but no weights, use it for weights
    if (is.null(probs) && is.null(weights) && !is.null(fpc)){
      pstr<-nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
      probs<-pstr[match(as.character(strata),as.character(fpc[,1]))]
      probs<-as.vector(probs)
    }


    certainty<-rep(FALSE,length(unique(strata)))
    names(certainty)<-as.character(unique(strata))
    if (any(nPSU==1)){
      ## lonely PSUs: are they certainty PSUs?
      if (!is.null(fpc)){
        certainty<- fpc$N < 1.01
        names(certainty)<-as.character(fpc$strata)
      } else if (all(as.vector(probs)<=1)){
        certainty<- !is.na(match(as.character(unique(strata)),as.character(strata)[probs > 0.99]))
        names(certainty)<-as.character(unique(strata))
      } else {
        warning("Some strata have only one PSU and I can't tell if they are certainty PSUs")
      }
      
    } 
    
    if (is.numeric(probs) && length(probs)==1)
        probs<-rep(probs, NROW(variables))
    
    if (length(probs)==0) probs<-rep(1,NROW(variables))
    
    if (NCOL(probs)==1) probs<-data.frame(probs)

    rval<-list(cluster=ids)
    rval$strata<-strata
    rval$has.strata<-has.strata
    rval$prob<- apply(probs,1,prod) 
    rval$allprob<-probs
    rval$call<-match.call()
    rval$variables<-variables
    rval$fpc<-fpc
    rval$certainty<-certainty
    rval$call<-sys.call()
    rval$nPSU<-nPSU
    class(rval)<-"survey.design"
    rval
  }

# |> DEFUNCT <| #
print.survey.design<-function(x,varnames=FALSE,design.summaries=FALSE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  # .svycheck(x)
  n<-NROW(x$cluster)
  if (x$has.strata) cat("Stratified ")
  un<-length(unique(x$cluster[,1]))
  if(n==un){
    cat("Independent Sampling design\n")
    is.independent<-TRUE
  } else {
    cat(NCOL(x$cluster),"- level Cluster Sampling design\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
    is.independent<-FALSE
  }
  print(x$call)
  if (design.summaries){
    cat("Probabilities:\n")
    print(summary(x$prob))
    if(x$has.strata){
      cat("Stratum sizes: \n")
      a<-rbind(obs=table(x$strata),
           design.PSU=x$nPSU,
               actual.PSU=if(!is.independent || !is.null(x$fpc))
               table(x$strata[!duplicated(x$cluster[,1])]))
      print(a)
    }
    if (!is.null(x$fpc)){
      if (x$has.strata) {
        cat("Population stratum sizes (PSUs): \n")
        print(x$fpc)
      } else {
        cat("Population size (PSUs):",x$fpc[,2],"\n")
      }
    }
  }
  if (varnames){
    cat("Data variables:\n")
    print(names(x$variables))
  }
  invisible(x)
}

# |> DEFUNCT <| #
"[.survey.design"<-function (x,i, ...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  
  if (!missing(i)){ 
    if (is.calibrated(x)){
      tmp<-x$prob[i,]
      x$prob<-rep(Inf, length(x$prob))
      x$prob[i,]<-tmp
    } else {
      x$variables<-"[.data.frame"(x$variables,i,...,drop=FALSE)
      x$cluster<-x$cluster[i,,drop=FALSE]
      x$prob<-x$prob[i]
      x$allprob<-x$allprob[i,,drop=FALSE]
      x$strata<-x$strata[i]
    }
  } else {
    x$variables<-x$variables[,...,drop=FALSE]
  }
  
  x
}

# |> DEFUNCT <| #
"[<-.survey.design"<-function(x, ...,value){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  if (inherits(value, "survey.design"))
    value<-value$variables
  x$variables[...]<-value
  x
}

# |> ALIVE <| #
dim.survey.design<-function(x,...){
    dim(x$variables)
}

# |> ALIVE <| #
model.frame.survey.design<-function(formula,...){
  formula$variables
}

# |> DEFUNCT <| #
na.fail.survey.design<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    tmp<-na.fail(object$variables,...)
    object
}

# |> DEFUNCT <| #
na.omit.survey.design<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  tmp<-na.omit(object$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object<-object[-omit,]
    object$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

# |> DEFUNCT <| #
na.exclude.survey.design<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    tmp<-na.exclude(object$variables,...)
    exclude<-attr(tmp,"na.action")
    if (length(exclude)){
           object<-object[-exclude,]
       object$variables<-tmp
       attr(object,"na.action")<-exclude
    }
    object
}

# |> DEFUNCT <| #
update.survey.design<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$variables[,newnames[j]]<-eval(dots[[j]],object$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}

# |> DEFUNCT <| #
subset.survey.design<-function(x,subset,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
        e <- substitute(subset)
        r <- eval(e, x$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
    x$call<-sys.call(-1)
    x
}

# |> DEFUNCT <| #
summary.survey.design<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  class(object)<-"summary.survey.design"
  object
}

# |> DEFUNCT <| #
print.summary.survey.design<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  y<-x
  class(y)<-"survey.design"
  print(y,varnames=TRUE,design.summaries=TRUE,...)
}    

# |> DEFUNCT <| #
postStratify.survey.design<-function(design, strata, population, partial=FALSE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  if(inherits(strata,"formula")){
    mf<-substitute(model.frame(strata, data=design$variables))
    strata<-eval.parent(mf)
  }
  strata<-as.data.frame(strata)

  sampletable<-xtabs(I(1/design$prob)~.,data=strata)
  sampletable<-as.data.frame(sampletable)

  if (inherits(population,"table"))
    population<-as.data.frame(population)
  else if (!is.data.frame(population))
    stop("population must be a table or dataframe")

  if (!all(names(strata) %in% names(population)))
    stop("Stratifying variables don't match")
  nn<- names(population) %in% names(strata)
  if (sum(!nn)!=1)
    stop("stratifying variables don't match")

  names(population)[which(!nn)]<-"Pop.Freq"
  
  both<-merge(sampletable, population, by=names(strata), all=TRUE)

  samplezero <- both$Freq %in% c(0, NA)
  popzero <- both$Pop.Freq %in% c(0, NA)
  both<-both[!(samplezero & popzero),]
  
  if (any(onlysample<- popzero & !samplezero)){
    print(both[onlysample,])
    stop("Strata in sample absent from population. This Can't Happen")
  }
  if (any(onlypop <- samplezero & !popzero)){
    if (partial){
      both<-both[!onlypop,]
      warning("Some strata absent from sample: ignored")
    } else {
      print(both[onlypop,])
      stop("Some strata absent from sample: use partial=TRUE to ignore them.")
    }
  } 

  reweight<-both$Pop.Freq/both$Freq
  both$label <- do.call("interaction", list(both[,names(strata)]))
  designlabel <- do.call("interaction", strata)
  index<-match(designlabel, both$label)

  attr(index,"oldweights")<-1/design$prob
  design$prob<-design$prob/reweight[index]
  attr(index,"weights")<-1/design$prob
  design$postStrata<-c(design$postStrata,list(index))
  
  ## Do we need to iterate here a la raking to get design strata
  ## and post-strata both balanced?
  design$call<-sys.call(-1)
  
  design
}

# |> DEFUNCT <| #
svyCprod<-function(x, strata, psu, fpc, nPSU, certainty=NULL, postStrata=NULL,
                   lonely.psu=getOption("RG.lonely.psu")){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  x<-as.matrix(x)
  n<-NROW(x)

  ## Remove post-stratum means, which may cut across PSUs
  if(!is.null(postStrata)){
    for (psvar in postStrata){
      if (inherits(psvar, "greg_calibration") || inherits(psvar, "raking"))
        stop("rake() and calibrate() not supported for old-style design objects")
      psw<-attr(psvar,"weights")
      psmeans<-rowsum(x/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
      x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
    }
  }

  ##First collapse over PSUs

  if (is.null(strata)){
    strata<-rep("1",n)
    if (!is.null(nPSU))
        names(nPSU)<-"1"
  }
  else
    strata<-as.character(strata) ##can't use factors as indices in for()'

  if (is.null(certainty)){
    certainty<-rep(FALSE,length(strata))
    names(certainty)<-strata
  }
  
  if (!is.null(psu)){
    x<-rowsum(x, psu, reorder=FALSE)
    strata<-strata[!duplicated(psu)]
    n<-NROW(x)
  }
  
  if (!is.null(nPSU)){
      obsn<-table(strata)
      dropped<-nPSU[match(names(obsn),names(nPSU))]-obsn
      if(sum(dropped)){
        xtra<-matrix(0,ncol=NCOL(x),nrow=sum(dropped))
        strata<-c(strata,rep(names(dropped),dropped))
          if(is.matrix(x))
       x<-rbind(x,xtra)
        else
       x<-c(x,xtra)
        n<-NROW(x)
      }
  } else obsn<-table(strata)

  if(is.null(strata)){
      x<-t(t(x)-colMeans(x))
  } else {
      strata.means<-drop(rowsum(x,strata, reorder=FALSE))/drop(rowsum(rep(1,n),strata, reorder=FALSE))
      if (!is.matrix(strata.means))
          strata.means<-matrix(strata.means, ncol=NCOL(x))
      x<- x- strata.means[ match(strata, unique(strata)),,drop=FALSE]
  }
  
  p<-NCOL(x)
  v<-matrix(0,p,p)
  
  ss<-unique(strata)
  for(s in ss){
      this.stratum <- strata %in% s
      
      ## original number of PSUs in this stratum 
      ## before missing data/subsetting
      this.n <-nPSU[match(s,names(nPSU))]
      
      this.df <- this.n/(this.n-1)    
      
      if (is.null(fpc))
          this.fpc <- 1
      else{
          this.fpc <- fpc[,2][ fpc[,1]==as.character(s)]
          this.fpc <- (this.fpc - this.n)/this.fpc
      }
      
      xs<-x[this.stratum,,drop=FALSE]

      this.certain<-certainty[names(certainty) %in% s]
      
      ## stratum with only 1 design cluster leads to undefined variance
      lonely.psu<-match.arg(lonely.psu, c("remove","adjust","fail",
                                          "certainty","average"))
      if (this.n==1 && !this.certain){
        this.df<-1
        if (lonely.psu=="fail")
          stop("Stratum ",s, " has only one sampling unit.")
        else if (lonely.psu!="certainty")
          warning("Stratum ",s, " has only one sampling unit.")
        if (lonely.psu=="adjust")
          xs<-strata.means[match(s,ss),,drop=FALSE]
      } else if (obsn[match(s,names(obsn))]==1 && !this.certain){
        ## stratum with only 1 cluster left after subsetting 
        warning("Stratum ",s," has only one PSU in this subset.")
        if (lonely.psu=="adjust")
          xs<-strata.means[match(s,ss),,drop=FALSE]
      }
      ## add it up
      if (!this.certain)
        v<-v+crossprod(xs)*this.df*this.fpc
    }
  if (lonely.psu=="average"){
    v<- v/(1-mean(obsn==1 & !certainty))
  }
  v
}



svymean<-function(x, design,na.rm=FALSE,...){
# .svycheck(design)
  UseMethod("svymean",design)
}

# |> DEFUNCT <| #
svymean.survey.design<-function(x,design, na.rm=FALSE,deff=FALSE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  }
  else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$variables)
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
            design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-svyCprod(x*pweights/psum,design$strata,design$cluster[[1]], design$fpc,
              design$nPSU,design$certainty, design$postStrata)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
    nobs<-NROW(design$cluster)
    vsrs<-svyvar(x,design,na.rm=na.rm)/nobs
    vsrs<-vsrs*(psum-nobs)/psum
    attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}

# |> DEFUNCT <| #
print.svystat<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    vv<-attr(x,"var")
    if (is.matrix(vv))
        m<-cbind(x,sqrt(diag(vv)))
    else
        m<-cbind(x,sqrt(vv))
    hasdeff<-!is.null(attr(x,"deff"))
    if (hasdeff) {
        m<-cbind(m,deff(x))
        colnames(m)<-c(attr(x,"statistic"),"SE","DEff")
    } else {
        colnames(m)<-c(attr(x,"statistic"),"SE")
    }
    printCoefmat(m)
}

# |> DEFUNCT <| #
as.data.frame.svystat<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  rval<-data.frame(statistic=coef(x),SE=SE(x))
  names(rval)[1]<-attr(x,"statistic")
  if (!is.null(attr(x,"deff")))
    rval<-cbind(rval,deff=deff(x))
  rval
}

# |> ALIVE <| #
coef.svystat<-function(object,...){
  attr(object,"statistic")<-NULL
  attr(object,"deff")<-NULL
  attr(object,"var")<-NULL
  unclass(object)
}

# |> ALIVE <| #
vcov.svystat<-function(object,...){
  as.matrix(attr(object,"var"))
}

deff <- function(object, ...) UseMethod("deff")

# |> ALIVE <| #
deff.default <- function(object, ...){
  rval<-attr(object,"deff")
  if (is.null(rval)) { 
      warning("object has no design effect information")
  } else rval<-diag(as.matrix(rval))
  rval
}


SE <- function(object,...){
  UseMethod("SE")
}

# |> ALIVE <| #
SE.default<-function(object,...){
  sqrt(diag(vcov(object,...)))
}

VAR <- function(object,...){
  UseMethod("VAR")
}

# |> DEFUNCT <| #
VAR.default<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  SE(object,...)^2
}

cv<-function(object,...) UseMethod("cv")

# |> ALIVE <| #
cv.default<-function(object, warn=TRUE, ...){
  rval<-SE(object)/coef(object)
  if (warn && any(coef(object)<0,na.rm=TRUE)) warning("CV may not be useful for negative statistics")
  rval
}

#--------------------------------- STOPPED HERE 14/004/2016 -------------------#

svytotal<-function(x,design,na.rm=FALSE,...){
# .svycheck(design)
  UseMethod("svytotal",design)
}

# |> DEFUNCT <| #
svytotal.survey.design<-function(x,design, na.rm=FALSE, deff=FALSE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
    if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  } else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, design$variables)

  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }

  N<-sum(1/design$prob)
  m <- svymean(x, design, na.rm=na.rm)
  total<-m*N
  attr(total, "var")<-v<-svyCprod(x/design$prob,design$strata,
                                  design$cluster[[1]], design$fpc,
                                  design$nPSU,design$certainty,design$postStrata)
  attr(total,"statistic")<-"total"
  if (is.character(deff) || deff){
    vsrs<-svyvar(x,design)*sum(weights(design)^2)
    vsrs<-vsrs*(N-NROW(design$cluster))/N
    attr(total,"deff")<-v/vsrs
  }
  return(total)
}

svyvar<-function(x, design, na.rm=FALSE,...){
# .svycheck(design)
  UseMethod("svyvar",design)
}

# |> DEFUNCT <| #
svyvar.survey.design<-function(x, design, na.rm=FALSE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    
    if (inherits(x,"formula"))
            x<-model.frame(x,design$variables,na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        
    n<-NROW(x)
    xbar<-svymean(x,design, na.rm=na.rm)
    if(NCOL(x)==1) {
            x<-x-xbar
            v<-svymean(x*x*n/(n-1),design, na.rm=na.rm)
            attr(v,"statistic")<-"variance"
            return(v)
    }
    x<-t(t(x)-xbar)
    p<-NCOL(x)
    a<-matrix(rep(x,p),ncol=p*p)
    b<-x[,rep(1:p,each=p)]
        ## Kish uses the n-1 divisor, so it affects design effects
    v<-svymean(a*b*n/(n-1),design, na.rm=na.rm)
    v<-matrix(v,ncol=p)
        attr(v,"statistic")<-"variance"
        v
    }



svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)

# Method not actually working on ReGenesees design objects (see .analytic method below)
# |> DEFUNCT <| #
svyquantile.survey.design<-function(x,design,quantiles,alpha=0.05,
                                    ci=FALSE, method="linear",f=1,
                                    interval.type=c("Wald","score","betaWald"),
                                    na.rm=FALSE,se=ci, ties=c("discrete","rounded"), ...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    if (inherits(x,"formula"))
      x<-model.frame(x ,model.frame(design), na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, model.frame(design,na.action=na.pass))
    
    if (na.rm){
        nas<-rowSums(is.na(x))
        design<-design[nas==0,]
        if (length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
        else
          x[nas>0,]<-0
      }

    w<-weights(design)
    
    computeQuantiles<-function(xx,p=quantiles){
      if (any(is.na(x))) return(NA*p)
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method=method,f=f,
                     yleft=min(xx),yright=max(xx)) 
      cdf(p)
    }
    
    computeQuantilesRounded<-function(xx,p=quantiles){
      if (any(is.na(xx))) return(NA*p)
      ww<-rowsum(w,xx,reorder=TRUE)
      xx<-sort(unique(xx))
      cum.w <- cumsum(ww)/sum(ww)
      cdf <- approxfun(cum.w, xx, method = method, f = f, 
                       yleft = min(xx), yright = max(xx))
      cdf(p)
    }
      
    
    
    computeScoreCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      
      U<-function(theta){ ((xx>theta)-(1-p))}
        
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }
      
      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }
    
    computePCI<-function(se,alpha,p){
      if (interval.type=="Wald"){
        p.up<-p+qnorm(alpha/2,lower.tail=FALSE)*se
        p.low<-p+qnorm(alpha/2,lower.tail=TRUE)*se
        c(p.low,p.up)
      } else if (interval.type=="betaWald"){
        n.eff <- (p*(1-p))/(se^2)
        n.eff <- n.eff * ( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
        p.up<-qbeta(1-alpha/2, n.eff*p+1, n.eff*(1-p))
        p.low<-qbeta(alpha/2,  n.eff*p, n.eff*(1-p)+1)
        c(p.low,p.up)
      }
      
    }
    
    computeWaldCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      theta0<-computeQuantiles(xx,p)
      U<- ((xx>theta0)-(1-p))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      approx(cum.w,xx[oo],xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx))$y 
      
    }
    
    computeWaldCIRounded<-function(xx,p){
        theta0<-computeQuantilesRounded(xx,p)
        U<- ((xx>theta0)-(1-p))
        ww<-rowsum(w,xx, reorder=TRUE)
        uxx <- sort(unique(xx))
        wtest<-svymean(U,design)
        p.ci<-computePCI(SE(wtest),alpha,p)
        p.low<-p.ci[1]
        p.up<-p.ci[2]
        oo<-order(xx)
        cum.w<-cumsum(ww)/sum(ww)
        approx(cum.w,uxx,xout=c(p.low,p.up), method=method,f=f,
               yleft=min(xx),yright=max(xx))$y 
        
      }

    ties<-match.arg(ties)
    computeQ<-switch(ties, discrete=computeQuantiles,rounded=computeQuantilesRounded)
    
    if (!is.null(dim(x)))
        rval<-t(matrix(apply(x,2,computeQ),nrow=length(quantiles),
                       dimnames=list(as.character(round(quantiles,2)),colnames(x))))
    else
      rval<-computeQ(x)
    
    if (!ci & !se) return(rval)
    
    interval.type<-match.arg(interval.type)
    
    computeCI<-switch(paste(interval.type,ties,sep="."), score.discrete=computeScoreCI,
                            score.rounded=stop("ties=\"rounded\" not available with interval.type=\"score\""),
                            Wald.rounded=computeWaldCIRounded,
                            betaWald.rounded=computeWaldCIRounded,
                            Wald.discrete=computeWaldCI,
                            betaWald.discrete=computeWaldCI)
    
    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
                   as.character(round(quantiles,2)),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    if (ci)
      rval<-list(quantiles=rval,CIs=cis)
    else
      rval<-list(quantiles=rval)
    
    if (is.null(dim(x)))
        ses<-(cis[2,]-cis[1,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    else
        ses<-(cis[2,,]-cis[1,,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    attr(rval,"SE")<-ses
    class(rval)<-"svyquantile"
    rval
  }
  
svyquantile.analytic <- function(x,design,quantiles,alpha=0.05,
                                 ci=FALSE, method="linear",f=1,
                                 interval.type=c("Wald","score","betaWald"),
                                 na.rm=FALSE,se=ci, ties=c("discrete","rounded"), ...){
##########################################
# svyquantile method for class analytic. #
##########################################
    if (inherits(x,"formula"))
      x<-model.frame(x ,model.frame(design), na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, model.frame(design,na.action=na.pass))
    
    ## Missing values treatment
    nas<-rowSums(is.na(x))
    if (na.rm && sum(nas)>0){
        # If domain has some non-NA values, use them for estimation:
        if (length(design[nas==0,]$prob) > 0) {
            design<-design[nas==0,]
            if (length(nas)>length(design$prob))
                x<-x[nas==0,,drop=FALSE]
            else
                x[nas>0,]<-0
        }
        # If domain has only NAs, cannot do anything (i.e. behave as na.rm=FALSE)
        else {
            na.rm <- FALSE
        }
    }

    ww <- weights(design)
    nw <- length(ww)
    # Discard observations with weights less than or equal to zero (if any)
    w <- ww[ww > 0]
    nw.pos <- length(w)
    if ( nw.pos < nw ){
        design <- design[ww > 0,]
        x <- x[ww > 0,,drop=FALSE]
        warning("Observations with weights less than or equal to zero (", nw - nw.pos,") have been discarded.")
       }

    # Check if, after removing non-positive weights and/or nas, there are still observations left
    no.obs <- (nrow(design) <= 0)
    
    computeQuantiles<-function(xx,p=quantiles){
      if (any(is.na(x))) return(NA*p)
      if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(NA*p)
        }
        
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method=method,f=f,
                     yleft=min(xx),yright=max(xx)) 
      cdf(p)
    }
    
    computeQuantilesRounded<-function(xx,p=quantiles){
      if (any(is.na(xx))) return(NA*p)
      if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(NA*p)
        }
      ww<-rowsum(w,xx,reorder=TRUE)
      xx<-sort(unique(xx))
      cum.w <- cumsum(ww)/sum(ww)
      cdf <- approxfun(cum.w, xx, method = method, f = f, 
                       yleft = min(xx), yright = max(xx))
      cdf(p)
    }
    
    computeScoreCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      
      U<-function(theta){ ((xx>theta)-(1-p))}
        
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }
      
      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }
    
    computePCI<-function(se,alpha,p){
      if (interval.type=="Wald"){
        p.up<-p+qnorm(alpha/2,lower.tail=FALSE)*se
        p.low<-p+qnorm(alpha/2,lower.tail=TRUE)*se
        c(p.low,p.up)
      } else if (interval.type=="betaWald"){
        n.eff <- (p*(1-p))/(se^2)
        n.eff <- n.eff * ( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
        p.up<-qbeta(1-alpha/2, n.eff*p+1, n.eff*(1-p))
        p.low<-qbeta(alpha/2,  n.eff*p, n.eff*(1-p)+1)
        c(p.low,p.up)
      }
      
    }
    
    computeWaldCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
        if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(c(NA,NA))
        }
      theta0<-computeQuantiles(xx,p)
      U<- ((xx>theta0)-(1-p))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      approx(cum.w,xx[oo],xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx))$y 
    }
    
    computeWaldCIRounded<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
        if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(c(NA,NA))
        }
      theta0<-computeQuantilesRounded(xx,p)
      U<- ((xx>theta0)-(1-p))
      ww<-rowsum(w,xx, reorder=TRUE)
      uxx <- sort(unique(xx))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(ww)/sum(ww)
      approx(cum.w,uxx,xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx))$y 
    }

    ties<-match.arg(ties)
    computeQ<-switch(ties, discrete=computeQuantiles,rounded=computeQuantilesRounded)
    
    if (!is.null(dim(x)))
        rval<-t(matrix(apply(x,2,computeQ),nrow=length(quantiles),
#                       dimnames=list(as.character(round(quantiles,2)),colnames(x))))
                       dimnames=list(format(round(quantiles, 3), digits = 3,nsmall=3), colnames(x))))
    else
      rval<-computeQ(x)
    
    if (!ci & !se) return(rval)
    
    interval.type<-match.arg(interval.type)
    
    computeCI<-switch(paste(interval.type,ties,sep="."), score.discrete=computeScoreCI,
                            score.rounded=stop("ties=\"rounded\" not available with interval.type=\"score\""),
                            Wald.rounded=computeWaldCIRounded,
                            betaWald.rounded=computeWaldCIRounded,
                            Wald.discrete=computeWaldCI,
                            betaWald.discrete=computeWaldCI)
    
    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
#                   as.character(round(quantiles,2)),
                   format(round(quantiles, 3), digits = 3,nsmall=3),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    if (ci)
      rval<-list(quantiles=rval,CIs=cis)
    else
      rval<-list(quantiles=rval)
    
    if (is.null(dim(x)))
        ses<-(cis[2,]-cis[1,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    else
        ses<-(cis[2,,]-cis[1,,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    attr(rval,"SE")<-ses
    class(rval)<-"svyquantile"
    rval
  }


svyquantile.cal.analytic <- function(x,design,quantiles,alpha=0.05,
                                 ci=FALSE, method="linear",f=1,
                                 interval.type=c("Wald","score","betaWald"),
                                 na.rm=FALSE,se=ci, ties=c("discrete","rounded"), ...){
##########################################################
# svyquantile method for class cal.analytic.             #
# Extension needed due to subsetting differences between #
# classes analytic and cal.analytic.                     #
# For calibrated designs it is not possible to discard   #
# negative weights in the same way achieved for analytic #
# objects.                                               #
##########################################################
    if (inherits(x,"formula"))
      x<-model.frame(x ,model.frame(design), na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, model.frame(design,na.action=na.pass))
    
    ## Missing values treatment
    nas<-rowSums(is.na(x))
    if (na.rm && sum(nas)>0){
        # If domain has some non-NA values, use them for estimation:
        if (length(design[nas==0,]$prob) > 0) {
            design<-design[nas==0,]
            if (length(nas)>length(design$prob))
                x<-x[nas==0,,drop=FALSE]
            else
                x[nas>0,]<-0
        }
        # If domain has only NAs, cannot do anything (i.e. behave as na.rm=FALSE)
        else {
            na.rm <- FALSE
        }
    }

    w <- weights(design)
    nw <- length(w)
    # "Discard" observations with weights less than zero (if any)
    # For calibrated objects "discard" means set the corresponding x value to zero
    w.pos <- w[w >= 0]
    nw.pos <- length(w.pos)
    if ( nw.pos < nw ){
        design <- design[w >= 0,]
        if (length(w)>length(design$prob)) {
              # thus the design IS NOT CALIBRATED
              # since has been actually subsetted
              x <- x[w >= 0,,drop=FALSE]
              stop("This should not happen! Please report this message to the author (zardetto@istat.it)")
            }
        else { 
              x[w < 0,] <- 0
            }
        warning("Observations with weights less than zero (", nw - nw.pos,") have been discarded.")
       }

    # Check if, after removing negative weights and/or nas, there are still observations left:
    # this cannot be assessed by counting rows... must find a different solution (temporarily
    # remove the check)
    no.obs <- FALSE

    computeQuantiles<-function(xx,p=quantiles){
      if (any(is.na(x))) return(NA*p)
      if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(NA*p)
        }
        
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method=method,f=f,
                     yleft=min(xx),yright=max(xx)) 
      cdf(p)
    }
    
    computeQuantilesRounded<-function(xx,p=quantiles){
      if (any(is.na(xx))) return(NA*p)
      if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(NA*p)
        }
      ww<-rowsum(w,xx,reorder=TRUE)
      xx<-sort(unique(xx))
      cum.w <- cumsum(ww)/sum(ww)
      cdf <- approxfun(cum.w, xx, method = method, f = f, 
                       yleft = min(xx), yright = max(xx))
      cdf(p)
    }
    
    computeScoreCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      
      U<-function(theta){ ((xx>theta)-(1-p))}
        
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }
      
      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),
                qlimit=qnorm(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }
    
    computePCI<-function(se,alpha,p){
      if (interval.type=="Wald"){
        p.up<-p+qnorm(alpha/2,lower.tail=FALSE)*se
        p.low<-p+qnorm(alpha/2,lower.tail=TRUE)*se
        c(p.low,p.up)
      } else if (interval.type=="betaWald"){
        n.eff <- (p*(1-p))/(se^2)
        n.eff <- n.eff * ( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
        p.up<-qbeta(1-alpha/2, n.eff*p+1, n.eff*(1-p))
        p.low<-qbeta(alpha/2,  n.eff*p, n.eff*(1-p)+1)
        c(p.low,p.up)
      }
      
    }
    
    computeWaldCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
        if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(c(NA,NA))
        }
      theta0<-computeQuantiles(xx,p)
      U<- ((xx>theta0)-(1-p))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      approx(cum.w,xx[oo],xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx))$y 
    }
    
    computeWaldCIRounded<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
        if ( !(nw.pos > 0) || no.obs ) {
          warning("No positive weights in (sub)population!")
          return(c(NA,NA))
        }
      theta0<-computeQuantilesRounded(xx,p)
      U<- ((xx>theta0)-(1-p))
      ww<-rowsum(w,xx, reorder=TRUE)
      uxx <- sort(unique(xx))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(ww)/sum(ww)
      approx(cum.w,uxx,xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx))$y 
    }

    ties<-match.arg(ties)
    computeQ<-switch(ties, discrete=computeQuantiles,rounded=computeQuantilesRounded)
    
    if (!is.null(dim(x)))
        rval<-t(matrix(apply(x,2,computeQ),nrow=length(quantiles),
#                       dimnames=list(as.character(round(quantiles,2)),colnames(x))))
                       dimnames=list(format(round(quantiles, 3), digits = 3,nsmall=3), colnames(x))))
    else
      rval<-computeQ(x)
    
    if (!ci & !se) return(rval)
    
    interval.type<-match.arg(interval.type)
    
    computeCI<-switch(paste(interval.type,ties,sep="."), score.discrete=computeScoreCI,
                            score.rounded=stop("ties=\"rounded\" not available with interval.type=\"score\""),
                            Wald.rounded=computeWaldCIRounded,
                            betaWald.rounded=computeWaldCIRounded,
                            Wald.discrete=computeWaldCI,
                            betaWald.discrete=computeWaldCI)
    
    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
#                   as.character(round(quantiles,2)),
                   format(round(quantiles, 3), digits = 3,nsmall=3),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    if (ci)
      rval<-list(quantiles=rval,CIs=cis)
    else
      rval<-list(quantiles=rval)
    
    if (is.null(dim(x)))
        ses<-(cis[2,]-cis[1,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    else
        ses<-(cis[2,,]-cis[1,,])/(2*qnorm(alpha/2,lower.tail=FALSE))
    attr(rval,"SE")<-ses
    class(rval)<-"svyquantile"
    rval
  }




# |> ALIVE <| #  
SE.svyquantile<-function(object,...){
    attr(object,"SE")
}

# |> DEFUNCT <| #
vcov.svyquantile<-function(object,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  se<-SE(object)
  if (is.null(se)) stop("no uncertainty information present")
  v<-matrix(NA,length(se),length(se))
  warning("Only diagonal of vcov() available")
  diag(v)<-se
  v
}

# |> ALIVE <| #
coef.svyquantile<-function(object,...){
  rval<-as.vector(object$quantiles)
  if(ncol(object$quantiles)==1)
    names(rval)<-rownames(object$quantiles)
  else if (nrow(object$quantiles)==1)
    names(rval)<-colnames(object$quantiles)
  else names(rval)<-t(outer(colnames(object$quantiles),
                            rownames(object$quantiles),
                            paste,sep=":"))
  rval
}

# |> DEFUNCT <| #
print.svyquantile<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
    print(list(quantiles=x$quantiles, CIs=x$CIs))
}


degf<-function(design,...) UseMethod("degf")

# |> DEFUNCT <| #
degf.survey.design2<-function(design,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  inset<- weights(design,"sampling")!=0
  length(unique(design$cluster[inset, 1])) - length(unique(design$strata[inset, 1]))
}

# |> ALIVE <| #
coef.svyratio<-function(object,...,drop=TRUE){
  if (!drop) return(object$ratio)
  cf<-as.vector(object$ratio)
  cross <- attr(object,"cross")
  if (cross) {
      nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    }
  else {
      nms<-names(object$ratio)
    }
  names(cf)<-nms
  cf
}

# |> ALIVE <| #
deff.svyratio<-function(object,...,drop=TRUE){
  if (!drop) return(attr(object,"deff"))
  deff <-as.vector(attr(object,"deff"))
  cross <- attr(object,"cross")
  if (cross) {
      nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    }
  else {
      nms<-names(object$ratio)
    }
  names(deff)<-nms
  deff
}

# |> ALIVE <| #
SE.svyratio<-function(object,...,drop=TRUE){
  if(!drop) return(sqrt(object$var))
  se<-as.vector(sqrt(object$var))
  cross <- attr(object,"cross")
  if (cross) {
      nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    }
  else {
      nms<-names(object$var)
    }
  names(se)<-nms
  se
}

# |> ALIVE <| #
cv.svyratio<-function(object,...){
# modified to print correct names
  cv <- as.vector(sqrt(object$var))/as.vector(object$ratio)
  cross <- attr(object,"cross")
  if (cross) {
      nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    }
  else {
      nms<-names(object$var)
    }
  names(cv)<-nms
  cv
}

# |> ALIVE <| #
vcov.svyratio <- function(object, ...){
  covmat<-object$vcov
  if (is.null(covmat)){
      covmat<-matrix(NaN,length(object$var),length(object$var))
      diag(covmat)<-as.vector(object$var)
  }
  cross <- attr(object,"cross")
  if (cross) {
      nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    }
  else {
      nms<-names(object$var)
    }
    dimnames(covmat)<-list(nms,nms)
    covmat
}

svyratio<-function(numerator,denominator, design,...){
# .svycheck(design)
  UseMethod("svyratio",design)
}

# |> DEFUNCT <| #
svyratio.survey.design<-function(numerator, denominator, design,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,design$variables)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,design$variables)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    allstats<-svytotal(all,design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-svyCprod(r*1/design$prob, design$strata, design$cluster[[1]], design$fpc,
                            design$nPSU, design$certainty,design$postStrata)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    rval$call<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }

# |> DEFUNCT <| #
print.svyratio_separate<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  cat("Stratified ratio estimate: ")
  if (!is.null(x$call))
    print(x$call)
  else if (!is.null(attr(x,"call")))
    print(attr(x$call))
  for(r in x$ratios) {
    print(r)
  }
  invisible(x)
}

# |> DEFUNCT <| #
print.svyratio<-function(x,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  cat("Ratio estimator: ")
  if (!is.null(x$call))
    print(x$call)
  else if(!is.null(attr(x,"call")))
    print(attr(x,"call"))
  cat("Ratios=\n")
  print(x$ratio)
  cat("SEs=\n")
  print(sqrt(x$var))
  invisible(NULL)
}

# |> DEFUNCT <| #
predict.svyratio<-function(object, total, se=TRUE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
  if (se)
    return(list(total=object$ratio*total,se=sqrt(object$var)*total))
  else
    return(object$ratio*total)
}

# |> DEFUNCT <| #
predict.svyratio_separate<-function(object, total, se=TRUE,...){
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )

  if (length(total)!=length(object$ratios))
    stop("Number of strata differ in ratio object and totals.")
  if (!is.null(names(total)) && !is.null(levels(object$strata))){
    if (!setequal(names(total), levels(object$strata)))
      warning("Names of strata differ in ratio object and totals")
    else if (!all(names(total)==levels(object$strata))){
      warning("Reordering supplied totals to make their names match the ratio object")
      total<-total[match(names(total),levels(object$strata))]
    }
  }
  totals<-mapply(predict, object=object$ratios, total=total,se=se,...,SIMPLIFY=FALSE)

  if(se){
    rval<-totals[[1]]$total
    v<-totals[[1]]$se^2
    for(ti in totals[-1]) {
      rval<-rval+ti$total
      v<-v+ti$se^2
    }
    list(total=rval,se=sqrt(v))
  } else {
    rval<-totals[[1]]
    for (ti in totals[-1]) rval<-rval+ti
    rval
  }

}
