##
##  Recursive estimation of linearisation variances
##  in multistage samples.
##

svydesign<-function(ids, probs = NULL, strata = NULL, variables = NULL, 
    fpc = NULL, data=NULL, nest = FALSE, check.strata = !nest, 
    weights = NULL,...){
    UseMethod("svydesign", data)
    }

    
weights.survey.design<-function(object,...){
  return(1/object$prob)
}

    
svydesign.default<-function(ids,probs=NULL,strata=NULL,variables=NULL, fpc=NULL,
                    data=NULL, nest=FALSE, check.strata=!nest,weights=NULL,...){

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

    na.failsafe<-function(message="missing values in object"){
      function(object,...){
        if (NCOL(object)==0)
          object
        else {
          ok <- complete.cases(object)
          if (all(ok)) 
            object
          else stop(message)
        }
      }
    }

     na.id<-na.failsafe("missing values in `id'")
     if(inherits(ids,"formula")) {
     mf<-substitute(model.frame(ids,data=data, na.action=na.id))
     ids<-eval.parent(mf)
         if (ncol(ids)==0) ## formula was ~1
           ids<-data.frame(id=1:nrow(ids))
       } else{
         if (is.null(ids))
           stop("Must provide ids= argument")
         else
           ids<-na.id(data.frame(ids))
       }

    na.prob<-na.failsafe("missing values in `prob'")
    if(inherits(probs,"formula")){
      mf<-substitute(model.frame(probs,data=data,na.action=na.prob))
      probs<-eval.parent(mf)
    }

    na.weight<-na.failsafe("missing values in `weights'")
    if(inherits(weights,"formula")){
      mf<-substitute(model.frame(weights,data=data,na.action=na.weight))
      weights<-eval.parent(mf)
     } else if (!is.null(weights))
         weights<-na.weight(data.frame(weights))
    if(!is.null(weights)){
      if (!is.null(probs))
         stop("Can't specify both sampling weights and probabilities")
       else
         probs<-as.data.frame(1/as.matrix(weights))
     }

      

    na.strata<-na.failsafe("missing values in `strata'")
    if (!is.null(strata)){
      if(inherits(strata,"formula")){
        mf<-substitute(model.frame(strata,data=data, na.action=na.strata))
        strata<-eval.parent(mf)
      }
      if (!is.list(strata))
        strata<-data.frame(strata=strata)
      has.strata<-TRUE
    } else {
      has.strata <-FALSE
      strata<-na.strata(as.data.frame(matrix(1, nrow=NROW(ids), ncol=NCOL(ids))))
    }

    
    if (inherits(variables,"formula")){
        mf<-substitute(model.frame(variables,data=data,na.action=na.pass))
        variables <- eval.parent(mf)
    } else if (is.null(variables)){
        variables<-data
    } else
        variables<-do.call("data.frame",variables)


    na.fpc<-na.failsafe("missing values in `fpc'")
    if (inherits(fpc,"formula")){
      mf<-substitute(model.frame(fpc,data=data,na.action=na.fpc))
      fpc<-eval.parent(mf)
    }
      
      
      ## force subclusters nested in clusters
      if (NCOL(ids)>1){
        N<-ncol(ids)
        for(i in 2:N){
          ids[,i]<-do.call("interaction", ids[,1:i,drop=FALSE])
        }
      }
      ## force clusters nested in strata
      if (nest && has.strata && NCOL(ids)){
        N<-NCOL(ids)
        NS<-NCOL(strata)
        for(i in 1:N)
          ids[,i]<-do.call("interaction",
                           c(strata[,1:min(i,NS),drop=FALSE], ids[,i,drop=FALSE]))
      }
      
    ## check if clusters nested in strata 
     if (check.strata && nest)
      warning("No point in check.strata=TRUE if nest=TRUE")
    if(check.strata && !is.null(strata) && NCOL(ids)){
       sc<-rowSums(table(ids[,1],strata[,1])>0)
       if(any(sc>1)) stop("Clusters not nested in strata at top level; you may want nest=TRUE.")
    }

      ## force substrata nested in clusters
      N<-ncol(ids)
      NS<-ncol(strata)
      if (N>1){
        for(i in 2:N)
          strata[,i]<-interaction(strata[,min(i,NS)], ids[,i-1])
      }
        
    ## Finite population correction: specified per observation
    ## Also incorporates design sample sizes formerly in nPSU

      if (!is.null(fpc) && !is.numeric(fpc) && !is.data.frame(fpc))
        stop("fpc must be a matrix or dataframe or NULL")

      fpc<-as.fpc(fpc,strata, ids)

      ## if FPC specified, but no weights, use it for weights
    if (is.null(probs) && is.null(weights)){
      if (is.null(fpc$popsize)){
        if (missing(probs) && missing(weights))
          warning("No weights or probabilities supplied, assuming equal probability")
        probs<-rep(1,nrow(ids))
      } else {
        probs<-1/weights(fpc, final=FALSE)
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
    rval$call<-sys.call(-1)
    class(rval)<-c("survey.design2","survey.design")
    rval
  }

onestrat<-function(x,cluster,nPSU,fpc, lonely.psu,stratum=NULL,stage=1,cal=cal){
  # Flag lonely strata in subdomain 
  is.lonely <- FALSE

# NOTE: Lonely PSUs treatment under the "adjust" option revised as suggested by Practical
#       Significance blog on 02/09/2022
  stratum_center <- attr(x, "center")
  if (is.null(stratum_center)) stratum_center <- 0  # The old way: centering in 0 (always
                                                    # good for linearized variables)

  if (is.null(fpc))
      f<-rep(1,NROW(x))
  else{
      f<-ifelse(fpc==Inf, 1, (fpc-nPSU)/fpc)
  }

  if (nPSU>1)
      scale<-f*nPSU/(nPSU-1)
  else
      scale<-f

  ## self-representing stratum
  if (all(f < 0.0000001)){
      rval <- matrix(0,NCOL(x),NCOL(x))
      attr(rval, "is.lonely") <- FALSE
      return(rval)
    }
  scale<-scale[!duplicated(cluster)]

  x<-rowsum(x,cluster)
  nsubset<-nrow(x)

  if (nsubset<nPSU) {
    ##can't be PPS, so scale must be a constant
    x<-rbind(x,matrix(0,ncol=ncol(x),nrow=nPSU-nrow(x)))
    scale<-rep(scale[1],NROW(x))
  }
  if (lonely.psu!="adjust" || nsubset>1 ||
      (nPSU>1 & !getOption("RG.adjust.domain.lonely"))) {
      stratum_center <- colMeans(x)
	}
  x<-sweep(x, 2, stratum_center, "-")  # Centering in Y^hat / all_PSUs when the lonely.psu
                                       # option is "adjust"

  if (nsubset==1 && nPSU>1){
      if (isTRUE(getOption("RG.warn.domain.lonely"))){
          warning("Stratum (",stratum,") has only one sampling unit at stage ",stage)
        }
      if (lonely.psu=="average" && getOption("RG.adjust.domain.lonely")) {
          scale<-NA
          is.lonely <- TRUE
        }
    }
  
  if (nPSU>1){
      rval <- crossprod(x*sqrt(scale))
      # Here is.lonely is set to current, rather than FALSE: handles lonelies in
      # subdomains when lonely.psu=="average" && getOption("RG.adjust.domain.lonely")
      attr(rval, "is.lonely") <- is.lonely
      return(rval)
  } else {
      rval<-switch(lonely.psu, 
                   certainty=crossprod(x*sqrt(scale)),
                   remove=crossprod(x*sqrt(scale)),
                   adjust=crossprod(x*sqrt(scale)),
                   # If lonely.psu='average' variance contributions from lonely
                   # strata are mapped to NA: this is subsequently used by onestage
                   # to identify lonely strata and to attach to them the average
                   # variance computed on non-lonely strata...
                   # THIS CREATES WRONG RESULTS WHEN INTEREST VARIABLES HAVE MISSING DATA
                   # BECAUSE THEY TOO ORIGINATE NA STRATUM VARIANCES AND THIS EVENTUALLY AFFECTS
                   # ALSO VARIABLES WITHOUT NAs ESTIMATED TOGHETER WITH THE FORMER...                     
                   average=NA*crossprod(x),
                   fail= stop("Stratum (",stratum,") has only one sampling unit at stage ",stage),
                   stop("Can't handle lonely.psu=",lonely.psu)
            )
      attr(rval, "is.lonely") <- TRUE
      return(rval)
  }
}


onestage<-function(x, strata, clusters, nPSU, fpc, lonely.psu=getOption("RG.lonely.psu"),stage=0, cal){

# NOTE: Lonely PSUs treatment under the "adjust" option revised as suggested by Practical
#       Significance blog on 02/09/2022
   if (lonely.psu=="adjust") {
     # Get the overall number of clusters selected at that stage (e.g., PSUs at stage 1)
     all_PSUs <- sum(tapply(nPSU, as.numeric(strata), FUN = head, 1))
     # Get the weighted total per cluster to use it to compute deviations for lonely PSUs
     center <- colSums(x) / all_PSUs
    } else {
     center <- 0
    }

  stratvars<-tapply(1:NROW(x), list(factor(strata)), function(index){
    onestrat(`attr<-`(x[index,,drop=FALSE], "center", center), clusters[index],
             nPSU[index][1], fpc[index], ##old was fpc[index][1],
             lonely.psu=lonely.psu,stratum=strata[index][1], stage=stage,cal=cal)
    }, simplify = FALSE)
  is.lonely <- sapply(stratvars, function(el) attr(el,"is.lonely"))
  p<-NCOL(x)

# BUG FIX: START
  expand.factor <- 1
  if (lonely.psu=="average"){
      nstrat<-length(unique(strata))
      nokstrat<-sum(!is.lonely)
      if (nokstrat > 0) {
          stratvars <- stratvars[!is.lonely]
          expand.factor <- nstrat/nokstrat
        }
      else {
          stratvars <- lapply(stratvars, function(el) {el[,]<-0; el})
        }
    }
# BUG FIX: END

  apply(array(unlist(stratvars),c(p,p,length(stratvars))),1:2,sum)*expand.factor
}


svyrecvar<-function(x, clusters,  stratas, fpcs, postStrata=NULL, design,
                    lonely.psu=getOption("RG.lonely.psu"),
                    one.stage=getOption("RG.ultimate.cluster")){
####################################################################
# MODIFIED version to accomodate cal.analytic object with their    #
# iterated calibration technique.                                  #
# NOTE: the new function has the new argument 'design'.            #
#       This argument is used twice: 1) in order to asses if the   #
#       design has variance PSUs, 2) in order to asses if it is a  #
#       cal.analytic object.                                       #
####################################################################

  if ( !getOption("RG.ultimate.cluster") && has.var.PSU(design) ){
      one.stage <- TRUE
      #  warning("Ultimate Cluster Approximation assumed for variance estimation on this design!\n
      #           (despite this option isn't currently active, see ?ReGenesees.options)")
    }

  x<-as.matrix(x)
  cal<-NULL
  
  ## Remove post-stratum means, which may cut across clusters
  ## Also center the data using any "g-calibration" models
  if(!is.null(postStrata)){
    for (psvar in postStrata){
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage==0){
          ## G-calibration at population level
          x<-qr.resid(psvar$qr,x/psvar$w)*psvar$w
        } else {
          ## G-calibration within clusters
          cal<-c(cal, list(psvar))
        }
      } else if (inherits(psvar, "analytic_calibration")){
        ####################################################
        # This is the code added to cope with cal.analytic #
        ####################################################
        nobs <- dim(design$variables)[1]
        ## If svyrecvar has been called on a subset get
        ## observations index:
        domain.index <- attr(design, "domain.index")
        ## else build it for the whole sample:
        if (is.null(domain.index)) domain.index <- 1:nobs
        if (psvar$stage!=0)
            stop("Analytic calibration must be at population level")
        ## residuals must have the same structure as x (thus a matrix)
        residuals <- x
        ## ...but filled by zeros
        residuals[, ] <- 0
        ## Now get the list of qr decompositions over partitions (if any)
        qr.list <- postStrata[[1]]$qr.list
        for (part in qr.list) {
             part.index <- part$group
             part.gwhalf <- part$gwhalf
             # Only partitions with non-empty overlap with the estimation domain matter
             if (any(domain.index %in% part.index)) {
                 # Check if any x value is NaN or infinite: this signals that the
                 # estimator gradient (recall x can actually be a linearized variable)
                 # is singular at the Taylor series expansion point.
                 # If this is the case, warn and avoid calling qr.resid (which
                 # would give an error, being Inf NA and NaN unmanageable by Fortran):
                 # return NaN "residuals" instead:
                 z.woodruff <- x[part.index, ,drop=FALSE]
                 if (any( is.infinite(z.woodruff) | is.nan(z.woodruff) )){
                     warning("Estimator gradient is singular at Taylor series expansion point!")
                     residuals[part.index, ] <- NaN
                     }
                 else {
                     residuals[part.index, ] <- qr.resid(part$qr,
                                                         z.woodruff/part.gwhalf) * part.gwhalf
                    }
            }
        }
        x <- residuals
        ####################################################        
        # End                                              #
        ####################################################        
      } else if (inherits(psvar, "raking")){
        ## raking by iterative proportional fitting
        for(iterations in 1:10){
          for(margin in psvar){
            psw<-attr(margin, "weights")
            x<- x - psw*apply(x/psw, 2, ave, margin)
          }
        }
      } else {
        ## ordinary post-stratification
        psw<-attr(psvar, "weights")
        oldw<-attr(psvar, "oldweights")
        if (is.null(oldw)) oldw<-rep(1,length(psw))
        psvar<-as.factor(psvar)
        psmeans<-rowsum(x*oldw/psw,psvar,reorder=TRUE)/as.vector(by(oldw,psvar,sum))
        x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
      }
    }
  }
  
  multistage(x, clusters,stratas,fpcs$sampsize, fpcs$popsize,
             lonely.psu=getOption("RG.lonely.psu"),
             one.stage=one.stage,stage=1,cal=cal)
}

multistage<-function(x, clusters,  stratas, nPSUs, fpcs,
                    lonely.psu=getOption("RG.lonely.psu"),
                     one.stage=FALSE,stage,cal){
  
  n<-NROW(x)
 
  
  v <- onestage(x,stratas[,1], clusters[,1], nPSUs[,1],
                fpcs[,1], lonely.psu=lonely.psu,stage=stage,cal=cal)
  
  if (one.stage!=TRUE && !is.null(fpcs) && NCOL(clusters)>1) {
  # FIX: old had as.numeric(cluster[,1]) but psus ids MAY BE character and such
  #      that as.numeric generates NAs 27/04/2016
  # v.sub<-by(1:n, list(as.numeric(clusters[,1])), function(index){
    v.sub<-by(1:n, list(clusters[,1]), function(index){
      ## residuals for G-calibration using population information
      ## only on clusters at this stage.
      for(cali in cal){
        if (cali$stage != stage)
          next
        j<-match(clusters[index,1],cali$index)
        if (length(unique(j))!=1)
          stop("Internal problem in g-calibration data: stage",stage,
               ", cluster", j)
        j<-j[[1]]
        x[index,]<-qr.resid(cali$qr[[j]], x[index,,drop=FALSE]/cali$w[[j]])*cali$w[[j]]
      }
      multistage(x[index,,drop=FALSE], clusters[index,-1,drop=FALSE],
                 stratas[index,-1,drop=FALSE], nPSUs[index,-1,drop=FALSE],
                 fpcs[index,-1,drop=FALSE],
                 lonely.psu=lonely.psu,one.stage=one.stage-1,
                 stage=stage+1,cal=cal)*nPSUs[index[1],1]/fpcs[index[1],1]
    })
    
    for(i in 1:length(v.sub))
      v<-v+v.sub[[i]]
  }
  dimnames(v)<-list(colnames(x),colnames(x))
  v
}


## fpc not given are zero: full sampling.
as.fpc<-function(df,strata,ids){

  count<-function(x) length(unique(x))
  
  sampsize<-matrix(ncol=ncol(ids),nrow=nrow(ids))
  #################################################
  # Next statement added to save computation time #
  #################################################
  ids <- as.matrix(ids)

  for(i in 1:ncol(ids))
    split(sampsize[,i],strata[,i])<-lapply(split(ids[,i],strata[,i]),count)
  
  if (is.null(df)){
    ## No fpc
    rval<-list(popsize=NULL, sampsize=sampsize)
    class(rval)="survey_fpc"
    return(rval)
  }
  
  fpc<-as.matrix(df)
  #Original xor was not correct when sample -> population
  if ( ( ispopsize <- any(df>1) ) && !all(df>=1) ) {
    big<-which(fpc > 1,arr.ind=TRUE) # -> FIX old was big<-which(fpc>=1,arr.ind=TRUE) 27/04/2016
    small<-which(fpc<1,arr.ind=TRUE)
    warning("record ",  big[1,1]," stage ",  big[1,2]," : fpc= ", fpc[  big[1,,drop=FALSE]])
    warning("record ",small[1,1]," stage ",small[1,2]," : fpc= ", fpc[small[1,,drop=FALSE]])      
    stop("Must have all fpc>=1 or all fpc<=1")
  }

# TO CHECK!!!  (lacking in survey, though documented:
#              "If \code{fpc} is specified but for fewer stages than \code{ids}
#              sampling is assumed to be complete for subsequent stages.")  
  if ( (nfpc <- ncol(fpc)) < (stages <- ncol(sampsize)) ) {
    fpc.new <- sampsize
    fpc.new[, 1:nfpc] <- fpc[, 1:nfpc]
    if (ispopsize) {
       fpc.new[, (nfpc+1):stages] <- sampsize[, (nfpc+1):stages]
    } else {
       fpc.new[, (nfpc+1):stages] <- 1
    }
    colnames(fpc.new) <- c(colnames(fpc), paste("sFPC.", (nfpc+1):stages, sep=""))
    fpc <- fpc.new
  }
# END
  if (ispopsize){
    popsize<-fpc
  } else {
    popsize<-sampsize/(fpc)
  }
  if (any(popsize<sampsize)){
    toobig<-which(popsize<sampsize,arr.ind=TRUE)
    warning("record ",toobig[1,1]," stage ",toobig[1,2]," : popsize= ",popsize[toobig[1,,drop=FALSE]],
        " sampsize= ", sampsize[toobig[1,,drop=FALSE]])
    stop("FPC implies >100% sampling in some strata")
  }
  if (!ispopsize && any(is.finite(popsize) & (popsize>1e10))){
    big<-which(popsize>1e10 & is.finite(popsize),arr.ind=TRUE)
    warning("FPC implies population larger than ten billion (record ",big[1,1]," stage ",big[1,2],")")
  }

  ## check that fpc is constant within strata.
  for(i in 1:ncol(popsize)){
    diff<-by(popsize[,i], list(strata[,i]), count)
    if (any(as.vector(diff)>1)){
        j<-which(as.vector(diff)>1)[1]
        stop("'fpc' values vary within strata: stratum ",names(diff)[j], " at stage ",i)
    }
  }
  
  rval<-list(popsize=popsize, sampsize=sampsize)
  class(rval)<-"survey_fpc"
  rval
}

"weights.survey_fpc"<-function(object,final=TRUE,...){
  if (is.null(object$popsize) || any(object$popsize>1e12))
    stop("Weights not supplied and can't be computed from fpc.")
  if (final) {
    pop<-apply(object$popsize,1,prod)
    samp<-apply(object$sampsize,1,prod)
    pop/samp
  } else {
    object$popsize/object$sampsize
  }
}

prcall<-function (xcall)
#########################################
# Prints the call stack which generated #
# a survey.design2 object.              #
# Note: Handles multiple cal.analytic   #
#       objects, both arising from the  #
#       command line and from the GUI.  #
######################################### 
{
    if (!is.list(xcall)) {
        if (is.call(xcall)){
            print(xcall)
            }
      else  {
            # cat(xcall, "\n")
            print(as.call(parse(text=xcall))[[1]])
            }
    }
    else {
        n <- length(xcall)
        for (i in 1L:n) {
            label <- paste(n - i + 1L, ": ", sep = "")
            if (is.call(xcall[[i]])){
                cat(label)
                print(xcall[[i]])
                }
            else{
                # cat(paste(label, xcall[[i]], sep = ""), sep = "\n")
                cat(label)
                print(as.call(parse(text=xcall[[i]]))[[1]])
                }
        }
    }
    invisible()
}
    

print.survey.design2<-function(x,varnames=FALSE,design.summaries=FALSE,...){
  n<-NROW(x$cluster)
  if (x$has.strata) cat("Stratified ")
  un<-length(unique(x$cluster[,1]))
  if(n==un){
    cat("Independent Sampling design")
    is.independent<-TRUE
    if (is.null(x$fpc$popsize))
      cat(" (with replacement)\n")
    else cat("\n")
  } else {
    cat(NCOL(x$cluster),"- level Cluster Sampling design")
    if (is.null(x$fpc$popsize))
      cat(" (with replacement)\n")
    else cat("\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
    is.independent<-FALSE
  }

  prcall(x$call)
  
  if (design.summaries){
    cat("Probabilities:\n")
    print(summary(x$prob))
    if(x$has.strata){
      if (NCOL(x$cluster)>1)
        cat("First-level ")
      cat("Stratum Sizes: \n")
      oo<-order(unique(x$strata[,1]))
      a<-rbind(obs=table(x$strata[,1]),
           design.PSU=x$fpc$sampsize[!duplicated(x$strata[,1]),1][oo],
               actual.PSU=table(x$strata[!duplicated(x$cluster[,1]),1]))
      print(a)
    }
    if (!is.null(x$fpc$popsize)){
      if (x$has.strata) {
        cat("Population stratum sizes (PSUs): \n")
        s<-!duplicated(x$strata[,1])
        a<-x$fpc$popsize[s,1]
        names(a)<-x$strata[s,1]
        print(a)
      } else {
        cat("Population size (PSUs):",x$fpc$popsize[1,1],"\n")
      }
    }
  }
  if (varnames){
    cat("Data variables:\n")
    print(colnames(x))
  }
  invisible(x)
}

print.analytic <- function(x, varnames = FALSE, design.summaries = FALSE, ...){

   if (inherits(x, "cal.analytic")) {
       cal <- "Calibrated, "
    }
   else {
        cal <- NULL    
    }
  cat(cal)
  n <- NROW(x$cluster)
  if (x$has.strata) {
      cat("Stratified ")
      nStrata <- length(unique(x$strata[, 1]))
    }
  un <- length(unique(x$cluster[,1]))
  if(n == un){
    cat("Independent Unit Sampling Design")
#   is.independent<-TRUE

#   If first stage fpc is always zero or Inf, print "(with replacement)",
#   as happened when PSU's fpc was not specified (17/11/2015)
    if ( is.null(x$fpc$popsize) || all(is.infinite(x$fpc$popsize[, 1])) )
      cat(" (with replacement)\n")
    else cat("\n")
    if (x$has.strata) {
         if (is.null(attr(x, "collapse.strata"))) {
             cat(paste("- [", nStrata , "] strata\n", sep = ""))
            }
         else {
             cat(paste("- [", nStrata , "] strata (collapsed)\n", sep = ""))
            }
        }
    cat(paste("- [",un,"] units\n", sep = ""))
  } else {
    cat(NCOL(x$cluster),"- Stage Cluster Sampling Design")
    if ( is.null(x$fpc$popsize) || all(is.infinite(x$fpc$popsize[, 1])) )
      cat(" (with replacement)\n")
    else cat("\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    if (x$has.strata) {
         if (is.null(attr(x, "collapse.strata"))) {
             cat(paste("- [", nStrata , "] strata\n", sep = ""))
            }
         else {
             cat(paste("- [", nStrata , "] strata (collapsed)\n", sep = ""))
            }
        }
    cat(paste("- [",paste(unlist(nn),collapse=", "),"] clusters\n", sep = ""))
#   is.independent<-FALSE
  }
  cat("\nCall:\n")
  prcall(x$call)

  if (design.summaries){
    if (is.null(cal)){
         cat("\nProbabilities:\n")
         print(summary(x$prob))
        }
    else {
         # As probabilities are ill defined for calibrated objects...
         cat("\nWeights:\n")
         print(summary(weights(x)))
        }
    if(x$has.strata){
      # if (NCOL(x$cluster)>1)
      cat("\nSample stratum sizes: \n")
      oo<-order(unique(x$strata[,1]))
      if (n == un){
         a<-rbind(obs=table(x$strata[,1]))
        }
      else {
         a<-rbind(Obs=table(x$strata[,1]),
                  PSUs=x$fpc$sampsize[!duplicated(x$strata[,1]),1][oo])
                  # Old code: design.PSUs could DIFFER from actual.PSUs only due
                  # to SUBSETTING on x. I find it not useful.
                  # design.PSUs=x$fpc$sampsize[!duplicated(x$strata[,1]),1][oo],
                  # actual.PSUs=table(x$strata[!duplicated(x$cluster[,1]),1]))
        }
      print(a)
    }
      # NOTE: If design is NOT stratified, the OVERALL sample size in terms
      #       of PSUs, SSUs is already provided by the call to print() above
    if (!is.null(x$fpc$popsize)){
      if (x$has.strata) {
        # cat("\nPopulation stratum sizes (PSUs): \n")
        cat("\nPopulation stratum sizes: \n")
        s<-!duplicated(x$strata[,1])
        a<-x$fpc$popsize[s,1]
        names(a)<-x$strata[s,1]
        a <- if (n == un) rbind(Obs = a) else rbind(PSUs = a) # Better vertical alignement with sample sizes
        print(a)
      } else {
        cat("\nPopulation size (PSUs):",x$fpc$popsize[1,1],"\n")
      }
    }
  }
  if (varnames){
    cat("\nData variables:\n")
    print(colnames(x))
  }
  invisible(x)
}
    
summary.survey.design2<-function(object,...){
  class(object)<-c("summary.survey.design2",class(object))
  object
}

print.summary.survey.design2<-function(x,...){
  y<-x
  class(y)<-c("survey.design2",class(x))
  print(y,varnames=TRUE,design.summaries=TRUE,...)
}    


summary.analytic <- function(object, ...){
  class(object) <- c("summary.analytic", class(object))
  object
}
 
print.summary.analytic <- function(x, ...){
  y <- x
  class(y) <- c("analytic", class(x))
  print(y, varnames = TRUE, design.summaries = TRUE, ...)
}    


.svycheck<-function(object){
  if (inherits(object,"survey.design") &&
      !is.null(object$nPSU))
    warning("This is an old-style design object. Please use as.svydesign2 to update it.")
}

as.svydesign2<-function(object){
  if (inherits(object,"survey.design2"))
    return(object)
  if (!inherits(object,"survey.design"))
    stop("This function is for updating old-style survey.design objects")
  

  count<-function(x) length(unique(x))
  
  strata<-data.frame(one=object$strata)
  if ((nc<-ncol(object$cluster))>1){
    for(i in 2:nc){
      strata<-cbind(strata,object$cluster[,i-1])
    }
  }
  
  sampsize<-matrix(ncol=nc,nrow=nrow(object$cluster))
  
  sampsize[,1]<-object$nPSU[match(object$strata, names(object$nPSU))]
  if (nc>1){
    for(i in 2:nc){
      split(sampsize[,i],strata[,i])<-lapply(split(object$cluster[,i],strata[,i]),count)
    }
  }
  
  if (!is.null(object$fpc)){
    popsize<-sampsize
    popsize[,1]<-object$fpc$N[match(object$strata,object$fpc$strata)]
  } else popsize<-NULL
  if (nc>1 && !is.null(object$fpc)){
    warning("Assuming complete sampling at stages 2 -",nc)
  }

  fpc<-list(popsize=popsize,sampsize=sampsize)
  class(fpc)<-"survey_fpc"
  
           
  object$fpc<-fpc
  object$strata<-strata
  object$nPSU<-NULL
  class(object)<-c("survey.design2","survey.design")
  object
  
}

    
"[.survey.design2"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
      if (is.calibrated(x) || !drop){
          ## Set weights to zero: no memory saving possible
          ## There should be an easier way to complement a subscript..
          if (is.logical(i))
              x$prob[!i]<-Inf
          else if (is.numeric(i) && length(i))
              x$prob[-i]<-Inf
          else {
              tmp<-x$prob[i,]
              x$prob<-rep(Inf, length(x$prob))
              x$prob[i,]<-tmp
          }
          index<-is.finite(x$prob)
          psu<-!duplicated(x$cluster[index,1])
          tt<-table(x$strata[index,1][psu])
          if (isTRUE(getOption("RG.warn.domain.lonely"))){
              if(any(tt==1)){
                 warning(sum(tt==1)," strata have only one PSU in this subset.")
                }
            }
      } else {
          ## subset everything.
          if (!is.null(x$variables)) ## phase 2 of twophase design
              x$variables<-"[.data.frame"(x$variables,i,..1,drop=FALSE)
          x$cluster<-x$cluster[i,,drop=FALSE]
          x$prob<-x$prob[i]
          x$allprob<-x$allprob[i,,drop=FALSE]
          x$strata<-x$strata[i,,drop=FALSE]
          x$fpc$sampsize<-x$fpc$sampsize[i,,drop=FALSE]
          x$fpc$popsize<-x$fpc$popsize[i,,drop=FALSE]
      }
      
  } else {
      if(!is.null(x$variables))
          x$variables<-x$variables[,..1,drop=FALSE]
  }
  
  x
}


"[.cal.analytic"<-function (x,i, ..., drop=TRUE){
#############################################
# Indexing method for cal.analytic objects. #
# This method is needed in svyby for domain #
# estimation.                               #
#############################################
  if (!missing(i)){ 
          ## Set weights to zero: no memory saving possible
          if (is.logical(i))
              x$prob[!i]<-Inf
          else if (is.numeric(i) && length(i))
              x$prob[-i]<-Inf
          else {
              tmp<-x$prob[i,]
              x$prob<-rep(Inf, length(x$prob))
              x$prob[i,]<-tmp
          }
          index<-is.finite(x$prob)
          ## Attach to the design object a new attribute
          attr(x,"domain.index") <- which(index)
          psu<-!duplicated(x$cluster[index,1])
          tt<-table(x$strata[index,1][psu])
          if (isTRUE(getOption("RG.warn.domain.lonely"))){
              if( any(tt==1) && (sum(index) < nrow(x)) ){
                 warning(sum(tt==1)," strata have only one PSU in this subset.")
                }
            }

  } else {
      if (!is.null(x$variables)) {
          x$variables<-x$variables[,..1,drop=FALSE]
          attr(x,"domain.index") <- 1:(dim(x$variables)[1])
        }
  }
  
  x
}




svytotal.survey.design2<-function(x,design, na.rm=FALSE, deff=FALSE,...){

  
    if (inherits(x,"formula")){
        ## do the right thing with factors
        mf<-model.frame(x,design$variables,na.action=na.pass)
        #PROVA! QUI CON MODIFICA POTREI TRATTARE FACTOR CON 1 SOLO LIVELLO!
        xx<-lapply(attr(terms(x),"variables")[-1],
                   function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
        cols<-sapply(xx,NCOL)
        x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
        scols<-c(0,cumsum(cols))
        for(i in 1:length(xx)){
            x[,scols[i]+1:cols[i]]<-xx[[i]]
        }
        colnames(x)<-do.call("c",lapply(xx,colnames))
    } else{
        if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        else {
            if(is.data.frame(x) && any(sapply(x,is.factor))){
                xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
                cols<-sapply(xx,NCOL)
                scols<-c(0,cumsum(cols))
                cn<-character(sum(cols))
                for(i in 1:length(xx))
                    cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
                x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
                for(i in 1:length(xx)){
                    x[,scols[i]+1:cols[i]]<-xx[[i]]
                }
                colnames(x)<-cn
            }
        }
    }
    x<-as.matrix(x)

    ## Missing values treatment
    nas<-rowSums(is.na(x))
    if (na.rm && sum(nas)>0){
        # If domain has some non-NA values, use them for estimation:
        if (length(design[nas==0,]$prob) > 0) {
            # Will use a model-based estimator, assuming MCAR: Y' = (Y.na/N.na)*N.all
            N.all<-sum(1/design$prob)
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

    N<-sum(1/design$prob)
    if (na.rm && sum(nas)>0) {
        # Will use a model-based estimator, assuming MCAR: Y' = (Y.na/N.na)*N.all
        x <- x*(N.all/N)
    }
    total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
    class(total)<-"svystat"

    attr(total, "var")<-v<-svyrecvar(x/design$prob,design$cluster,
                                     design$strata, design$fpc,
                                     postStrata=design$postStrata, design = design)
    attr(total,"statistic")<-"total"

    if (is.character(deff) || deff){
    ###################################
    # Here z.svyvar instead of svyvar #
    ###################################
      nobs<-NROW(design$cluster)
      if (deff=="replace")
        vsrs<-z.svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2/nobs
      else
        vsrs<-z.svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
      attr(total, "deff")<-v/vsrs
    }
    

  return(total)
}


svytotal.cal.analytic<-function(x,design, na.rm=FALSE, deff=FALSE,...){
#############################################
# svytotal method for cal.analytic objects. #
# NOTE: svyrecvar is called with the new    #
#       argument design.                    #
#############################################

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
    } else{
        if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        else {
            if(is.data.frame(x) && any(sapply(x,is.factor))){
                xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
                cols<-sapply(xx,NCOL)
                scols<-c(0,cumsum(cols))
                cn<-character(sum(cols))
                for(i in 1:length(xx))
                    cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
                x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
                for(i in 1:length(xx)){
                    x[,scols[i]+1:cols[i]]<-xx[[i]]
                }
                colnames(x)<-cn
            }
        }
    }
    x<-as.matrix(x)

    
    ## Missing values treatment
    nas<-rowSums(is.na(x))
    if (na.rm && sum(nas)>0){
        # If domain has some non-NA values, use them for estimation:
        if (length(design[nas==0,]$prob) > 0) {
            # Will use a model-based estimator, assuming MCAR: Y' = (Y.na/N.na)*N.all
            N.all<-sum(1/design$prob)
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

    N<-sum(1/design$prob)
    if (na.rm && sum(nas)>0) {
        # Will use a model-based estimator, assuming MCAR: Y' = (Y.na/N.na)*N.all
        x <- x*(N.all/N)
    }
    total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
    class(total)<-"svystat"

    attr(total, "var")<-v<-svyrecvar(x/design$prob,design$cluster,
                                     design$strata, design$fpc,
                                     postStrata=design$postStrata, design = design)
    attr(total,"statistic")<-"total"

    if (is.character(deff) || deff){
    ## If svytotal has been called on a subset get nobs from the domain index
    ## else compute it for the whole sample:
      if (is.null(di <- attr(design, "domain.index"))) nobs<-NROW(design$cluster) else nobs <- length(di)
      ###################################
      # Here z.svyvar instead of svyvar #
      ###################################
      if (deff=="replace") {
          vsrs<-z.svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2/nobs
          }
      else {
          if (N < nobs) {
              vsrs<-NA*v
              warning("Sample size greater than population size: are weights correctly scaled?")
              }
          else {
              vsrs<-z.svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
               }
            }
      # Check for possible anomalies arising from negative calibrated weights (if any)
      if (inherits(vsrs, "matrix")){
           # this happens for factor interest variables
           anom <- ( !is.na(vsrs) & (vsrs < 0) & (row(vsrs)==col(vsrs)) )
      }
      else {
           # this happens for numeric interest variables
           anom <- (!is.na(vsrs) & (vsrs < 0))
      }
      if (any(anom)) {
           vsrs[anom] <- NA
           warning("Negative Deff estimate due to negative calibration weights: returning NA")
      }
      attr(total, "deff")<-v/vsrs
    }
  return(total)
}

svymean.survey.design2<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
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
  else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else if(is.data.frame(x) && any(sapply(x,is.factor))){
          xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
          cols<-sapply(xx,NCOL)
          scols<-c(0,cumsum(cols))
          cn<-character(sum(cols))
          for(i in 1:length(xx))
              cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
          x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
          for(i in 1:length(xx)){
              x[,scols[i]+1:cols[i]]<-xx[[i]]
          }
          colnames(x)<-cn
      }
    }
  x<-as.matrix(x)
  
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
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  
  v<-svyrecvar(x*pweights/psum,design$cluster,design$strata, design$fpc,
              postStrata=design$postStrata, design = design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
  ###################################
  # Here z.svyvar instead of svyvar #
  ###################################
      nobs<-NROW(design$cluster)
      if(deff=="replace"){
        vsrs<-z.svyvar(x,design,na.rm=na.rm)/(nobs)
      } else {
        if(psum<nobs) {
          vsrs<-NA*v
          warning("Sample size greater than population size: are weights correctly scaled?")
        } else{
          vsrs<-z.svyvar(x,design,na.rm=na.rm)*(psum-nobs)/(psum*nobs)
        }
      }
      attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}


svymean.cal.analytic<-function(x,design, na.rm=FALSE,deff=FALSE,...){
#############################################
# svymean method for cal.analytic objects.  #
# NOTE: svyrecvar is called with the new    #
#       argument design.                    #
#############################################
  
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
  else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else if(is.data.frame(x) && any(sapply(x,is.factor))){
          xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
          cols<-sapply(xx,NCOL)
          scols<-c(0,cumsum(cols))
          cn<-character(sum(cols))
          for(i in 1:length(xx))
              cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
          x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
          for(i in 1:length(xx)){
              x[,scols[i]+1:cols[i]]<-xx[[i]]
          }
          colnames(x)<-cn
      }
    }
  x<-as.matrix(x)
  
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
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  
  v<-svyrecvar(x*pweights/psum,design$cluster,design$strata, design$fpc,
              postStrata=design$postStrata, design = design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
  ## If svymean has been called on a subset get nobs from the domain index
  ## else compute it for the whole sample:
      if (is.null(di <- attr(design, "domain.index"))) nobs<-NROW(design$cluster) else nobs <- length(di)
      ##################################
      # Here z.svyvar intead of svyvar #
      ##################################
      if(deff=="replace"){
        vsrs<-z.svyvar(x,design,na.rm=na.rm)/(nobs)
      } else {
        if(psum < nobs) {
          vsrs<-NA*v
          warning("Sample size greater than population size: are weights correctly scaled?")
        } else {
          vsrs<-z.svyvar(x,design,na.rm=na.rm)*(psum-nobs)/(psum*nobs)
        }
      }
      # Check for possible anomalies arising from negative calibrated weights (if any)
      if (inherits(vsrs, "matrix")){
           # this happens for factor interest variables
           anom <- ( !is.na(vsrs) & (vsrs < 0) & (row(vsrs)==col(vsrs)) )
      }
      else {
           # this happens for numeric interest variables
           anom <- (!is.na(vsrs) & (vsrs < 0))
      }
      if (any(anom)) {
           vsrs[anom] <- NA
           warning("Negative Deff estimate due to negative calibration weights: returning NA")
      }
      attr(average, "deff")<-v/vsrs
  }
  return(average)
}


svyratio.survey.design2<-function(numerator=formula, denominator, design, 
    separate=FALSE,na.rm=FALSE,formula,covmat=FALSE,cross=FALSE,deff=FALSE,...){

    if (separate){
      strats<-sort(unique(design$strata[,1]))
      if (!design$has.strata)
        warning("Separate and combined ratio estimators are the same for unstratified designs")
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$strata[,1] %in% s),
                                   separate=FALSE,...)
                     attr(tmp,"call")<-bquote(Stratum==.(s))
                     tmp}))
      names(rval$ratios)<-strats
   
      class(rval)<-c("svyratio_separate")
      rval$call<-sys.call()
      rval$strata<-strats
      return(rval)
    }
  
    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,design$variables,na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,design$variables,na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    numerator<-as.matrix(numerator)
    denominator<-as.matrix(denominator)
    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)

    ## Missing values treatment
    nas<-!complete.cases(all)
    if (na.rm && any(nas)){
      design<-design[!nas,]
      # If domain has some non-NA values, use them for estimation:
      if (length(design$prob) > 0) {
          if (NROW(design$cluster) == NROW(all)){
            ## subset by zero weights
            all[nas,]<-1
            numerator[nas,]<-0
            denominator[nas,]<-1
          } else {
            ## subset by actually dropping rows
            all<-all[!nas,,drop=FALSE]
            numerator<-numerator[!nas,,drop=FALSE]
            denominator<-denominator[!nas,,drop=FALSE]
          }
        }
      else na.rm <- FALSE
    }
    #######################################
    # Here z.svytotal instead of svytotal #
    #######################################
    allstats <- z.svytotal(all,design)

    if (ask.deff <- (is.character(deff) || deff)){
         nobs<-NROW(design$cluster)
         N<-sum(1/design$prob)
        }

    if (cross){
        # Compute estimates and error for all pairs
        rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))

        # Initialize matrix of Variances
        vars<-matrix(ncol=nd,nrow=nn)

        if (ask.deff){
             # Initialize matrix of Deffs
             deffs<-matrix(ncol=nd,nrow=nn)
            }

        for(i in 1:nn){
             for(j in 1:nd){
                 r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)

                 vars[i,j]<-svyrecvar(r*1/design$prob, design$cluster, design$strata, design$fpc,
                                      postStrata=design$postStrata, design = design)

                 if (ask.deff){
                 ###################################
                 # Here z.svyvar instead of svyvar #
                 ###################################
                         if (deff=="replace"){
                             vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2/nobs
                            }
                         else {
                             vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
                            }
                         deffs[i,j] <- vars[i,j]/vsrs
                    }
                }
            }

        if (covmat){
            ii<-rep(1:nn,nd)
            jj<-rep(1:nd,each=nn)
            allr<-sweep(numerator[,ii]-t(as.vector(rval$ratio)*t(denominator[,jj,drop=FALSE])),
                        2, colSums(denominator[,jj,drop=FALSE]/design$prob),"/")
            vcovmat<-svyrecvar(allr*1/design$prob, design$cluster, design$strata, design$fpc,
                               postStrata=design$postStrata, design = design)
            colnames(vcovmat)<-colnames(denominator)[ii]
            rval$vcov<-vcovmat
        }
        colnames(vars)<-colnames(denominator)
        rownames(vars)<-colnames(numerator)
        if (ask.deff){
             colnames(deffs)<-colnames(denominator)
             rownames(deffs)<-colnames(numerator)
            }
    }
    else{
        # Compute parallel estimates and error
        # Totals first:
        numT <- allstats[1:nn]
        denT <- allstats[nn+1:nd]
        # Let mapply do the work of recycling (if needed)... 
        rval<-list(ratio=mapply("/", numT, denT))
        nRatios <- length(rval$ratio)
        # Now let's work with recycled, equal length totals...
        numT <- rep(numT, length.out = nRatios)
        denT <- rep(denT, length.out = nRatios)
        # ... as well as with recycled, equal length matrices storing the variables...
         # first, pick the names
         num.names <- rep(colnames(numerator), length.out = nRatios)
         den.names <- rep(colnames(denominator), length.out = nRatios)
         Rnames <- paste(num.names,"/",den.names, sep="")
         names(rval$ratio) <- Rnames
         # then, extend matrices
         numerator <- numerator[, rep(1:nn, length.out = nRatios), drop=FALSE]
         denominator <- denominator[, rep(1:nd, length.out = nRatios), drop=FALSE]
        
        # vars now is a vector, initialized to NA...
        vars <- vector(mode = "logical", length = nRatios)
        vars[] <- NA

        if (ask.deff){
             # Initialize vector of Deffs
             deffs <- vector(mode = "logical", length = nRatios)
             deffs[] <- NA
            }
        
        for(i in 1:nRatios){
             r<-(numerator[,i]-rval$ratio[i]*denominator[,i])/sum(denominator[,i]/design$prob)

             vars[i]<-svyrecvar(r*1/design$prob, design$cluster, design$strata, design$fpc,
                                postStrata=design$postStrata, design = design)

             if (ask.deff){
             ###################################
             # Here z.svyvar instead of svyvar #
             ###################################
                 if (deff=="replace"){
                     vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2/nobs
                    }
                 else {
                     vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
                    }
                 deffs[i] <- vars[i]/vsrs
                }
            }
        names(vars) <- Rnames
        if (ask.deff){
             names(deffs) <- Rnames
            }
        if (covmat){
            ii<-rep(1:nRatios)
            jj<-rep(1:nRatios)
            allr<-sweep(numerator[,ii]-t(as.vector(rval$ratio)*t(denominator[,jj,drop=FALSE])),
                        2, colSums(denominator[,jj,drop=FALSE]/design$prob),"/")
            vcovmat<-svyrecvar(allr*1/design$prob, design$cluster, design$strata, design$fpc,
                               postStrata=design$postStrata, design = design)
            colnames(vcovmat) <- rownames(vcovmat) <- Rnames
            rval$vcov<-vcovmat
        }
    }
    
    rval$var<-vars
    attr(rval,"cross") <- cross
    if (ask.deff){
         attr(rval, "deff") <- deffs
        }
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
  }


svyratio.cal.analytic<-function(numerator=formula, denominator, design,
                       separate=FALSE,na.rm=FALSE,formula,covmat=FALSE,cross=FALSE,deff=FALSE,...){
#############################################
# svyratio method for cal.analytic objects. #
# NOTE: svyrecvar is called with the new    #
#       argument design.                    #
#############################################

    if (separate){
      strats<-sort(unique(design$strata[,1]))
      if (!design$has.strata)
        warning("Separate and combined ratio estimators are the same for unstratified designs")
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$strata[,1] %in% s),
                                   separate=FALSE,...)
                     attr(tmp,"call")<-bquote(Stratum==.(s))
                     tmp}))
      names(rval$ratios)<-strats
   
      class(rval)<-c("svyratio_separate")
      rval$call<-sys.call()
      rval$strata<-strats
      return(rval)
    }
  
    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,design$variables,na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,design$variables,na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    numerator<-as.matrix(numerator)
    denominator<-as.matrix(denominator)
    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)

    ## Missing values treatment
    nas<-!complete.cases(all)
    if (na.rm && any(nas)){
      design<-design[!nas,]
      # If domain has some non-NA values, use them for estimation:
      if (length(design$prob) > 0) {
          if (NROW(design$cluster) == NROW(all)){
            ## subset by zero weights
            all[nas,]<-1
            numerator[nas,]<-0
            denominator[nas,]<-1
          } else {
            ## subset by actually dropping rows
            all<-all[!nas,,drop=FALSE]
            numerator<-numerator[!nas,,drop=FALSE]
            denominator<-denominator[!nas,,drop=FALSE]
          }
        }
      else na.rm <- FALSE
    }
    #######################################
    # Here z.svytotal instead of svytotal #
    #######################################
    allstats <- z.svytotal(all,design)

    if (ask.deff <- (is.character(deff) || deff)){
         N<-sum(1/design$prob)
        }

    if (cross){
        # Compute estimates and error for all pairs
        rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))

        # Initialize matrix of Variances
        vars<-matrix(ncol=nd,nrow=nn)

        if (ask.deff){
             # Initialize matrix of Deffs
             deffs<-matrix(ncol=nd,nrow=nn)
            }

        for(i in 1:nn){
             for(j in 1:nd){
                 r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        
                 vars[i,j]<- v <- svyrecvar(r*1/design$prob, design$cluster, design$strata, design$fpc,
                                            postStrata=design$postStrata, design = design)

                 if (ask.deff){
                     ## If svytotal has been called on a subset get nobs from the domain index
                     ## else compute it for the whole sample:
                     if (is.null(di <- attr(design, "domain.index"))) {
                         nobs<-NROW(design$cluster)
                        }
                     else {
                         nobs <- length(di)
                        }
                     ###################################
                      # Here z.svyvar instead of svyvar #
                     ###################################
                     if (deff=="replace") {
                         vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2/nobs
                        }
                     else {
                         if (N < nobs) {
                             vsrs<-NA*v
                             warning("Sample size greater than population size: are weights correctly scaled?")
                            }
                         else {
                             vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
                            }
                        }
                     # Check for possible anomalies arising from negative calibrated weights (if any)
                     if (inherits(vsrs, "matrix")){
                          # this happens for factor interest variables
                          anom <- ( !is.na(vsrs) & (vsrs < 0) & (row(vsrs)==col(vsrs)) )
                     }
                     else {
                          # this happens for numeric interest variables
                          anom <- (!is.na(vsrs) & (vsrs < 0))
                     }
                     if (any(anom)) {
                          vsrs[anom] <- NA
                          warning("Negative Deff estimate due to negative calibration weights: returning NA")
                     }
                     deffs[i,j] <- vars[i,j]/vsrs
                    }
                }
            }
        if (covmat){
            ii<-rep(1:nn,nd)
            jj<-rep(1:nd,each=nn)
            allr<-sweep(numerator[,ii]-t(as.vector(rval$ratio)*t(denominator[,jj,drop=FALSE])),
                        2, colSums(denominator[,jj,drop=FALSE]/design$prob),"/")
            vcovmat<-svyrecvar(allr*1/design$prob, design$cluster, design$strata, design$fpc,
                               postStrata=design$postStrata, design = design)
            colnames(vcovmat)<-colnames(denominator)[ii]
            rval$vcov<-vcovmat
        }
        colnames(vars)<-colnames(denominator)
        rownames(vars)<-colnames(numerator)
        if (ask.deff){
             colnames(deffs)<-colnames(denominator)
             rownames(deffs)<-colnames(numerator)
            }
    }
    else{
        # Compute parallel estimates and error
        # Totals first:
        numT <- allstats[1:nn]
        denT <- allstats[nn+1:nd]
        # Let mapply do the work of recycling (if needed)... 
        rval<-list(ratio=mapply("/", numT, denT))
        nRatios <- length(rval$ratio)
        # Now let's work with recycled, equal length totals...
        numT <- rep(numT, length.out = nRatios)
        denT <- rep(denT, length.out = nRatios)
        # ... as well as with recycled, equal length matrices storing the variables...
         # first, pick the names
         num.names <- rep(colnames(numerator), length.out = nRatios)
         den.names <- rep(colnames(denominator), length.out = nRatios)
         Rnames <- paste(num.names,"/",den.names, sep="")
         names(rval$ratio) <- Rnames
         # then, extend matrices
         numerator <- numerator[, rep(1:nn, length.out = nRatios), drop=FALSE]
         denominator <- denominator[, rep(1:nd, length.out = nRatios), drop=FALSE]
        
        # vars now is a vector, initialized to NA...
        vars <- vector(mode = "logical", length = nRatios)
        vars[] <- NA

        if (ask.deff){
             # Initialize vector of Deffs
             deffs <- vector(mode = "logical", length = nRatios)
             deffs[] <- NA
            }

        for(i in 1:nRatios){
             r<-(numerator[,i]-rval$ratio[i]*denominator[,i])/sum(denominator[,i]/design$prob)

             vars[i]<- v <- svyrecvar(r*1/design$prob, design$cluster, design$strata, design$fpc,
                                      postStrata=design$postStrata, design = design)

             if (ask.deff){
                 ## If svytotal has been called on a subset get nobs from the domain index
                 ## else compute it for the whole sample:
                 if (is.null(di <- attr(design, "domain.index"))) {
                     nobs<-NROW(design$cluster)
                    }
                 else {
                     nobs <- length(di)
                    }
                 ###################################
                    # Here z.svyvar instead of svyvar #
                 ###################################
                 if (deff=="replace") {
                     vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2/nobs
                    }
                 else {
                     if (N < nobs) {
                         vsrs<-NA*v
                         warning("Sample size greater than population size: are weights correctly scaled?")
                        }
                     else {
                         vsrs<-z.svyvar(r,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
                        }
                    }
                 # Check for possible anomalies arising from negative calibrated weights (if any)
                 if (inherits(vsrs, "matrix")){
                      # this happens for factor interest variables
                      anom <- ( !is.na(vsrs) & (vsrs < 0) & (row(vsrs)==col(vsrs)) )
                 }
                 else {
                      # this happens for numeric interest variables
                      anom <- (!is.na(vsrs) & (vsrs < 0))
                 }
                 if (any(anom)) {
                      vsrs[anom] <- NA
                      warning("Negative Deff estimate due to negative calibration weights: returning NA")
                 }
                 deffs[i] <- vars[i]/vsrs
                }
            }
        names(vars) <- Rnames
        if (ask.deff){
             names(deffs) <- Rnames
            }
        if (covmat){
            ii<-rep(1:nRatios)
            jj<-rep(1:nRatios)
            allr<-sweep(numerator[,ii]-t(as.vector(rval$ratio)*t(denominator[,jj,drop=FALSE])),
                        2, colSums(denominator[,jj,drop=FALSE]/design$prob),"/")
            vcovmat<-svyrecvar(allr*1/design$prob, design$cluster, design$strata, design$fpc,
                               postStrata=design$postStrata, design = design)
            colnames(vcovmat) <- rownames(vcovmat) <- Rnames
            rval$vcov<-vcovmat
        }
    }

    rval$var<-vars
    attr(rval,"cross") <- cross
    if (ask.deff){
         attr(rval, "deff") <- deffs
        }
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
  }

