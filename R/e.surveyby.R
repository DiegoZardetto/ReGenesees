##
##  tables of statistics.
##
svyby<-function(formula, by, design,...) UseMethod("svyby",design)


svyby.default<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                keep.names=TRUE,verbose=FALSE,vartype=c("se","ci","ci","cv","cvpct","var"), ci.lev=0.95,
                drop.empty.groups=TRUE, covmat=FALSE, return.replicates=FALSE){
######################################################
# MODIFIED things:                                   #
# 1. Unwrap was not correctly working for svyratio   #
# 2. Treat 'formula' argument to cope with svylin    #
#    (which has an expression as first argument)     #
# 3. If formula is not an object of class formula,   #
#    svyby was not correctly working whenever design #
#    happened to be a calibrated object              # 
######################################################

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, model.frame(design), na.action=na.pass)
  else
    byfactors<-as.data.frame(by)

# No replicated designs! Thus no covmat argument...
#  if(covmat || return.replicates){
#    if (!inherits(design,"svyrep.design"))
#      stop("covmat=TRUE not implemented for this design type")
#  }
  if(covmat){
     stop("covmat=TRUE not implemented for domain estimation!")
    }

  ## all combinations that actually occur in this design
  byfactor<-do.call("interaction", byfactors)
  dropped<- weights(design,"sampling")==0
  uniquelevels<-unique(byfactor[!dropped])
  uniques <- match(uniquelevels, byfactor)

####################################
# Following <if> clause has been   #
# ADDED to correctly manage svylin #
####################################
if (!(inherits(formula,"expression") && identical(FUN, svylin))){
  ask.svylin <- FALSE

  ## some people insist on using vectors rather than formulas
  ## so I suppose we should be nice to them
  if (!inherits(formula, "formula")){
      if (NROW(formula)!=length(byfactor))
          stop("'formula' is the wrong length")
      if (!(is.data.frame(formula) ||
            is.matrix(formula) ||
            is.vector(formula))){
          stop("invalid type for 'formula'")
      }
  }
####################################
# End                              #
####################################
}
####################################
# Following <else> clause has been #
# ADDED to correctly manage svylin #
####################################
else {
    ask.svylin <- TRUE
####################################
# End                              #
####################################
}

  if(missing(vartype)) vartype <- "se"
  vartype <- match.arg(vartype,several.ok=TRUE)
  o.vartypes <- eval(formals(sys.function())$vartype)
  nvartype <- which(o.vartypes %in% vartype)
  if(any(is.na(nvartype))) stop("invalid vartype")
  vartype <- o.vartypes[nvartype]                     # Here the original ordering matters (for SE methods etc.) 
  ci.l.tag <- paste("CI.l(", round(100 * ci.lev, 1), "%)", sep = "")
  ci.u.tag <- paste("CI.u(", round(100 * ci.lev, 1), "%)", sep = "")

##########################################
# Following <if else> clause has been    #
# ADDED to correctly manage svyquantile  #
# when a value for ci.lev is passed      #
# (whichever vartype needs ci=TRUE)      #
##########################################
  if (deparse(substitute(FUN))=="svyquantile") { 
    ci.alpha <- 1 - ci.lev
    ci.quant <- TRUE
    } 
  else {
    ci.alpha <- NULL
    ci.quant <- FALSE
    }
  
  if (keep.var){
  ## NOTE: This must happen always!
  ##       Because ReGenesees do NOT allow users to set keep.var = FALSE
      unwrap <-function(x){
        rval<-c(coef(x))
        nvar<-length(rval)
        # variances piece by piece
        se <- c(SE=SE(x))
        ci.l <- confint(x, level=ci.lev)[,1]
        names(ci.l) <- paste(ci.l.tag, names(rval), sep=".")
        ci.u <- confint(x, level=ci.lev)[,2]
        names(ci.u) <- paste(ci.u.tag, names(rval), sep=".")
        cv <- c(CV=cv(x,warn=FALSE))
        cvpct <- c(`CV%`=cv(x,warn=FALSE)*100)
        var <- c(VAR=SE(x)^2)
        # put variances together and keep only those requested
        variances <- c(se, ci.l, ci.u, cv, cvpct, var)[rep((nvartype-1)*(nvar),each=nvar)+(1:nvar)]
        rval<-c(rval, variances)
        if(!is.null(attr(x,"deff")))
          rval<-c(rval,DEff=deff(x))
        rval
      }

      ## In dire need of refactoring (or rewriting)
      ## but it seems to work.
      results<-lapply(uniques,
                      function(i){
                        if(verbose) print(as.character(byfactor[i]))
                        ########################################
                        # Following <else or> clause into <if> #
                        # ADDED to correctly manage calibrated #
                        # design when formula isn't a formula  # 
                        ########################################
                        if (inherits(formula,"formula") || is.calibrated(design))
                          data<-formula
                        #######################################
                        # Following <else if> clause has been #
                        # ADDED to correctly manage svylin    #
                        #######################################
                        else if (ask.svylin)
                          data <- formula
                        else
                          data<-subset(formula, byfactor %in% byfactor[i])
                         # No replicated designs!
                         # if (covmat || return.replicates) {
                         #  FUN(data,
                         #      design[byfactor %in% byfactor[i],],
                         #      deff=deff,alpha=ci.alpha,ci=ci.quant,...,return.replicates=TRUE)
                         # } else {
                          FUN(data,
                              design[byfactor %in% byfactor[i],],
                              deff=deff,alpha=ci.alpha,ci=ci.quant,...)
                         # }
                      })

      ## NOTE 13/05/2021: Below the added if clause, one assumes results to have
      ##                  constant dim for all 'by' subpops. This could be NOT
      ##                  TRUE for svystatB, due to different aliasing in 'by'
      ##                  subpops. The if clause should handle this corner case.
      if (identical(FUN, svylinB)) {
         results2df <- lapply(results, as.data.frame)
         results2df.dims <- lapply(results2df, dim)
         # If dims are not ALL the same, cannot unwrap to eventually produce a
         # a dataframe in output!
         # Therefore: build a list and return it as output!
         if (length(unique(results2df.dims)) > 1) {
             # Get the dataframe of unique by levels
             by.df <- byfactors[uniques, , drop = FALSE]
             # Use it to build meaningful list names
             by.names <- apply(by.df, 1, FUN = function(x) paste(colnames(by.df), x, collapse=":", sep="."))
             names(results2df) <- by.names
             # Order the list the same way rows of ordinary (i.e. matrix like)
             # output objects would be ordered
             results2df <- results2df[order(byfactor[uniques])]
             # Attach minimal attributes to the list
             attr(results2df, "call") <- sys.call()
             # Define a new class for svyby output that are unfortunately lists
             # NOTE: Perhaps one should think to special accessor functions
             #       (coef, SE, ..., summary) for this class!!!!!!!!!!!!!!!!!!!!
             class(results2df) <- c("svyby.list", class(results2df))
             # Return the list
             return(results2df)
            }
        }

      rval<-t(sapply(results, unwrap))
      # No replicated designs!
      # if (covmat || return.replicates) {
      #   replicates<-do.call(cbind,lapply(results,"[[","replicates"))
      #   covmat.mat<-svrVar(replicates,design$scale,design$rscales)
      # }
    } else {
      ## NOTE: This should never happen!
      ##       Because ReGenesees do NOT allow users to set keep.var = FALSE
      unwrap2 <- function(x){
          if(!is.null(attr(x, "deff")))
              c(statistic = unclass(x),
                DEff = deff(x))
          else c(statistic = unclass(x))
      }
      rval<-sapply(uniques,
                   function(i) {
                     if(verbose) print(as.character(byfactor[i]))
                     if (inherits(formula,"formula"))
                           data<-formula
                    #######################################
                    # Following <else if> clause has been #
                    # ADDED to correctly manage svylin    #
                    #######################################
                     else if (ask.svylin)
                            data <- formula
                     else
                       data<-subset(formula, byfactor %in% byfactor[i])
                     unwrap2(FUN(data,
                                 design[byfactor %in% byfactor[i],],
                                 deff=deff,alpha=ci.alpha,ci=ci.quant,...))
                    }
                   )
      if (is.matrix(rval)) rval<-t(rval)
  }

  nr<-NCOL(rval)
# nstats<-nr/(1+ keep.var*length(vartype) +deff) # Here there were problems if deff="replace"
  with.deff <- if ( isTRUE(deff) || (deff=="replace") ) 1 else 0   # modified to deal with deff="replace"
  nstats<-nr/(1+ keep.var*length(nvartype) + with.deff)    # modified to deal with deff="replace"
                                                           # bug fixed: "vartype"->"nvartype" to correctly manage
                                                           # "ci" which counts twice
              
  if (nr>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)

  expand.index<-function(index,reps,x=FALSE){
    ns<-max(index)
    if (x){
      i<-matrix(1:(ns*reps),ncol=reps)
      rval<-t(i[index,])
      
    } else{
      i<-matrix(1:(ns*reps), ncol=reps, nrow=ns, byrow=TRUE)
      rval<- i[index,]
    }
    as.vector(rval)
  }

  if(drop.empty.groups){
      if (keep.names)
          rownames(rval)<-paste(byfactor[uniques])
      rval<-rval[order(byfactor[uniques]),]
      if(covmat){
        i<-expand.index(order(byfactor[uniques]),nstats)
        covmat.mat<-covmat.mat[i,i]
      }
  } else {
      a<-do.call("expand.grid", lapply(byfactors,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor[uniques], levels(byfactor)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor)
      if(covmat){
        tmp<-matrix(ncol=nrow(a)*nstats,nrow=nrow(a)*nstats)
        i<-expand.index(match(byfactor[uniques], levels(byfactor)),nstats,TRUE)
        tmp[i,i]<-covmat.mat
        covmat.mat<-tmp
      }
  }
                  
  attr(rval,"svyby")<-list(margins=1:NCOL(byfactors),nstats=nstats,
                           vars=if(keep.var) length(nvartype) else 0,  # bug fixed: "vartype"->"nvartype" to
                                                                       # correctly manage "ci" which counts twice   
                           deffs=as.logical(with.deff),                # modified to deal with deff="replace"
                           statistic=deparse(substitute(FUN)),
                           variables= names(rval)[-(1:NCOL(byfactors))][1:nstats],
                           vartype=vartype,
                           ci.lev=ci.lev
                           )
  if (!keep.names)
    rownames(rval)<-1:NROW(rval)

  if(covmat)
    attr(rval,"var")<-covmat.mat
  # No replicated designs!
  # if (return.replicates)
  #   attr(rval,"replicates")<-replicates
  attr(rval,"call")<-sys.call()
  class(rval)<-c("svyby","data.frame")
  rval
}

SE.svyby <-function(object,...){
# Modified to correctly manage names etc.
    aa<-attr(object,"svyby")
    if (!aa$vars) stop("Object does not contain variances")
    vartype<-attr(object,"svyby")$vartype
    if (pos<-match("se",vartype,0))
        se <- object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats), drop = FALSE]
    else if (pos<-match("var",vartype,0)) {
        se <- sqrt(object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats), drop = FALSE])
        names(se) <- sub("VAR.", "SE.", x = names(se), fixed = TRUE)
    }
    else if (pos<-match("cv",vartype,0)){
        se <- object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats), drop = FALSE]*coef(object)
        names(se) <- sub("CV.", "SE.", x = names(se), fixed = TRUE)
    }
    else if (pos<-match("cvpct",vartype,0)){
        se <- object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats), drop = FALSE]*coef(object)/100
        names(se) <- sub("CV%.", "SE.", x = names(se), fixed = TRUE)
    }
    else stop("This can't happen")
    return(se)       
}


VAR.svyby <-function(object,...){
  rval <- SE(object,...)^2
  if (!is.null(nam.old <- colnames(rval))) colnames(rval) <- sub("SE.", "VAR.", nam.old, fixed = TRUE)
  # If substitution did not happen, force it  24/02/2016
  # NOTE: This is known to happen only for svystatQ when length(probs) == 1 and by = NULL
  if (identical(nam.old, colnames(rval))) colnames(rval) <- sub("SE", "VAR", colnames(rval), fixed = TRUE)
  rval
}

cv.svyby <-function(object, warn=FALSE, ...){
  rval<-SE(object)/coef(object)
  if (!is.null(nam.old <- names(rval))) names(rval) <- sub("SE.", "CV.", x = nam.old, fixed = TRUE)
  # If substitution did not happen, force it  24/02/2016
  # NOTE: This is known to happen only for svystatQ when length(probs) == 1 and by = NULL
  if (identical(nam.old, names(rval))) names(rval) <- sub("SE", "CV", names(rval), fixed = TRUE)
  if (warn && any(coef(object)<0,na.rm=TRUE)) warning("CV may not be useful for negative statistics")
  rval
}


coef.svyby<-function (object, ...)
{
    aa <- attr(object, "svyby")
    rval <- object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
    }
    names(rval)<-outer(rownames(object),
                       gsub("statistics\\.","",aa$variables), paste, sep=":")
    rval
}

deff.svyby<-function(object,...){
    aa<-attr(object,"svyby")
    if (!aa$deffs) stop("object does not have design effect information")
    object[,max(aa$margins)+aa$nstats*(1+aa$vars)+(1:aa$nstats), drop = FALSE]
}

vcov.svyby<-function(object,...){
  rval<-attr(object,"var")
  if(is.null(rval)){
    # warning("Only diagonal elements of vcov() available")
    se <- SE(object)
    if (is.data.frame(se)) se<-as.vector(as.matrix(se))
    if(length(se)>1)
      rval<-diag(se^2)
    else
      rval<-as.matrix(se^2)
  }
  nms <- names(coef(object))
  dimnames(rval)<-list(nms,nms)
  rval
}

confint.svyquantile<-function(object,parm=NULL,level=NULL,...){
# if (!is.null(level)) stop("need to re-run svyquantile to specify level") # no point, due to new svyby arg ci.lev
  ci<-t(matrix(as.vector(object$CIs),nrow=2))
  colnames(ci)<-dimnames(object$CIs)[[1]]
  rownames(ci)<-outer(dimnames(object$CIs)[[2]],
                      dimnames(object$CIs)[[3]],paste,sep="_")
  if (is.null(parm)) 
    ci
  else 
    ci[parm,,drop=FALSE]
}

extract.ci.by <- function(x) {
# This is for svystatQ
if (!(data.class(x) %in% c("svyby", "svystatQ.by")))
   stop("Input object is not of class svyby")
if ( attr(x,"svyby")$statistic != "svyquantile" )
   stop("Fun argument of svyby was not svyquantile")
cols <- names(x)
cis.pos  <- c(grep("CI.l", cols), grep("CI.u", cols))
x[,cis.pos, drop = FALSE]
}

confint.svyby <- function(object, parm, level = 0.95, ...)
{
    if (missing(level)) level <- attr(object,"svyby")$ci.lev
    # if FUN was svyquantile...
    if ( attr(object,"svyby")$statistic == "svyquantile" ){
        #... check that level is actually equal to the original one
        lev.original <- as.numeric(attr(object,"svyby")$ci.lev)
        if (!isTRUE(all.equal(as.numeric(level), lev.original)))
            stop("Need to re-run svystatQ to specify a confidence level different from the original one (",
                 lev.original,")!")
        return(extract.ci.by(object))
        }
    confint.default(object, level = level, ...)
}
