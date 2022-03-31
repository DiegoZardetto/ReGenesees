`svyDelta` <- function(expr, design1, design2, by = NULL,
                       des.INDEP = FALSE, rho.STRAT = c("Full", "noJump", "noStrat"),
                       vartype = c("se", "cv", "cvpct", "var"), conf.int = FALSE, conf.lev = 0.95) {

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design1.expr <- if (directly) substitute(design1)
design2.expr <- if (directly) substitute(design2)

## Only most essential check on input type (more in inner functions...)
# Designs
if (!inherits(design1, "analytic")) {
     stop("Object 'design1' must inherit from class analytic")
    }
if (!inherits(design2, "analytic")) {
     stop("Object 'design2' must inherit from class analytic")
    }
# Estimator expression
if (!inherits(expr, "expression")){
    stop("The estimator of the measure of change must be specified as an expression")
    }
# Vartype and confidence intervals
if (missing(vartype)) vartype <- "se"
vartype <- match.arg(vartype, several.ok = TRUE)
vartype <- unique(vartype)
if (conf.int) {
    vartype <- c(vartype, "ci")
    l.conf.tag <- paste("CI.l(", round(100 * conf.lev, 1), "%)", sep = "")
    u.conf.tag <- paste("CI.u(", round(100 * conf.lev, 1), "%)", sep = "")
    if (!is.numeric(conf.lev))
        stop("conf.lev must be numeric")
    if (conf.lev < 0 || conf.lev > 1) 
        stop("conf.lev must fall inside [0,1]")
    }
# By variables
if (!is.null(by)) {
     if (!inherits(by, "formula")) {
         stop("If specified, 'by' must be supplied as a formula")
        }
     by.vars <- all.vars(by)
     na.Fail(design1$variables, by.vars)
     na.Fail(design2$variables, by.vars)
#    # The type check below is commented because usually other ReGenesees
#    # estimation functions do *not* perform it. Consider de-commenting, if needed
#      typetest1 <- sapply(by.vars, function(v) is.factor(design1$variables[, v]))
#      if (!all(typetest1)) 
#          stop("By variables must be factors")
#      typetest2 <- sapply(by.vars, function(v) is.factor(design2$variables[, v]))
#      if (!all(typetest2)) 
#          stop("By variables must be factors")
    }

## Stratification
# Approach to handle strata in estimation of rho
rho.STRAT <- match.arg(rho.STRAT)
rho.UNSTRAT <- isTRUE(rho.STRAT == "noStrat")
rho.NOSTRATJUMP <- isTRUE(rho.STRAT == "noJump")

has.strata1 <- design1$has.strata
has.strata2 <- design2$has.strata
if (has.strata1 && has.strata2){
     has.strata <- TRUE
    } else {
     if (!has.strata1 && !has.strata2) {
         has.strata <- FALSE
        } else {
         if (!isTRUE(rho.UNSTRAT) && !des.INDEP) {
             stop("Objects 'design1' and 'design2' must be either both stratified or both unstratified")
            }
        }
    }

# If samples are declared to have been sampled independently, no rho estimation
# needed and thus can set has.strata arbitrarily (will not be used):
if (des.INDEP) {
     has.strata <- FALSE
    }

## In case of stratification, should we use it to compute rho? And how?
# if (has.strata1 || has.strata2) {
#      if (isTRUE(rho.UNSTRAT) && !des.INDEP) {
#          # Not sure wheter a warning is appropriate for a user tunable option...
#          warning("Stratification of input design objects disregarded in estimating sampling correlations")
#         }
#      if (isTRUE(rho.NOSTRATJUMP) && !des.INDEP) {
#          # Not sure wheter a warning is appropriate for a user tunable option...
#          warning("Units that changed stratum from 'design1' to 'design2' (if any) disregarded in estimating sampling correlations")
#         }
#     }
use.strata <- has.strata && !isTRUE(rho.UNSTRAT)

## Clustering
is.element1 <- is.element(design1)
is.element2 <- is.element(design2)
if (is.element1 && is.element2){
     is.element <- TRUE
    } else {
     if (!is.element1 && !is.element2) {
         is.element <- FALSE
        } else {
         if (!des.INDEP) {
             stop("Objects 'design1' and 'design2' must be either both element sampling designs or both cluster sampling designs")
            }
        }
    }

# If samples are declared to have been sampled independently, no rho estimation
# needed and thus can set is.element arbitrarily (will not be used):
if (des.INDEP) {
     is.element <- TRUE
    }

## Domain estimation
if (!is.null(by)) {
     stat <- svyby.svydelta(expr = expr, by = by, design1 = design1, design2 = design2,
                   has.strata = use.strata, is.element = is.element, are.indep = des.INDEP, no.strat.jump = rho.NOSTRATJUMP,
                   keep.names = TRUE, verbose = FALSE,
                   vartype = vartype, ci.lev = conf.lev,
                   drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm)
     class(stat) <- c("svyDelta.by", class(stat))
     attr(stat,"design1") <- design1.expr
     attr(stat,"design2") <- design2.expr
     attr(stat, "call") <- sys.call()
     return(stat)
    } else {
     design1 <- des.addvars(design1, FaKe.by = factor(1))
     design2 <- des.addvars(design2, FaKe.by = factor(1))
     stat <- svyby.svydelta(expr = expr, by = ~FaKe.by, design1 = design1, design2 = design2,
                   has.strata = use.strata, is.element = is.element, are.indep = des.INDEP, no.strat.jump = rho.NOSTRATJUMP,
                   keep.names = TRUE, verbose = FALSE,
                   vartype = vartype, ci.lev = conf.lev,
                   drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm)

     # Get the details attribute
     details <- attr(stat, "details")
     # Now change stat's format (recall stat is a dataframe)    
     # Drop fake 'by' column
     stat.f <- stat[, -which(names(stat) == "FaKe.by"), drop=FALSE]
     # Count how many stats (numbers) inside
     nstats <- ncol(stat.f)
     # Count how many variability measures inside (recall "ci" counts twice but already in attr(stat, "svyby")$vars)
     nvariances <- attr(stat, "svyby")$vars + attr(stat, "svyby")$deffs
     # Count how many variables inside
     nvariables <- attr(stat, "svyby")$nstats
     # Consistency check
       # DEBUG 12/06/2020: BUG: if clauses have to use !isTRUE(all.equal()). Fixed
     if (!isTRUE(all.equal(nstats, nvariables*(1 + nvariances)))) stop("Format modification failed")
     # Store properly the relevant stats into a matrix
     stat.mat <- matrix(nrow = nvariables, ncol = (1 + nvariances))
     stat.mat[, ] <- as.matrix(stat.f[, ])
     # Rownames are the variables names
     rownames(stat.mat) <- attr(stat, "svyby")$variables
     # Colnames are Estimator type + Variability measures
     o.vartypes <- c("se","ci","ci","cv","cvpct","var")
     varia.name  <- o.vartypes[o.vartypes %in% attr(stat, "svyby")$vartype]
     colnames(stat.mat) <- c("Delta", varia.name, if (attr(stat, "svyby")$deffs) "DEff")
     # Cast into a data.frame
     stat.df <- as.data.frame(stat.mat)
     attr(stat.df, "origin") <- stat
     class(stat.df) <- c("svyDelta", class(stat.df))
     # Better column names
     tmp.names <- names(stat.df)
     tmp.names[match("se", tmp.names)]    <- "SE"
     tmp.names[match("cv", tmp.names)]    <- "CV"
     tmp.names[match("cvpct", tmp.names)] <- "CV%"
     tmp.names[match("var", tmp.names)]   <- "VAR"
     if (conf.int) {
         tmp.names[match("ci", tmp.names)] <- l.conf.tag
         tmp.names[match("ci", tmp.names)] <- u.conf.tag
     }
     names(stat.df) <- tmp.names

     # Attach the details data frame
     attr(stat.df, "details") <- details

     attr(stat.df,"design1") <- design1.expr
     attr(stat.df,"design2") <- design2.expr
     attr(stat.df, "call") <- sys.call()
     stat.df
     }
}

## Utility identify element vs cluster sampling designs
`is.element` <- function(design) {
     n <- NROW(design$cluster)
     un <- length(unique(design$cluster[, 1]))
     if (n == un) {
         ans <- TRUE
        } else {
         ans <- FALSE
        }
     return(ans)
    }

## Accessor functions for objects created by svyDelta
`coef.svyDelta` <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

`vcov.svyDelta` <- function(object, ...){
    by.object <- attr(object, "origin")
    rval <- attr(by.object,"var")
    if(is.null(rval)){
       # warning("Only diagonal elements of vcov() available")
       se <- SE(by.object)
    if (is.data.frame(se)) se <- as.vector(as.matrix(se))
    if (length(se) > 1)
        rval<-diag(se^2)
    else
        rval<-as.matrix(se^2)
  }
  nms <- names(coef(object))
  dimnames(rval)<-list(nms,nms)
  rval
}

`SE.svyDelta` <- function(object, ...) SE(attr(object, "origin"))

`VAR.svyDelta` <- function(object, ...) VAR(attr(object, "origin"))

`cv.svyDelta` <- function(object, ...) cv(attr(object, "origin"))

`deff.svyDelta` <- function(object, ...) deff(attr(object, "origin"))

`confint.svyDelta` <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}


# Accessor function to the details data frame of classes
`details` <- function(object, print.call = TRUE, ...) {
  if (! ( inherits(object, "svyDelta") || inherits(object, "svyDelta.by") ) ) {
     stop("Input object was not created by function svyDelta") 
    }
  print(object)
  if (isTRUE(print.call)) {
     cat("\nCall:\n")
     prcall(attr(object, "call"))
    }
  cat("\n")
  cat("Details:\n")
  cat("\n")
  details <- attr(object, "details")
  if ("msg" %in% colnames(details)) {
     details <- unique(details[["msg"]])
     cat(details, "\n", sep="")
    } else {
     print(details)
    }
  return(invisible(details))
}


`svyby.svydelta` <- function(expr, by, design1, design2, has.strata, is.element, are.indep, no.strat.jump, ..., 
                             keep.names=TRUE, verbose=FALSE,
                             vartype=c("se","ci","ci","cv","cvpct","var"), ci.lev=0.95,
                             drop.empty.groups=TRUE, covmat=FALSE){
######################################################################
# Workhorse function for domain estimation, to serve svyDelta needs. #
######################################################################
  if(covmat){
     stop("covmat=TRUE not implemented for domain estimation!")
    }

     ## Order combinations as svyby would do
     ## design1
     byfactors1 <- model.frame(by, model.frame(design1), na.action = na.pass)
     ## All combinations that actually occur in this design
     byfactor1<-do.call("interaction", byfactors1)
     dropped1<- weights(design1,"sampling")==0
     uniquelevels1<-unique(byfactor1[!dropped1])
     ndomlev1 <- length(uniquelevels1)
     ## design2
     byfactors2 <- model.frame(by, model.frame(design2), na.action = na.pass)
     ## All combinations that actually occur in this design
     byfactor2<-do.call("interaction", byfactors2)
     dropped2<- weights(design2,"sampling")==0
     uniquelevels2<-unique(byfactor2[!dropped2])
     ndomlev2 <- length(uniquelevels2)     
     ## Restrict to by-levels that are common to design1 and design2
     uniquelevels <- intersect(uniquelevels1, uniquelevels2)
     ndomlev12 <- length(uniquelevels)
     if (ndomlev12 == 0) {
         stop("Objects 'design1' and 'design2' have no 'by' domains in common")
        } else {
         if ( ndomlev12 < max(ndomlev1, ndomlev2) ) {
             warning("Objects 'design1' and 'design2' have only a subset of 'by' domains in common")
            }
        } 
     ## All common-domain sub-designs of design1
     uniques1 <- match(uniquelevels, byfactor1)
#     designs1 <- lapply(uniques1, function(i) design1[byfactor1 %in% byfactor1[i], ])
     ## All common-domain sub-designs of design1
     uniques2 <- match(uniquelevels, byfactor2)
#     designs2 <- lapply(uniques2, function(i) design2[byfactor2 %in% byfactor2[i], ])

  if(missing(vartype)) vartype <- "se"
  vartype <- match.arg(vartype,several.ok=TRUE)
  o.vartypes <- eval(formals(sys.function())$vartype)
  nvartype <- which(o.vartypes %in% vartype)
  if(any(is.na(nvartype))) stop("invalid vartype")
  vartype <- o.vartypes[nvartype]                     # Here the original ordering matters (for SE methods etc.) 
  ci.l.tag <- paste("CI.l(", round(100 * ci.lev, 1), "%)", sep = "")
  ci.u.tag <- paste("CI.u(", round(100 * ci.lev, 1), "%)", sep = "")

  
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
        rval
    }

      ## In dire need of refactoring (or rewriting)
      ## but it seems to work.
      results <- mapply(uniques1, uniques2,
                        FUN = function(i, j){
                                if(verbose) print(as.character(byfactor1[i]))
                                svydelta(expr, design1[byfactor1 %in% byfactor1[i], ], design2[byfactor2 %in% byfactor2[j], ],
                                         has.strata = has.strata, is.element = is.element, are.indep = are.indep, no.strat.jump = no.strat.jump, ...
                                        )
                            }, SIMPLIFY = FALSE
                    )

      ## Get estimates and errors
      rval<-t(sapply(results, unwrap))

      ## Get details data frames and rbind them
      details <- Reduce(rbind, lapply(results, function(el) attr(el, "details")))
      ## Add domain columns
      if (deparse(by) != "~FaKe.by") { 
             details <- cbind(byfactors1[uniques1,,drop=FALSE], details)
            }

  nr<-NCOL(rval)
  nstats<-nr/(1+ length(nvartype))  # bug fixed: "vartype"->"nvartype" to correctly manage
                                    # "ci" which counts twice

  if (nr>1)
    rval<-cbind(byfactors1[uniques1,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors1[uniques1,,drop=FALSE], statistic=rval)

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
      if (keep.names) {
          rownames(rval) <- paste(byfactor1[uniques1])
          if (deparse(by) != "~FaKe.by") { 
             rownames(details) <- paste(byfactor1[uniques1])
            }
        }
      # reorder rows
      rval<-rval[order(byfactor1[uniques1]),]
      details <- details[order(byfactor1[uniques1]), , drop = FALSE]

      if(covmat){
        i<-expand.index(order(byfactor1[uniques1]),nstats)
        covmat.mat<-covmat.mat[i,i]
      }
  } else { # never the case in ReGenesees, in practice
      a<-do.call("expand.grid", lapply(byfactors1,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor1[uniques1], levels(byfactor1)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor1)
      if(covmat){
        tmp<-matrix(ncol=nrow(a)*nstats,nrow=nrow(a)*nstats)
        i<-expand.index(match(byfactor1[uniques1], levels(byfactor1)),nstats,TRUE)
        tmp[i,i]<-covmat.mat
        covmat.mat<-tmp
      }
  }
                  
  attr(rval,"svyby")<-list(margins=1:NCOL(byfactors1),nstats=nstats,
                           vars = length(nvartype),                    # bug fixed: "vartype"->"nvartype" to
                                                                       # correctly manage "ci" which counts twice   
                           deffs=FALSE,                                # modified to deal with deff="replace"
                           statistic="svydelta",
                           variables= names(rval)[-(1:NCOL(byfactors1))][1:nstats],
                           vartype=vartype,
                           ci.lev=ci.lev
                           )

  if (!keep.names) {
     rownames(rval)<-1:NROW(rval)
     rownames(details)<-1:NROW(details)
    }

  # Details data frame
  attr(rval, "details") <- details

  if(covmat)
    attr(rval,"var")<-covmat.mat
  attr(rval,"call")<-sys.call()
  class(rval)<-c("svyby","data.frame")
  rval
}
