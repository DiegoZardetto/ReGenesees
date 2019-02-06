##
##  Linearization of multiple regression coefficients.
##  These estimators are SMOOTH ANALYTIC functions of HT
##  or CALIBRATED estimators of totals of numeric variables.
##


linB <- function(model, design, na.rm = TRUE){
###############################################################
# Given a linear regression model formula ('model') computes: #
# 1) The estimate of the regression coefficients              #
# 2) The Woodruff transform of these estimators               #
#                                                             #
# NOTE: Full NA handling.                                     #
# NOTE: Non-positive weights handling: corresponding          #
#       observations give no contribution to estimates and    #
#       sampling errors.                                      #
###############################################################
# Check model:
  # has right class?
if (!inherits(model, "formula"))
    stop("Linear regression model must be specified as a formula")
  # has a response term?
if (length(model)<3)
    stop("Linear regression model must have a response term!")

## Extract all symbols names from model
varnames <- all.vars(model)
## Extract names of design variables
desnames <- names(design$variables)
## Are there unknown symbols in model?
unknown.vars <- varnames[!(varnames %in% desnames)]
if (length(unknown.vars) > 0)
    stop("Unknown variables in regression model formula: ", paste(unknown.vars, collapse = ", "))

## Check model variable types
typetest <- sapply(varnames, function(y) is.numeric(design$variables[, y]) |
                                         is.factor(design$variables[, y])   )
if (!all(typetest))
    stop("Variables appearing in regression model must be numeric or factor")

## Get basic structures: 1. model frame
##                       2. model matrix
##                       3. model response
##                       4. weights
# 1. Model frame
mf <- model.frame(model, data = design$variables, na.action = na.pass)
# 2. Model matrix
mm <- model.matrix(model, data = mf, na.action = na.pass)
# 3. Model response
  # get
resp <- model.response(mf)
  # check type
if (!is.numeric(resp))
    stop("Regression model response term must be numeric!")
# 4. Sampling unit weights
ww <- weights(design)

# Check for non-positive weights: if any, "drop" corresponding cases
ww.ko <- (ww <= 0)
N.ww.ko <- sum(ww.ko)
if (N.ww.ko > 0) {
     warning("Detected ", N.ww.ko,
             " observations with weight <= 0: they will be ignored!", sep="")
     ww[ww.ko] <- 0
    }

## Missing values treatment
 # Check for NAs
nas <- !complete.cases(mf)
N.nas <- sum(nas)
# Store the vector of indices of KO cases (i.e. with w<=0 or NAs)
# it will be an attribute of the result
ko.obs <- (ww.ko | nas)
if (N.nas > 0){
     var.na <- names(mf)[sapply(mf, function(x) any(is.na(x)))]
     if (!na.rm){
         stop("Missing values in: ", paste(var.na, collapse = ", "),
              "\n(you may want to specify na.rm=TRUE)")
        }
      else {
         ## Due to NAs must "reduce" basic structures 1.-4. (indeed foreign
         ## language routines, e.g. qr, cannot cope with NAs)
         warning("Missing values in: ", paste(var.na, collapse = ", "),
                 "; corresponding observations (", N.nas, ") will be ignored!")
         mf.r <- mf[!nas, , drop=FALSE]
         mm.r <- mm[!nas, , drop=FALSE]
         resp.r <- resp[!nas]
         ww.r <- ww[!nas]

         # Identify 0 weights among those of non-NAs observations
         ww.r.ko <- (ww.r <= 0)
         # Take weights sqrt: now it's safe (no more ww < 0)
         whalf.r <- sqrt(ww.r)
         # Compute weighted model matrix and its QR decomposition
         # now it'safe (no more NAs)
         mm.whalf.r <- mm.r * whalf.r
         tqr.r <- qr(mm.whalf.r)
 
         # Compute regression coefficients
         Beta <- qr.coef(tqr.r, resp.r * whalf.r)

         # Collinearity check: redundant variables will be aliased
         alias <- any(aliased <- is.na(Beta))
         if (alias) {
             var.aliased <- names(Beta)[aliased]
             warning("Variables ", paste(var.aliased, collapse = ", "), " have been aliased due to collinearity!")
             # Reduce the weighted model matrix and re-compute its QR decomposition
             mm.r <- mm.r[, !aliased, drop=FALSE]
             mm.whalf.r <- mm.r * whalf.r
             tqr.r <- qr(mm.whalf.r)
             # Re-compute regression coefficients
             Beta <- qr.coef(tqr.r, resp.r * whalf.r)
            }

         # Compute residuals
         ee.r <- qr.resid(tqr.r, resp.r * whalf.r) / whalf.r
         # Take care of NaN arising from 0 whalf.r (if any)
         ee.r[ww.r.ko] <- 0

         # Compute (generalized, if needed) inverse of t(X)%*%W%*%X
         T.r <- crossprod(mm.whalf.r)
         Tm1.r <- try(solve(T.r), silent = TRUE)

         # Residual collinearity after aliasing (this should happen only
         # for numerical reasons...)
         if (collin.res <- inherits(Tm1.r, "try-error")) {
             warning("Model matrix is singular: switching to Moore-Penrose generalized inverse.")
             ## No longer needed: ReGenesees now IMPORTS MASS
             # require(MASS)
             Tm1.r <- ginv(T.r)
            }

         # Compute Beta Woodruff transform (see my notes on paper)
         z.Beta.r <- tcrossprod(mm.r * as.numeric(ee.r), Tm1.r)

         ## Re-expand reduced z.Beta.r to original dimensions
          # prepare a matrix with right dim 
          z.Beta <- mm
          # set all entries to 0
          z.Beta[, ] <- 0
          # Copy meaningful values in the right places (the non-NA ones)
          z.Beta[!nas, ] <- z.Beta.r
         # Set appropriate column names
         dimnames(z.Beta)[[2]] <- names(Beta)

         ## Build result
         out <- list(Beta = Beta, z.Beta = z.Beta)
        }
    }
else{
     ## No NAs case: everything is ok.
     # Take weights sqrt: now it's safe (no more ww<0)
     whalf <- sqrt(ww)
     # Compute weighted model matrix and its QR decomposition
     mm.whalf <- mm * whalf
     tqr <- qr(mm.whalf)

     # Compute regression coefficients
     Beta <- qr.coef(tqr, resp * whalf)

     # Collinearity check: redundant variables will be aliased
     alias <- any(aliased <- is.na(Beta))
     if (alias) {
         var.aliased <- names(Beta)[aliased]
         warning("Variables ", paste(var.aliased, collapse = ", "), " have been aliased due to collinearity!")
         # Reduce the weighted model matrix and re-compute its QR decomposition
         mm <- mm[, !aliased, drop=FALSE]
         mm.whalf <- mm * whalf
         tqr <- qr(mm.whalf)
         # Re-compute regression coefficients
         Beta <- qr.coef(tqr, resp * whalf)
        }

     # Compute residuals
     ee <- qr.resid(tqr, resp * whalf) / whalf
     # Take care of NaN arising from 0 whalf (i.e. non-positive weights, if any)
     ee[ww.ko] <- 0

     # Compute (generalized, if needed) inverse of t(X)%*%W%*%X
     T <- crossprod(mm.whalf)
     Tm1 <- try(solve(T), silent = TRUE)

     # Residual collinearity after aliasing (this should happen only
     # for numerical reasons...)
     if (collin.res <- inherits(Tm1, "try-error")) {
         warning("Model matrix is singular: switching to Moore-Penrose generalized inverse.")
         ## No longer needed: ReGenesees now IMPORTS MASS
         # require(MASS)
         Tm1 <- ginv(T)
        }

     # Compute Beta Woodruff transform
     z.Beta <- tcrossprod(mm * as.numeric(ee), Tm1)
     # Set appropriate column names
     dimnames(z.Beta)[[2]] <- names(Beta)

     ## Build result
     out <- list(Beta = Beta, z.Beta = z.Beta)
    }
## Return value
# Save ko observation indices (if any)
if (any(ko.obs)) {
     attr(out, "ko.obs") <- which(ko.obs)
    }
return(out)
}


svylinB <- function(model, design, na.rm=FALSE, ...){
  UseMethod("svylinB", design)
}

svylinB.survey.design2 <- function(model, design, na.rm=FALSE, deff=FALSE, ...){
B <- linB(model = model, design = design, na.rm = na.rm)
Beta <- B$Beta
z.Beta <- B$z.Beta
class(Beta)<- c("svylinB", "svystat")

## ko cases treatment: "drop" them from design (actually put weight to 0)
if (!is.null(ko <- attr(B,"ko.obs"))){
    design$prob[ko] <- Inf
    }

attr(Beta, "var") <- v <- svyrecvar(z.Beta/design$prob, design$cluster,
                                    design$strata, design$fpc,
                                    postStrata=design$postStrata, design = design)
dimnames(attr(Beta, "var"))<-list(names(Beta),names(Beta))
attr(Beta,"statistic")<-"regcoef"
if (is.character(deff) || deff){
    ###################################
    # Here z.svyvar instead of svyvar #
    ###################################
    N<-sum(1/design$prob)
    nobs<-NROW(design$cluster)
    if (deff=="replace")
        vsrs<-z.svyvar(z.Beta,design,na.rm=na.rm)*sum(weights(design))^2/nobs
    else
        vsrs<-z.svyvar(z.Beta,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
    attr(Beta, "deff")<-v/vsrs
    dimnames(attr(Beta, "deff")) <- list(names(Beta),names(Beta))
    }
return(Beta)
}

svylinB.cal.analytic <- function(model, design, na.rm=FALSE, deff=FALSE, ...){
B <- linB(model = model, design = design, na.rm = na.rm)
Beta <- B$Beta
z.Beta <- B$z.Beta
class(Beta)<- c("svylinB", "svystat")

## ko cases treatment: "drop" them from design (actually put weight to 0)
if (!is.null(ko <- attr(B,"ko.obs"))){
    design$prob[ko] <- Inf
    }

attr(Beta, "var") <- v <- svyrecvar(z.Beta/design$prob, design$cluster,
                                    design$strata, design$fpc,
                                    postStrata=design$postStrata, design = design)
dimnames(attr(Beta, "var"))<-list(names(Beta),names(Beta))
attr(Beta,"statistic")<-"regcoef"
if (is.character(deff) || deff){
    N<-sum(1/design$prob)
    ## If svylinB has been called on a subset (actually this CANNOT be the case,
    ## as no 'by' argument is provided) get nobs from the domain index
    ## else compute it for the whole sample:
    if (is.null(di <- attr(design, "domain.index"))) nobs<-NROW(design$cluster) else nobs <- length(di)
    ###################################
    # Here z.svyvar instead of svyvar #
    ###################################
    if (deff=="replace") {
        vsrs<-z.svyvar(z.Beta,design,na.rm=na.rm)*sum(weights(design))^2/nobs
    }
    else {
        if (N < nobs) {
            vsrs<-NA*v
            warning("Sample size greater than population size: are weights correctly scaled?")
        }
        else {
            vsrs<-z.svyvar(z.Beta,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
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
    attr(Beta, "deff")<-v/vsrs
    dimnames(attr(Beta, "deff")) <- list(names(Beta),names(Beta))
}
return(Beta)
}
