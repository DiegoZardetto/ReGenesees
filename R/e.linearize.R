##
##  Linearization of a complex estimator (assumed to be a SMOOTH ANALYTIC function
##  of HT or CALIBRATED estimators of totals of numeric variables).
##


linearize <- function(expr, design, na.rm = FALSE, ...){
#################################################################################
# Given a complex estimator, expressed as an analytic function of HT estimators #
# (possibly including some known fixed parameter, e.g. a denominator total for  #
# a Ratio Estimator of a Total) computes:                                       #
# 1) The estimate of the complex estimator                                      #
# 2) The Woodruff transform of the estimator                                    #
#                                                                               #
# NOTE: Must write some documentation for the function's arguments              #
#                                                                               #
# NOTE: Code has been slighly revised and simplified: '...' argument is no      #
#       longer used to pass known fixed parameters. These are instead simply    #
#       searched in the standard way R searches binding to unbound symbols      #
#       inside function bodies. 29/02/2015                                      #
#       I decided to let the '...' argument survive in unexported function      #
#       linearize for possible future extensions. 29/02/2016                    #
#################################################################################
if (!inherits(expr, "expression")){
    stop("Estimator to be linearized must be specified as an expression")
    }
if (length(expr)>1){
    stop("Can specify just a single Complex Estimator at a time!")
    }
## Extract all symbols names from expr (including parameters to be passed by '...')
varnames <- all.vars(expr)
## Check if expr refers to 'ones' convenience variable
if ("ones" %in% varnames) design$variables[["ones"]] <- 1

## Extract parameters names from '...'
# other <- list(...)
# othernames <- names(other)

## Extract names of design variables
desnames <- names(design$variables)

## Are there unknown symbols in expr?
# unknown.vars <- varnames[!(varnames %in% c(desnames, othernames))]
# if (length(unknown.vars) > 0)
   # stop("Unknown variables in estimator expression: ", paste(unknown.vars, collapse = ", "))
## Are there conflicting symbols in expr, i.e. symbols contained both in design and '...'?
# conflict.vars <- desnames[desnames %in% othernames]
# if (length(conflict.vars) > 0) {
    # warning("Symbols <",
            # paste(conflict.vars, collapse = ", "),
            # "> refer both to survey variables and passed parameters: will be treated as survey variables",
            # immediate. = TRUE)
    # othernames <- othernames[othernames!=conflict.vars]
# }

## Extract from expr the "internal variables" i.e. symbols referencing design variables
# in.vars <- varnames[!(varnames %in% othernames)]
in.vars <- varnames[varnames %in% desnames]
##  Check if there are actually "internal variables"
if (length(in.vars) == 0) {
    stop("Specified estimator isn't a function of observed variables")
}

## Check "internal variables" type
typetest <- sapply(in.vars, function(y) is.numeric(design$variables[, y]))
if (!all(typetest)) 
    stop("Estimator formula cannot contain non numeric variables")

# Check for missing values in interest variables (one at a time, since only a
# single estimator is involved)
for (var in in.vars){
     NA.estvars(design = design, estvars = var, na.rm = na.rm)
    }

## Extract from expr the "parameters" i.e. symbols NOT referencing design
## variables
parameters <- varnames[!(varnames %in% desnames)]

## Parameters must exist and be numeric scalars, check it out
is.scalar <- function(x) {
     x <- get(x)
     is.vector(x, mode = "numeric") && identical(length(x), 1L)
    }
scalartest <- sapply(parameters, is.scalar)
if (!all(scalartest)) 
    stop("Estimator formula can contain only scalar numeric parameters")

## Build a convenience formula with "internal variables"...
in.vars.formula <- as.formula(paste("~",paste(in.vars, collapse="+"),sep=""))
## ... and its model frame
mf <- model.frame(in.vars.formula, design$variables, na.action = na.pass)
## Convert model frame into matrix
x <- as.matrix(mf)

## Missing values treatment
nas <- rowSums(is.na(x))
if (na.rm && sum(nas)>0){
    design.old <- design
    design<-design[nas==0,]
    # If domain has some non-NA values, use them for estimation:
    if (length(design$prob) > 0) {
        if (length(nas)>length(design$prob))  # i.e. design was not cal
            x<-x[nas==0,,drop=FALSE]
        else                                  # i.e. design was cal
            x[nas>0,]<-0
        }
    # If domain has only NAs, cannot do anything:
    else {
         na.rm <- FALSE
         design <- design.old
        }
    }

## Compute estimators of totals of "internal variables"
estimators <- as.list(coef(z.svytotal(x , design = design)))
# env <- c(estimators, other)
env <- estimators
## Symbolically compute both expr and its first order derivatives w.r.t the "internal variables"
f.plus.grad <- deriv(expr = expr, namevec = in.vars, function.arg = FALSE)
## Evaluate symbolic expressions
f.plus.grad.val <- eval(f.plus.grad, envir = env)
 # 1. Gradient
grad.val <- attr(f.plus.grad.val, "gradient")
 # 2. Linearized variable
z <- x%*%t(grad.val)
 # 3. Estimate
attr(z,"estimate") <- as.numeric(f.plus.grad.val)
## If call had na.rm = TRUE and there were some non-Na values,
## store NA indexes to later subset the design object
## (see 'svylin' function below):
attr(z,"nas") <- if (na.rm && sum(nas)>0) nas
## Return result
z
}


svylin <- function(expr, design, na.rm = FALSE, ...){
# .svycheck(design)
  UseMethod("svylin",design)
}

svylin.survey.design2 <- function(expr, design, na.rm = FALSE, deff = FALSE, ...){
z <- linearize(expr = expr, design = design, na.rm = na.rm, ...)
estimate <- attr(z,"estimate")
names(estimate) <- as.character(expr)
class(estimate)<- c("svylin", "svystat")

## Missing values treatment
if (na.rm && !is.null(nas <- attr(z,"nas"))){
    design<-design[nas==0,]
    }

attr(estimate, "var") <- v <- svyrecvar(z/design$prob, design$cluster,
                                        design$strata, design$fpc,
                                        postStrata=design$postStrata, design = design)
dimnames(attr(estimate, "var"))<-list(names(estimate),names(estimate))
attr(estimate,"statistic")<-"complex"
if (is.character(deff) || deff){
    ###################################
    # Here z.svyvar instead of svyvar #
    ###################################
    N<-sum(1/design$prob)
    nobs<-NROW(design$cluster)
    if (deff=="replace")
        vsrs<-z.svyvar(z,design,na.rm=na.rm)*sum(weights(design))^2/nobs
    else
        vsrs<-z.svyvar(z,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
    attr(estimate, "deff")<-v/vsrs
    dimnames(attr(estimate, "deff")) <- list(names(estimate),names(estimate))
    }
return(estimate)
}

svylin.cal.analytic <- function(expr, design, na.rm = FALSE, deff = FALSE, ...){
z <- linearize(expr = expr, design = design, na.rm = na.rm, ...)
estimate <- attr(z,"estimate")
names(estimate) <- as.character(expr)
class(estimate)<- c("svylin", "svystat")

## Missing values treatment
if (na.rm && !is.null(nas <- attr(z,"nas"))){
    design<-design[nas==0,]
    }

attr(estimate, "var") <- v <- svyrecvar(z/design$prob, design$cluster,
                                        design$strata, design$fpc,
                                        postStrata=design$postStrata, design = design)
dimnames(attr(estimate, "var"))<-list(names(estimate),names(estimate))
attr(estimate,"statistic")<-"complex"
if (is.character(deff) || deff){
    N<-sum(1/design$prob)
    ## If svylin has been called on a subset get nobs from the domain index
    ## else compute it for the whole sample:
    if (is.null(di <- attr(design, "domain.index"))) nobs<-NROW(design$cluster) else nobs <- length(di)
    ###################################
    # Here z.svyvar instead of svyvar #
    ###################################
    if (deff=="replace") {
        vsrs<-z.svyvar(z,design,na.rm=na.rm)*sum(weights(design))^2/nobs
    }
    else {
        if (N < nobs) {
            vsrs<-NA*v
            warning("Sample size greater than population size: are weights correctly scaled?")
        }
        else {
            vsrs<-z.svyvar(z,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
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
    attr(estimate, "deff")<-v/vsrs
    dimnames(attr(estimate, "deff")) <- list(names(estimate),names(estimate))
}
return(estimate)
}
