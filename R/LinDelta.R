##
##  Linearization of the estimator of a Measure of Change computed using two
##  (in general not independent) samples s1, s2.
##
##  The estimator is assumed to be a SMOOTH ANALYTIC function of HT or
##  CALIBRATED estimators of totals of numeric variables derived from s1 and s2.
##

# References
# - [Berger, Priam 16]
# Berger, Y. G., Priam, R. (2016). A simple variance estimator of change for
# rotating repeated surveys: an application to the European Union Statistics on
# Income and Living Conditions household surveys. Journal of the Royal
# Statistical Society: Series A (Statistics in Society), 179(1), 251-272.


`LinDelta` <- function(expr, design1, design2, ...){
#################################################################################
# Given a Measure of Change, Delta, expressed as an analytic function of HT or  #
# CALIBRATION estimators of totals of variables coming from two (in general NOT #
# INDEPENDENT) samples, design1 and design2, computes:                          #
#                                                                               #
# 1) The estimate of the Measure of Change, Delta                               #
#                                                                               #
# 2) The Woodruff transforms (w.r.t. design1 and design2) of the estimator      #
#                                                                               #
# 3) The variance of the estimator of the total of the Woodruff transforms      #
#    w.r.t. to their reference designs (design1 and design2)                    #
#                                                                               #
# NOTE: Any measure of change Delta will be reduced, upon linearization, to the #
#       standard difference form:                                               #
#                                                                               #
#       LDelta = Ly.2 - Ly.1                                                    #
#                                                                               #
#       where Ly.1 and Ly.2 are the Woodruff transforms of the Measure of       #
#       Change w.r.t. design1 and design2, with suitable signs.                 #
#                                                                               #
# NOTE: In case any of the designs is calibrated, its Woodruff transform will   #
#       include the appropriate residual of the variable w.r.t. the calibration #
#       model.                                                                  #
#################################################################################
if (!inherits(expr, "expression")){
    stop("Estimator to be linearized must be specified as an expression")
    }
if (length(expr)>1){
    stop("Can specify just a single Complex Estimator at a time!")
    }
## Extract all symbols names from expr (including parameters to be passed by '...')
varnames <- all.vars(expr)

## Split them according to the design object they refer to
## design1
varnames1 <- varnames[endsWith(varnames, ".1")]
varnames1 <- substr(varnames1, 1, nchar(varnames1) - 2)
## design2
varnames2 <- varnames[endsWith(varnames, ".2")]
varnames2 <- substr(varnames2, 1, nchar(varnames2) - 2)

## Check if expr refers to 'ones' convenience variable
## design1
if ("ones" %in% varnames1) design1$variables[["ones"]] <- 1
## design2
if ("ones" %in% varnames2) design2$variables[["ones"]] <- 1

## Extract names of design variables
## design1
desnames1 <- names(design1$variables)
## design2
desnames2 <- names(design2$variables)

## Extract from expr the "internal variables" i.e. symbols referencing design variables
## design1
in.vars1 <- varnames1[varnames1 %in% desnames1]
## design2
in.vars2 <- varnames2[varnames2 %in% desnames2]

##  Check if there are actually "internal variables"
## design1
if (length(in.vars1) == 0) {
    stop("Specified estimator isn't a function of variables observed in design1")
}
## design2
if (length(in.vars2) == 0) {
    stop("Specified estimator isn't a function of variables observed in design2")
}

## Check "internal variables" type
## design1
typetest1 <- sapply(in.vars1, function(y) is.numeric(design1$variables[, y]))
if (!all(typetest1)) 
    stop("Estimator formula cannot contain non numeric variables")
## design2
typetest2 <- sapply(in.vars2, function(y) is.numeric(design2$variables[, y]))
if (!all(typetest2)) 
    stop("Estimator formula cannot contain non numeric variables")

# Check for missing values in interest variables (one at a time, since only a
# single estimator is involved)
## design1
for (var in in.vars1){
     NA.estvars(design = design1, estvars = var, draconian = TRUE)
    }
## design2
for (var in in.vars2){
     NA.estvars(design = design2, estvars = var, draconian = TRUE)
    }

## Extract from expr the "parameters" i.e. symbols NOT referencing design
## variables
## design1
parameters1 <- varnames1[!(varnames1 %in% desnames1)]
## design2
parameters2 <- varnames2[!(varnames2 %in% desnames2)]

## Parameters must exist and be numeric scalars, check it out
is.scalar <- function(x) {
     x <- get(x)
     is.vector(x, mode = "numeric") && identical(length(x), 1L)
    }
## design1
scalartest1 <- sapply(parameters1, is.scalar)
if (!all(scalartest1)) 
    stop("Estimator formula can contain only scalar numeric parameters")
## design2
scalartest2 <- sapply(parameters2, is.scalar)
if (!all(scalartest2)) 
    stop("Estimator formula can contain only scalar numeric parameters")

## Build a convenience formula with "internal variables"...
## design1
in.vars.formula1 <- as.formula(paste("~",paste(in.vars1, collapse="+"),sep=""))
## design2
in.vars.formula2 <- as.formula(paste("~",paste(in.vars2, collapse="+"),sep=""))

## ... and its model frame
## design1
mf1 <- model.frame(in.vars.formula1, design1$variables, na.action = na.pass)
## design2
mf2 <- model.frame(in.vars.formula2, design2$variables, na.action = na.pass)

## Convert model frame into matrix
## design1
x1 <- as.matrix(mf1)
## design2
x2 <- as.matrix(mf2)

## Compute estimators of totals of "internal variables"
## design1
estimators1 <- as.list(coef(z.svytotal(x1, design = design1)))
## design2
estimators2 <- as.list(coef(z.svytotal(x2, design = design2)))

# Build suitable environments for resolving symbolic derivatives
## design1
env1 <- estimators1
## design2
env2 <- estimators2

## Re-build univocal variable names
## design1
namevec1 <- paste(in.vars1, "1", sep = ".")
names(env1) <- paste(names(env1), "1", sep = ".")
## design2
namevec2 <- paste(in.vars2, "2", sep = ".")
names(env2) <- paste(names(env2), "2", sep = ".")

## Symbolically compute both expr and its first order derivatives w.r.t the "internal variables"
f.plus.grad <- deriv(expr = expr, namevec = c(namevec1, namevec2), function.arg = FALSE)

## Evaluate symbolic expressions
f.plus.grad.val <- eval(f.plus.grad, envir = c(env1, env2))
 # 1. Gradient
grad.val <- attr(f.plus.grad.val, "gradient")

 # 2. Linearized variables
## design1
y.1 <- design1$variables$y.1 <- x1%*%t(grad.val[, namevec1, drop = FALSE])
   ## If calibrated, get the residual w.r.t. the calibration model. This will
   ## eventually feed the RHO calculation
   if (is.calibrated(design1)) {
         y.1 <- design1$variables$ey.1 <- as.numeric(get.residuals(design1, ~y.1, scale = "no"))
        }

## design2
y.2 <- design2$variables$y.2 <- x2%*%t(grad.val[, namevec2, drop = FALSE])
   ## If calibrated, get the residual w.r.t. the calibration model. This will
   ## eventually feed the RHO calculation
   if (is.calibrated(design2)) {
         y.2 <- design2$variables$ey.2 <- as.numeric(get.residuals(design2, ~y.2, scale = "no"))
        }

 # 3. Estimate
Delta <- as.numeric(f.plus.grad.val)

## Final design1 AND design2 Woodruff transforms
attr(Delta, "lin.variables.names") <- c("Ly.1", "Ly.2")

 # 4. Variance estimates of linearized variables from design1 and design2
 #    NOTE: Of course HERE we don't need the calibration residuals even if
 #          design is calibrated, as the residuals will be computed by inner
 #          functions (e.g. svyrecvar)... 
## design1
V1 <- VAR(svystatTM(design1, ~y.1))
## design2
V2 <- VAR(svystatTM(design2, ~y.2))
## design1 AND design2
variances <- c(V1, V2)
names(variances) <- c("V1", "V2")
attr(Delta, "variances") <- variances

 # 3. Designs (with variables to feed RHO)
design1$variables$Ly.1 <- y.1
design2$variables$Ly.2 <- y.2 
attr(Delta, "designs") <- list("design1" = design1, "design2" = design2)

## Return result
Delta
}


`svydelta` <- function(expr, design1, design2, has.strata, is.element, are.indep, no.strat.jump, ...){
################################################################################
# Compute RHO along the lines of Berger and Priam 2016.                        #
################################################################################

## Get estimate, linearized variables, variances, etc.
Delta <- LinDelta(expr = expr, design1 = design1, design2 = design2, ...)
## Estimate
estimate <- as.numeric(Delta)
names(estimate) <- as.character(expr)

## Variances (of linearized variables)
V1 <- attr(Delta, "variances")[["V1"]]
V2 <- attr(Delta, "variances")[["V2"]]
## Designs (with additional variables)
design1 <- attr(Delta, "designs")$design1
design2 <- attr(Delta, "designs")$design2
## Data from designs, taking care of subsetting
## If svydelta has been called on a subset AND design is calibrated get the
## domain index and subset the data, else keep the whole sample (as it is
## already subsetted):
di1 <- attr(design1, "domain.index")
if (!is.null(di1)) {
     data1 <- design1$variables[di1,, drop = FALSE]
    } else {
     data1 <- design1$variables
    }
di2 <- attr(design2, "domain.index")
if (!is.null(di2)) {
     data2 <- design2$variables[di2,, drop = FALSE]
    } else {
     data2 <- design2$variables
    }

# Are design1 and design2 independent?
if (!are.indep) {
    ## DESIGNS ARE NOT INDEPENDENT
    ## HERE MUST DERIVE RHO AND ALL THE REST (CovY1Y2, etc.) --------- START - #

## Get names of first stage identifier and weights
## design1
ids.char1 <- all.vars(attr(design1, "ids"))[1]
w.char1 <- all.vars(attr(design1, "weights"))[1]
## design2
ids.char2 <- all.vars(attr(design2, "ids"))[1]
w.char2 <- all.vars(attr(design2, "weights"))[1]

## In case of stratification, get the *design strata* variable (i.e. measured at
## sampling time by the sampling frame)... paying attention to handle well
## *collapsed strata* designs (if any): avoid the design$strata structure and
## go for the strata column of the embedded data
## !!! REFLECT ON WHETHER IT IS OK TO DO SO, OR IT IS INSTEAD BETTER TO USE
##     COLLAPSED STRATA FOR RHO !!!

if (has.strata) {
     ## design1
     # is.clps1 <- !is.null(attr(design1, "collapse.strata"))
     # get strata variable name
     strata.char1 <- all.vars(attr(design1, "strata"))
     # get strata variable values
     # if (!is.clps1) {
     #     strata1 <- design1$strata[, 1]
     #    } else {
     #     strata1 <- design1$variables[[strata.char1]]
     #   }
     ## design2
     # is.clps2 <- !is.null(attr(design2, "collapse.strata"))
     # get strata variable name
     strata.char2 <- all.vars(attr(design2, "strata"))
     # get strata variable values
     # if (!is.clps2) {
     #     strata2 <- design2$strata[, 1]
     #    } else {
     #     strata2 <- design2$variables[[strata.char2]]
     #   }
    }

## Restrict to the columns actually needed for RHO (* = optional)
## id, y1, w1, ones1, strata1*
## id, y2, w2, ones2, strata2*

## design1
data1$ones1 <- 1
if (!has.strata) {
     data1 <- data1[, c(ids.char1, "Ly.1", w.char1, "ones1")]
     colnames(data1) <- c("id", "y1", "w1", "ones1")
    } else {
     data1 <- data1[, c(ids.char1, "Ly.1", w.char1, "ones1", strata.char1)]
     colnames(data1) <- c("id", "y1", "w1", "ones1", "strata1")
    }
## design2
data2$ones2 <- 1
if (!has.strata) {
     data2 <- data2[, c(ids.char2, "Ly.2", w.char2, "ones2")]
     colnames(data2) <- c("id", "y2", "w2", "ones2")
    } else {
     data2 <- data2[, c(ids.char2, "Ly.2", w.char2, "ones2", strata.char2)]
     colnames(data2) <- c("id", "y2", "w2", "ones2", "strata2")
    }


if (!is.element) {
## CLUSTER sampling designs - START ##
# HERE MUST AGGREGATE data1 and data2 BY PSUs, so to obtain a "element" samples
# whose elements are PSUs and whose Y variables are estimated PSU totals

# Prepare weighted Ys
# data1
data1$wy1 <- data1$y1 * data1$w1
# data2
data2$wy2 <- data2$y2 * data2$w2

# Aggregate at PSU level
# data1
if (!has.strata) {
   # data1 <- aggregate(. ~ id, data = data1, FUN = sum) # DEBUG: wrong col names
     data1 <- aggregate(data1[, c("y1", "w1", "ones1", "wy1")], by = list(id = data1$id), FUN = sum)
    } else {
   # data1 <- aggregate(. ~ strata1 + id, data = data1, FUN = sum) # DEBUG: wrong col names
     data1 <- aggregate(data1[, c("y1", "w1", "ones1", "wy1")], by = list(strata1 = data1$strata1, id = data1$id), FUN = sum)
    }
colnames(data1)[colnames(data1) == "ones1"] <- "nPSU1"
data1$ones1 <- 1
# data2
if (!has.strata) {
   # data2 <- aggregate(. ~ id, data = data2, FUN = sum) # DEBUG: wrong col names
     data2 <- aggregate(data2[, c("y2", "w2", "ones2", "wy2")], by = list(id = data2$id), FUN = sum)
    } else {
   # data2 <- aggregate(. ~ strata2 + id, data = data2, FUN = sum) # DEBUG: wrong col names
     data2 <- aggregate(data2[, c("y2", "w2", "ones2", "wy2")], by = list(strata2 = data2$strata2, id = data2$id), FUN = sum)
    }
colnames(data2)[colnames(data2) == "ones2"] <- "nPSU2"
data2$ones2 <- 1

## Merge data1 and data2
data <- merge(data1, data2, all = TRUE)

# Deal with NAs in data
# Model Zs (these use real strata if designs are stratified, else numeric
# dummies that actually account for rotation/overlap)
# data1
data$Z1 <- data$ones1
data$Z1[is.na(data$Z1)] <- 0
# NOTE: To save memory for stratified designs will use a more parsimonious
#       multivariate regression model than proposed by Berger & Priam's Z1 * Z2
data$z1 <- data$Z1
if (has.strata) {
     data$Z1 <- data$strata1
     data$Z1 <- as.character(data$Z1)
     data$Z1[is.na(data$Z1)] <- "NotIn1"
     data$Z1 <- factor(data$Z1)
    }
# data2
data$Z2 <- data$ones2
data$Z2[is.na(data$Z2)] <- 0
# NOTE: To save memory for stratified designs will use a more parsimonious
#       multivariate regression model than proposed by Berger & Priam's Z1 * Z2
data$z2 <- data$Z2
if (has.strata) {
     data$Z2 <- data$strata2
     data$Z2 <- as.character(data$Z2)
     data$Z2[is.na(data$Z2)] <- "NotIn2"
     data$Z2 <- factor(data$Z2)
    }

# Weighted Ys
# data1
data$wy1[is.na(data$wy1)] <- 0
# data2
data$wy2[is.na(data$wy2)] <- 0

## CLUSTER sampling designs - END   ##

} else {

## ELEMENT sampling designs - START ##

## Merge data1 and data2
data <- merge(data1, data2, all = TRUE)

# Deal with NAs in data
# Model Zs (these use real strata if designs are stratified, else numeric
# dummies that actually account for rotation/overlap)
# data1
data$Z1 <- data$ones1
data$Z1[is.na(data$Z1)] <- 0
# NOTE: To save memory for stratified designs, if no.strat.jump = TRUE, will use
#       a more parsimonious multivariate regression model than proposed by
#       Berger & Priam's Z1 * Z2
data$z1 <- data$Z1
if (has.strata) {
     data$Z1 <- data$strata1
     data$Z1 <- as.character(data$Z1)
     data$Z1[is.na(data$Z1)] <- "NotIn1"
     data$Z1 <- factor(data$Z1)
    }
# data2
data$Z2 <- data$ones2
data$Z2[is.na(data$Z2)] <- 0
# NOTE: To save memory for stratified designs, if no.strat.jump = TRUE, will use
#       a more parsimonious multivariate regression model than proposed by
#       Berger & Priam's Z1 * Z2
data$z2 <- data$Z2
if (has.strata) {
     data$Z2 <- data$strata2
     data$Z2 <- as.character(data$Z2)
     data$Z2[is.na(data$Z2)] <- "NotIn2"
     data$Z2 <- factor(data$Z2)
    }

# Ys
# data1
data$y1[is.na(data$y1)] <- 0
# data2
data$y2[is.na(data$y2)] <- 0

# Weights
# data1
data$w1[is.na(data$w1)] <- 0
# data2
data$w2[is.na(data$w2)] <- 0

# Prepare weighted Ys
# data1
data$wy1 <- data$y1 * data$w1 
# data2
data$wy2 <- data$y2 * data$w2 

## ELEMENT sampling designs - END   ##
}

# Fit Berger and Priam multivariate regression model to estimate RHO
# NOTE: If has.strata == FALSE the *OVERLAP is COMPLETE*, then Z1 = Z2 and
#       therefore the inclusion in the model of their interaction term would
#       generate an error. Conclusion: the model should be handled differently
#       in case of complete overlap (which, by the way, Berger and Priam define
#       not neatly as g = 1 - using their g = nc / n1 -, which could happen even
#       if the samples are not exactly the same, when s1 is a subset of s2...)

# Empirical overlap rate (at stage 1, in general)
# First stage units
# design1
n1 <- NROW(data1)
# design2
n2 <- NROW(data2)
# common to design1 and design2
nc <- sum(data$ones1 * data$ones2, na.rm = TRUE)
# merge of design1 and design2 (nm = n1 + n2 - nc)
nm <- NROW(data)
# Overlap rate
overlap.rate <- nc / mean(c(n1, n2))

# Is the overlap full?
full.overlap <- ( n1 == n2 ) && ( nm == n1 )

# Unstratified case
if (!has.strata) {
     if (full.overlap) {
          mm <- lm(formula = cbind(wy1, wy2) ~ z1 - 1, data = data)
        } else {
          mm <- lm(formula = cbind(wy1, wy2) ~ (z1 * z2) - 1, data = data)
        }
    }
# Stratified case
if (has.strata) {
     # NOTE: If has.strata it is still meaningful to keep strata interactions
     #       even if overlap is full, to take care of dynamic strata (e.g.
     #       stratum jumpers from 1 to 2 in the *domain*)
     #
     # NOTE: However, one has to take care of the special case when the user
     #       asks for domain estimates at *stratum* level, as this means that
     #       some of the two Z factors can happen to have just 1 level, which
     #       would result into an error of lm()
     if ( (length(levels(data$Z1)) == 1) ||  (length(levels(data$Z2)) == 1) ) {
         # NOTE: This condition means the user asked for domain estimates at
         #       *stratum* level
          mm <- lm(formula = cbind(wy1, wy2) ~ (z1 * z2) - 1, data = data)
        } else {
          if (no.strat.jump) {
              # Save memory and processing time: obtain correct results if there
              # are no stratum-jumpers from design1 to design2
              mm <- lm(formula = cbind(wy1, wy2) ~ (z1 * z2) : (Z1 + Z2) - 1, data = data)
            }
          else {
              # Full complexity: fully accommodates dynamic stratification, but
              # likely unaffordable for highly stratified sampling design
              # (especially of units)
              mm <- lm(formula = cbind(wy1, wy2) ~ (Z1 * Z2) - 1, data = data)
            }
        }
    }

# Get the residual variance covariance matrix from the model fit
VV <- estVar(mm)

# Estimate the correlation between Y1 and Y2
# Model covariance
VV12 <- mean(VV[row(VV) != col(VV)])
# Product of the square roots of the model variances 
VV11sqVV22sq <- sqrt(prod(diag(VV)))
# Correlation (standard difference form LDelta)
RHO <- - VV12 / VV11sqVV22sq

# Estimate of the covariance between Y1 and Y2 (standard difference form LDelta)
covY1Y2 <- RHO * sqrt(V1) * sqrt(V2)

# Estimate of the variance of the measure of change Delta
VDelta <- V1 + V2 - 2 * covY1Y2

## Build a diagnostic data frame
details <- data.frame("n1"       = n1,
                      "n2"       = n2,
                      "nc"       = nc,
                      "overlap"  = overlap.rate,
                      "V1"       = V1,
                      "V2"       = V2,
                      "Vind"     = V1 + V2,
                      "rho"      = RHO,
                      "CoV"      = covY1Y2,
                      "V"        = VDelta)

    ## HERE MUST DERIVE RHO AND ALL THE REST (CovY1Y2, etc.) ----------- END - #
    } else {
    ## DESIGNS ARE INDEPENDENT
     VDelta <- V1 + V2
     details <- data.frame("msg" = "Input design objects were declared to be independent")
    }

## Attach key attributes to the estimate
# Estimated variance of Delta
attr(estimate, "var") <- as.matrix(VDelta)
dimnames(attr(estimate, "var"))<-list(names(estimate),names(estimate))

# Details on the calculation of the variance (especially the covariance term)
attr(estimate,"details") <- details

# Kind of statistic returned
attr(estimate,"statistic") <- "Delta"

# Class of return object
class(estimate) <- c("svydelta", "svystat")

return(estimate)
}

`coef.svydelta` <- function (object, ...) 
{
    out <- as.numeric(unclass(object))
    names(out) <- names(object)
    out
}
