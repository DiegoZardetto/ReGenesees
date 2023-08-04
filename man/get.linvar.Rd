\name{get.linvar}
\alias{get.linvar}
\title{Linearized Variable(s) of Complex Estimators by Domains}
\description{
Computes the linearized variable(s) of a Complex Estimator in subpopulations (domains). The Complex Estimator can be any analytic function of Horvitz-Thompson or Calibration estimators.
}
\usage{
get.linvar(design, expr, by = NULL, stack = FALSE, na.rm = FALSE)
}
\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata. Can be either uncalibrated or calibrated.}
  \item{expr}{R expression defining the Complex Estimator (see \code{\link{svystatL}} for the basic syntax; see also section \sQuote{Details}).}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) the estimator refers to the whole population.}
  \item{stack}{If \code{FALSE} (the default), the function will return a \code{matrix} of domain-specific linearized variables, with one column per domain. If \code{TRUE}, the domain-specific linearized variables will be stacked without information loss and returned as a single column. The \emph{stacked} output format is only possible for \emph{uncalibrated} designs (see \sQuote{Details}).}
  \item{na.rm}{Allow missing values in the variables entering the estimator expression? The default is
              \code{FALSE}.}
}
\details{
This function has been designed mainly for programmers or advanced users willing to build upon \pkg{ReGenesees}: typical users are not expected to feel much need of it.

The \pkg{ReGenesees} package adopts the Taylor-series linearization technique to estimate the sampling variance of nonlinear (smooth) estimators. Using this technique, it can handle any Complex Estimator that can be expressed as an analytic function of Horvitz-Thompson or Calibration estimators. This includes even user-defined estimators, through function \code{\link{svystatL}} [Zardetto, 15].

In the Taylor-series linearization approach, estimating the sampling variance of a Complex Estimator amounts to estimating the sampling variance of the Horvitz-Thompson estimator of the total of the \emph{linearized variable} of the Complex Estimator under the sampling design at hand. Function \code{get.linvar} computes that linearized variable. If one is interested in estimating the sampling variance of a Complex Estimator by \emph{domains} (i.e. for subpopulations), then - owing to the Taylor linearization technique - a different linearized variable has to be calculated for each domain. Therefore, in the most general case, the output of function \code{get.linvar} will be a \code{matrix} of linearized variables, with one column per domain.

If the input object \code{design} is calibrated, then the linearized variable produced by function \code{get.linvar} correctly takes into account the theoretical properties of Calibration estimators (see, e.g., [Deville, 1999]). In particular, \code{get.linvar} computes the needed g-weighted residuals (see, e.g., equations (8) and (9) of [Zardetto, 15]) by leveraging function \code{\link{get.residuals}}.

The mandatory argument \code{expr}, which identifies the Complex Estimator, must be an object of class \code{expression}. You can specify just a single Complex Estimator at a time, i.e. \code{length(expr)} must be equal to \code{1}. Any analytic function of estimators of Totals and Means derived from the input design object \code{design} is allowed.

The basic syntax follows the same rules described for function \code{\link{svystatL}}. Inside \code{expr}, the estimator of the Total of a variable is simply represented by the \emph{name} of the variable itself. The reserved name \code{ones} can be used to reference an artificial variable (which will be created on-the-fly, if not already present) whose value is \code{1} for each sampling unit, so that its Total estimator actually estimates the size of the population in terms of elementary units. Therefore, e.g., expression \code{y/ones} represents the estimator of the Mean of variable \code{y}. Variables referenced inside \code{expr} must be \code{numeric} and belong to \code{design}.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatL} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatL} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced.

As already noted, in the most general case, the output of function \code{get.linvar} will be a \code{matrix} of linearized variables, with one column per domain. However, the linearized variable of any \emph{domain estimator} that is a function of Horvitz-Thompson estimators is identically equal to zero outside the domain. Thanks to this property, and only for functions of HT estimators, the columns of the output matrix can be stacked into one single column without information loss. Users may request this convenient output format by specifying \code{stack = TRUE}. Unfortunately, the same result cannot be obtained for functions of Calibration estimators, because the g-weighted residuals are generally non-zero outside the active estimation domain. Therefore, argument \code{stack} does nothing if \code{design} is calibrated.
}

\value{
A matrix of linearized variables.
}

\author{Diego Zardetto}

\references{
Deville, J. C. (1999) \emph{\dQuote{Variance Estimation for Complex Statistics and Estimators: Linearization and Residual Techniques.}}. Survey Methodology 25: 193-203.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\seealso{
\code{\link{svystatL}} to calculate estimates and sampling errors of Complex Estimators, \code{\link{get.residuals}} to compute residuals of interest variables w.r.t. the calibration model adopted to build a calibrated object.
}

\examples{
# Load sbs data:
data(sbs)

# Create a design object to be calibrated:
sbsdes<-e.svydesign(data= sbs, ids= ~id, strata= ~strata, weights= ~weight, fpc= ~fpc)

##################################################################
# Just some checks on the consistency of the numerical results   #
# obtained by ReGenesees with well known theoretical properties. #
##################################################################

## In the Taylor-series linearization approach, estimating the sampling variance of a
## Complex Estimator amounts to estimating the sampling variance of the Horvitz-Thompson
## estimator of the total of the linearized variable of the Complex Estimator under the
## sampling design at hand.

#############
## Check 1 ##
#############
# Global estimate (i.e., no domains) of a nonlinear function of HT estimators:
# - Average value added per employee

# Compute the linearized variable of the estimator:
z <- get.linvar(sbsdes, expression(va.imp2/emp.num))

# Have a look: 
head(z)

# Now add that variable to the input design object:
sbsdes <- des.addvars(sbsdes, z = z)

# Now compute the estimated standard error of the HT total of the linearized variable:
svystatTM(sbsdes, ~z)

# And compare it to the estimated standard error of the complex estimator:
svystatL(sbsdes, expression(va.imp2/emp.num))

# ...which - in this case - can equivalently be obtained exploiting its ratio nature:
svystatR(sbsdes, num = ~va.imp2, den= ~emp.num)

# Test the equality of the results:
all.equal(svystatTM(sbsdes, ~z)$SE, 
          svystatL(sbsdes, expression(va.imp2/emp.num))$SE)
# OK

#############
## Check 2 ##
#############
# Domain estimates of a nonlinear function of HT estimators:
# - Average value added per employee by NACE macro-sectors

# Compute the linearized variable(s) of the domain estimator(s):
z.dom <- get.linvar(sbsdes, expression(va.imp2/emp.num), by = ~nace.macro)

# Have a look: 
head(z.dom)
tail(z.dom)

# Note the staggered matrix-format of the result. Since you are dealing with HT estimators
# (namely sbsdes is not calibrated), you can opt for the more convenient (but equivalent)
# stacked format:
z.dom <- get.linvar(sbsdes, expression(va.imp2/emp.num), by = ~nace.macro, stack = TRUE)

# Have a look: 
head(z.dom)
tail(z.dom)

# Now add that variable to the input design object:
sbsdes <- des.addvars(sbsdes, z.dom = z.dom)

# Now compute the estimated standard error of the HT total of the linearized variable by
# domains:
svystatTM(sbsdes, ~z.dom, by= ~nace.macro)

# And compare it to the estimated standard error of the complex estimator:
svystatL(sbsdes, expression(va.imp2/emp.num), by= ~nace.macro)

# ...which - in this case - can equivalently be obtained exploiting its ratio nature:
svystatR(sbsdes, num = ~va.imp2, den= ~emp.num, by= ~nace.macro)

# Test the equality of the results:
all.equal(svystatTM(sbsdes, ~z.dom, by= ~nace.macro)$SE,
          svystatL(sbsdes, expression(va.imp2/emp.num), by= ~nace.macro)$SE)
# OK

###################################################################
# Now run the same checks for functions of Calibration Estimators #
###################################################################

# Build a population totals template and fill it with actual known totals computed from
# the sampling frame (sbs.frame):
pop <- pop.template(sbsdes, calmodel= ~ent:emp.cl + y:nace.macro-1, partition= ~region)
pop <- fill.template(sbs.frame, pop)

# Calibrate the weights:
sbscal <- e.calibrate(sbsdes, pop)
g.range(sbscal)

#############
## Check 3 ##
#############
# Global estimate (i.e., no domains) of a nonlinear function of CAL estimators:
# - Average value added per employee

# Compute the linearized variable of the estimator (the name gez should remind you of
# the g-weighted residuals of the HT linearized variable w.r.t. the calibration model):
gez <- get.linvar(sbscal, expression(va.imp2/emp.num))

# Have a look: 
head(gez)

# Now add that variable to the initial *uncalibrated* design object:
sbsdes <- des.addvars(sbsdes, gez = gez)

# Now compute the estimated standard error of the *HT total* of the linearized variable:
svystatTM(sbsdes, ~gez)

# And compare it to the estimated standard error of the complex *calibration* estimator:
svystatL(sbscal, expression(va.imp2/emp.num))

# ...which - in this case - can equivalently be obtained exploiting its ratio nature:
svystatR(sbscal, num = ~va.imp2, den= ~emp.num)

# Test the equality of the results:
all.equal(svystatTM(sbsdes, ~gez)$SE,
          svystatL(sbscal, expression(va.imp2/emp.num))$SE)
# OK

#############
## Check 4 ##
#############
# Domain estimates of a nonlinear function of CAL estimators:
# - Average value added per employee by NACE macro-sectors

# Compute the linearized variable(s) of the domain estimator(s):
gez.dom <- get.linvar(sbscal, expression(va.imp2/emp.num), by = ~nace.macro)

# Have a look: 
head(gez.dom)
tail(gez.dom)

# As you see the matrix no longer has the staggered structure of the HT case. In fact, the
# calibration residuals are in general not 0 outside the active domain. Therefore, you
# cannot opt for the more convenient stacked format!

# Let's thus add all the column of the linearized variable(s) matrix to the initial
# *uncalibrated* design object:
sbsdes$variables <- cbind(sbsdes$variables, gez.dom)

# Now compute the estimated standard error of the *HT totals* of the *domain-specific*
# linearized variables:
svystatTM(sbsdes, ~gez_Agriculture + gez_Industry + gez_Commerce + gez_Services)

# And compare it to the estimated standard errors of the complex *calibration* estimator
# by domains:
svystatL(sbscal, expression(va.imp2/emp.num), by= ~nace.macro)

# ...which - in this case - can equivalently be obtained exploiting its ratio nature:
svystatR(sbscal, num = ~va.imp2, den= ~emp.num, by= ~nace.macro)

# Test the equality of the results:
all.equal(
  svystatTM(sbsdes, ~gez_Agriculture + gez_Industry + gez_Commerce + gez_Services)$SE,
  svystatL(sbscal, expression(va.imp2/emp.num), by= ~nace.macro)$SE)
# OK

}
\keyword{survey}