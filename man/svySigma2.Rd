\name{svySigma2}
\alias{svySigma2}
\alias{coef.svySigma2}
\alias{SE.svySigma2}
\alias{VAR.svySigma2}
\alias{cv.svySigma2}
\alias{confint.svySigma2}
\title{Estimation of the Population Variance of a Variable}
\description{
Computes estimates and sampling errors of the Population Variance of a numeric variable (in subpopulations too).
}
\usage{
svySigma2(design, y, by = NULL,
          fin.pop = TRUE,
          vartype = c("se", "cv", "cvpct", "var"),
          conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
          na.rm = FALSE)

\method{coef}{svySigma2}(object, ...)
\method{SE}{svySigma2}(object, ...)
\method{VAR}{svySigma2}(object, ...)
\method{cv}{svySigma2}(object, ...)
\method{confint}{svySigma2}(object, ...)
}
\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{y}{Formula identifying the \emph{numeric} interest variable.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{fin.pop}{If \code{TRUE} (the default) the estimation target is the \emph{finite} population formula of the variance, i.e. the one with the \bold{N - 1} denominator. If \code{FALSE} the estimation target is the population variance with the \bold{N} denominator.} 
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.int}{Compute confidence intervals for the estimates? The default is
                 \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{deff}{Should the design effect be computed? The default is \code{FALSE} (see \sQuote{Details}).}
  \item{na.rm}{Should missing values (if any) be removed from the variable of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
  \item{object}{An object of class \code{svySigma2}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}
\details{
Function \code{svySigma2} computes estimates and sampling errors of the Population Variance of a numeric variable. These estimates play an important role in many contexts, including sample size guesstimation and power calculations.

As the Population Variance is a \emph{complex} estimator, \code{svySigma2} automatically linearizes it to estimate its sampling variance. Automatic linearization is performed as function \code{\link{svystatL}} would do, along the lines illustrated in [Zardetto, 15]. This, of course, also entails the usage of the residuals technique when the input \code{design} object is calibrated (i.e. of class \code{cal.analytic}).

The mandatory argument \code{y} identifies the variable of interest. The \code{design} variable referenced by \code{y} must be \code{numeric}.

If variable \code{y} is binary (i.e. has only values \code{0} and \code{1}), the estimated Population Variance coincides with the classical Bernoulli expression \code{p*(1 - p)}, where \code{p} is the estimated proportion of population units with \code{y = 1} (see \sQuote{Examples}).

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svySigma2} refer to the whole population. If specified, estimation domains must be defined by a formula, following the usual syntactic and semantic rules (see e.g. \code{\link{svystatTM}}).

Argument \code{fin.pop} allows the users to select which population variance formula they prefer as estimation target. If \code{fin.pop = TRUE} (the default) the \emph{finite} population version of the variance formula will be used, namely the one with the \bold{N - 1} denominator [Sarndal, Swensson, Wretman 92]. If \code{fin.pop = FALSE} the population variance formula with the \bold{N} denominator will be used.

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates. By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the sampling variance of the estimator under the actual sampling design and the sampling variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}. See \code{\link{svystatTM}} for further details.

Missing values (\code{NA}) in interest variables should be avoided. If \code{na.rm=FALSE} (the default) they generate NAs in estimates (or even an error, if \code{design} is calibrated). If \code{na.rm=TRUE}, observations containing NAs are dropped, and estimates get computed on non missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random}: should this not be the case, computed estimates would be \emph{biased}.
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values.
}

\author{
Diego Zardetto
}

\references{
Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.

Kish, L. (1995). \emph{\dQuote{Methods for design effects}}. Journal of Official Statistics, Vol. 11, pp. 55-77.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\seealso{
Function \code{\link{svySigma}} to estimate the Population Standard Deviation of a numeric variable. Estimators of Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}. Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, and all of the above \code{\link{svystat}}.
}

\examples{
## Load sbs data and create a design object:
data(sbs)
sbsdes <- e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
          fpc=~fpc)

# Estimation of the population variance of value added (variable 'va.imp2'):
svySigma2(sbsdes, ~va.imp2, vartype = "cvpct", conf.int = TRUE, deff = TRUE)

# Compare with the true value computed from the sampling frame ('sbs.frame'):
var(sbs.frame$va.imp2)

# The same as above, by classes of macro-class of economic activity ('nace.macro'):
svySigma2(sbsdes, ~va.imp2, ~nace.macro, vartype = "cvpct", conf.int = TRUE)

# Compare with the true value computed from the sampling frame ('sbs.frame'):
tapply(sbs.frame$va.imp2, sbs.frame$nace.macro, var)

## An example with a binary variable
# Load household data and create a design object:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Build the indicator variable of the 'widowed' marital status:
des<-des.addvars(des, is.widowed = as.numeric(marstat == "widowed"))

# Estimate and store the population proportion of widowed people:
svystatTM(des, ~is.widowed, estimator = "Mean")
# which of course is equal to what one would get directly:
svystatTM(des, ~marstat, estimator = "Mean")

# Store only the estimated proportion
p.widowed <- coef(svystatTM(des, ~is.widowed, estimator = "Mean"))

# Now estimate the population variance of the binary variable 'is.widowed' *with
# fin.pop = FALSE*, and verify that it *exactly* equals the Bernoulli expression
# p.widowed * (1 - p.widowed)
svySigma2(des, ~is.widowed, fin.pop = FALSE, conf.int = TRUE)
p.widowed * (1 - p.widowed)

# ...as it must be.
}

\keyword{survey}