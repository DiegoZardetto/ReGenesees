\name{extractors}
\alias{extractors}
\alias{SE}
\alias{SE.default}
\alias{VAR}
\alias{VAR.default}
\alias{cv}
\alias{cv.default}
\alias{deff}
\alias{deff.default}
\title{Extractor Functions for Variability Statistics}
\description{
These functions extract standard errors (\code{SE}), variances (\code{VAR}), coefficients of variation (\code{cv}) and design effects (\code{deff}) from an object which has been returned by a survey statistic function (e.g. \code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}}, \code{\link{svystat}}, \ldots).
}
\usage{
SE(object, ...)
VAR(object, ...)
cv(object, ...)
deff(object, ...)
}
\arguments{
  \item{object}{An object containing survey statistics.}
  \item{\dots}{Arguments for future expansion.}
}
\details{
With the exception of \code{deff}, all extractor functions can be used on any object returned by a survey statistic function: the correct answer will be obtained whatever the call that generated the object. For getting the design effect, \code{object} must have been built with option \code{deff = TRUE}.  
}
\value{
A data structure (typically inheriting from classes \code{matrix} or \code{data.frame}) storing the requested information.
}
\note{
Package \pkg{ReGenesees} provides extensions of methods \code{\link{coef}} and \code{\link{confint}} (originally from package \pkg{stats}) that can be used to extract estimates and confidence intervals respectively.
}

\author{
Diego Zardetto
}

\seealso{
Function \code{\link{coef}} to extract estimates and function \code{\link{confint}} to extract confidence intervals. Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
# Creation of a design object:
data(sbs)
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
     fpc=~fpc)

# Estimation of the average value added at the
# nation level (by default one gets the SE):
VA.avg <- svystatTM(des,~va.imp2,estimator="Mean")
VA.avg

# Extractions of some variance statistics from the
# object above:
## 1) SE
SE(VA.avg)
## 2) CV
cv(VA.avg)
## 3) VAR
VAR(VA.avg)

# Design effects have to be requested in advance,
# i.e. the following invocation produces an error:
\dontrun{
deff(VA.avg)
}
# ...while the following works:
VA.avg <- svystatTM(des,~va.imp2,estimator="Mean",deff=TRUE)
deff(VA.avg)

# Further examples:
## extract the statistic:
coef(VA.avg)
## extract the confidence interval at 90%
## confidence level (the default would be 95%):
confint(VA.avg, conf.lev=0.9)
}

\keyword{survey}