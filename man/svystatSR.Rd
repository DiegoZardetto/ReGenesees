\name{svystatSR}
\alias{svystatSR}
\alias{coef.svystatSR}
\alias{SE.svystatSR}
\alias{VAR.svystatSR}
\alias{cv.svystatSR}
\alias{deff.svystatSR}
\alias{confint.svystatSR}

\title{Estimation of Share Ratios in Subpopulations}

\description{
Calculates estimates, standard errors and confidence intervals for Ratios between Shares of a numeric variables in subpopulations.
}

\usage{
svystatSR(design, y, classes, by = NULL,
          vartype = c("se", "cv", "cvpct", "var"),
          conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
          na.rm = FALSE)

\method{coef}{svystatSR}(object, ...)
\method{SE}{svystatSR}(object, ...)
\method{VAR}{svystatSR}(object, ...)
\method{cv}{svystatSR}(object, ...)
\method{deff}{svystatSR}(object, ...)
\method{confint}{svystatSR}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{y}{Formula defining the interest variable.}
  \item{classes}{Formula defining the population groups among which ratios of \code{y} shares must be estimated.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.int}{Compute confidence intervals for the estimates? The default is
                 \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{deff}{Should the design effect be computed? The default is \code{FALSE} (see \sQuote{Details}).}
  \item{na.rm}{Should missing values (if any) be removed from the variables of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
  \item{object}{An object of class \code{svystatSR}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Ratios between Shares of a numeric variable, using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors are calculated using the Taylor linearization technique.

Ratios of Shares are a special case of Ratios. Therefore, at the price of some additional and heavy data preparation effort, ratios of shares could also be estimated using function \code{\link{svystatR}}. However, \code{svystatSR} makes estimation by far easier, in particular when share ratios have to be estimated for many population groups and/or within many domains.

The mandatory argument \code{classes} identifies population groups whose ratios of \code{y} shares have to be estimated. Note that ratios of shares will be estimated and returned \emph{for all the ordered pairs of population groups} defined by \code{classes}. Therefore, if \code{classes} defines \code{G} groups, \code{svystatSR} will have to compute estimates and sampling errors for \code{G * (G - 1)} share ratios. To prevent combinatorial explosions (e.g. \code{G = 20} would generate \code{380} share ratios), \code{classes} formula can reference just a \emph{single} \code{design} variable, which must be a \code{factor}.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatSR} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatSR} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates.  By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the variance of the estimator under the actual sampling design and the variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}.\cr
Being Ratios nonlinear estimators, the design effect is estimated on the linearized version of the estimator (that is: for the estimator of the total of the linearized variable, aka "Woodruff transform").\cr
When dealing with domain estimation, the design effects referring to a given subpopulation are currently computed by taking the ratios between the actual variance estimates and those that would have been obtained if a simple random sampling were carried out \emph{within} that subpopulation. This is the same as the \code{srssubpop} option for Stata's function \code{estat}.

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

European Commission, Eurostat, (2013). \emph{\dQuote{Handbook on precision requirements and variance estimation for ESS households surveys: 2013 edition}}, Publications Office. doi: 10.2785/13579
}

\seealso{
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
# Load household data:
data(data.examples)

# Create a design object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
        weights=~weight)

# Add convenience variable 'ones' to estimate counts:
des<-des.addvars(des,ones=1)

### Simple examples to illustrate the syntax:
# Population sex ratios:
svystatSR(des, y=~ones, classes=~sex, vartype="cvpct")

# Population sex ratios within provinces:
svystatSR(des, y=~ones, classes=~sex, by=~procod, vartype="cvpct")

# Ratios of population shares for 5 age classes:
# NOTE: This yields 5*(5-1)=20 ratios
svystatSR(des, y=~ones, classes=~age5c, vartype="cvpct")


### One more complicated example:
#######################################################################
# Ratios between shares of income held by people for income quintiles #
#######################################################################
    # First: estimate income quintiles
inc.Q5 <- svystatQ(des, y=~income, probs=seq(0.2, 0.8, 0.2), ties="rounded")
inc.Q5

    # Second: add a convenience factor variable classifying people by income
    # quintiles 
des<-des.addvars(des, quintile = cut(income, breaks = c(0, coef(inc.Q5), Inf),
                                     labels = 1:5, include.lowest=TRUE)
                    )

    # Third: estimate income shares by income quintiles
QS5 <- svystatSR(des, y=~income, classes=~quintile, vartype="cvpct")
QS5

### Therefore, for instance, the *S80/S20 income quintile share ratio* is:
S80.20 <- QS5["quintile5/quintile1",]
S80.20

### NOTE: Procedure above yields *correct point estimates* of income quintile
###       share ratios, while *variance estimation is approximated* since
###       we neglected the sampling variability of the estimated quintiles.
}
\keyword{survey}