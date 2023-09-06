\name{svystatS}
\alias{svystatS}
\alias{coef.svystatS}
\alias{SE.svystatS}
\alias{VAR.svystatS}
\alias{cv.svystatS}
\alias{deff.svystatS}
\alias{confint.svystatS}

\title{Estimation of Shares in Subpopulations}

\description{
Calculates estimates, standard errors and confidence intervals for Shares of a numeric variable within subpopulations.
}

\usage{
svystatS(design, y, classes, by = NULL,
         vartype = c("se", "cv", "cvpct", "var"),
         conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
         na.rm = FALSE)

\method{coef}{svystatS}(object, ...)
\method{SE}{svystatS}(object, ...)
\method{VAR}{svystatS}(object, ...)
\method{cv}{svystatS}(object, ...)
\method{deff}{svystatS}(object, ...)
\method{confint}{svystatS}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{y}{Formula defining the interest variable.}
  \item{classes}{Formula defining the population groups whose \code{y} shares must be estimated.}
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
  \item{object}{An object of class \code{svystatS}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Shares of a numeric variable, using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors are calculated using the Taylor linearization technique.

Shares are a special case of Ratios. Therefore, at the price of some additional (and possibly heavy) data preparation effort, shares could also be estimated using function \code{\link{svystatR}}. However, \code{svystatS} makes estimation by far easier, in particular when shares have to be estimated for many population groups and/or within many domains.

The mandatory argument \code{y} identifies the variable of interest, that is the variable for which estimates of shares have to be calculated. The \code{design} variable referenced by \code{y} must be \code{numeric}.

The mandatory argument \code{classes} identifies population groups whose shares of \code{y} have to be estimated. The \code{design} variables referenced by \code{classes} must be of class \code{factor}. Groups can be identified by crossing factors, e.g. statement \code{classes = ~C1:C2} selects as groups the subpopulations determined by crossing the levels of factors \code{C1} and \code{C2}.

The optional argument \code{by} specifies the variables defining the "estimation domains", that is the subpopulations within which shares of \code{y} by \code{classes} must be estimated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatS} refer to the whole population. Estimation domains must be defined by a formula: for instance the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatS} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

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
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
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
# Shares of income for sex classes:
svystatS(des, y=~income, classes=~sex, vartype="cvpct")

# Shares of income for sex and 5 age classes:
svystatS(des, y=~income, classes=~age5c:sex, vartype="cvpct")

# Shares of income for sex classes within region domains:
svystatS(des, y=~income, classes=~sex, by=~regcod, vartype="cvpct")

# Shares of income for sex classes within domains defined by crossing region and
# 5 age classes:
svystatS(des, y=~income, classes=~sex, by=~age5c:regcod, vartype="cvpct")

# MARGINAL, CONDITIONAL and JOINT relative frequencies (see also ?svystatTM)
  # MARGINAL: e.g. proportions of people by provinces:
  svystatS(des, y=~ones, classes=~procod, vartype="cvpct")
  # CONDITIONAL: e.g. proportions of people by sex within provinces:
  svystatS(des, y=~ones, classes=~sex, by=~procod, vartype="cvpct")
  # JOINT: e.g. proportions of people cross-classified by sex and procod:
  svystatS(des, y=~ones, classes=~sex:procod, vartype="cvpct")


### One more complicated example:
########################################################
# Shares of income held by people for income quintiles #
########################################################
    # First: estimate income quintiles
inc.Q5 <- svystatQ(des, y=~income, probs=seq(0.2, 0.8, 0.2), ties="rounded")
inc.Q5

    # Second: add a convenience factor variable classifying people by income
    # quintiles 
des<-des.addvars(des, quintile = cut(income, breaks = c(0, coef(inc.Q5), Inf),
                                     labels = 1:5, include.lowest=TRUE)
                    )

    # Third: estimate income shares by income quintiles
svystatS(des, y=~income, classes=~quintile, vartype="cvpct")

### NOTE: Procedure above yields *correct point estimates* of income shares by
###       income quintiles, while *variance estimation is approximated* since
###       we neglected the sampling variability of the estimated quintiles.
}
\keyword{survey}