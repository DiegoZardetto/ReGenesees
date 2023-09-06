\name{svystatQ}
\alias{svystatQ}
\alias{coef.svystatQ}
\alias{SE.svystatQ}
\alias{VAR.svystatQ}
\alias{cv.svystatQ}
\alias{confint.svystatQ}
\title{Estimation of Quantiles in Subpopulations}
\description{
Calculates estimates, standard errors and confidence intervals for Quantiles of numeric variables in subpopulations.
}
\usage{
svystatQ(design, y, probs = c(0.25, 0.5, 0.75), by = NULL,
         vartype = c("se", "cv", "cvpct", "var"),
         conf.lev = 0.95, na.rm = FALSE,
         ties=c("discrete", "rounded"))

\method{coef}{svystatQ}(object, ...)
\method{SE}{svystatQ}(object, ...)
\method{VAR}{svystatQ}(object, ...)
\method{cv}{svystatQ}(object, ...)
\method{confint}{svystatQ}(object, ...)
}
\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{y}{Formula defining the interest variable.}
  \item{probs}{Vector of probability values to be used to calculate the quantiles estimates. The default value selects estimates of quartiles.} 
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{na.rm}{Should missing values (if any) be removed from the variable of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
  \item{ties}{How should duplicated observed values be treated? Select \code{'discrete'} for a genuinely discrete
              interest variable and \code{'rounded'} for a continuous one.}
  \item{object}{An object of class \code{svystatQ}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}
\details{
This function calculates weighted estimates for the Quantiles of a quantitative variable using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise.

Standard errors are calculated using the so-called "Woodruff method" [Woodruff 52][Sarndal, Swensson, Wretman 92]: (i) first a confidence interval (at a given confidence level 1-a) is constructed for the relative frequency of units with values below the estimated quantile, (ii) then the inverse of the estimated cumulative relative frequency distribution (ECDF) is used to map this interval to a confidence interval for the quantile, (iii) lastly the desired standard error is estimated by dividing the length of the obtained confidence interval by the value 2*qnorm(1-a/2). Notice that the procedure above builds, in general, asymmetric confidence intervals around the estimated quantiles.

The mandatory argument \code{y} identifies the variable of interest, that is the variable for which estimates of quantiles have to be calculated. The \code{design} variable referenced by \code{y} must be \code{numeric}.

The optional argument \code{probs} specifies the probability values (\code{0.001<=probs[i]<=0.999}) corresponding to the quantiles one wants to estimate; the default option selects quartiles.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatQ} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatQ} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

Missing values (\code{NA}) in interest variables should be avoided. If \code{na.rm=FALSE} (the default) they generate NAs in estimates (or even an error, if \code{design} is calibrated). If \code{na.rm=TRUE}, observations containing NAs are dropped, and estimates get computed on non missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random}: should this not be the case, computed estimates would be \emph{biased}.

Argument \code{ties} addresses the problem of how to treat duplicated observed values (if any) when computing the ECDF. Option \code{'discrete'} (the default) is appropriate when the variable of interest is genuinely discrete, while \code{'rounded'} is a better choice for a continuous variable, i.e. when duplicates stem from rounding. In the first case the ECDF will show a vertical step corresponding to a duplicated value, in the second a smoother shape will be achieved by linear interpolation.
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values.
}

\author{
Diego Zardetto
}

\references{
Woodruff, R.S. (1952) \emph{\dQuote{Confidence Intervals for Medians and Other Position Measures}}, Journal of the American Statistical Association, Vol. 47, No. 260, pp. 635-646.

Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.
}

\seealso{
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
# Creation of a design object:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Estimate of the deciles of the income variable for
# the whole population:
svystatQ(des,~income,probs=seq(0.1,0.9,0.1),ties="rounded")


# Another design object:
data(sbs)
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
     fpc=~fpc)

# Estimation of the median value added 
# for economic activity macro-sectors:
svystatQ(des,~va.imp2,probs=0.5,by=~nace.macro,
         ties="rounded",vartype="cvpct")

# Estimation of the Interquartile Range (IQR) of the number
# of employees for economic activity macro-sectors:
apply(svystatQ(des,~emp.num,probs=c(0.25,0.75),by=~nace.macro)[,2:3],1,diff)
}

\keyword{survey}