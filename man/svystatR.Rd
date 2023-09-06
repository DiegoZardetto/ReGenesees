\name{svystatR}
\alias{svystatR}
\alias{coef.svystatR}
\alias{SE.svystatR}
\alias{VAR.svystatR}
\alias{cv.svystatR}
\alias{deff.svystatR}
\alias{confint.svystatR}

\title{Estimation of Ratios in Subpopulations}

\description{
Calculates estimates, standard errors and confidence intervals for Ratios between Totals in subpopulations.
}

\usage{
svystatR(design, num, den, by = NULL, cross = FALSE,
         vartype = c("se", "cv", "cvpct", "var"),
         conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
         na.rm = FALSE)

\method{coef}{svystatR}(object, ...)
\method{SE}{svystatR}(object, ...)
\method{VAR}{svystatR}(object, ...)
\method{cv}{svystatR}(object, ...)
\method{deff}{svystatR}(object, ...)
\method{confint}{svystatR}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{num}{Formula defining the numerator variables for the ratios.}
  \item{den}{Formula defining the denominator variables for the ratios.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{cross}{Should ratios be estimated for all the pairs of variables in \code{'num'} and \code{'den'}? The default is \code{FALSE}, meaning that ratios get estimated parallel-wise (see \sQuote{Details}).}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.int}{Compute confidence intervals for the estimates? The default is
                 \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{deff}{Should the design effect be computed? The default is \code{FALSE} (see \sQuote{Details}).}
  \item{na.rm}{Should missing values (if any) be removed from the variables of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
  \item{object}{An object of class \code{svystatR}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Ratios between Totals using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors are calculated using the Taylor linearization technique.

The mandatory argument \code{num} (\code{den}) identifies the variables whose totals appear as numerators (denominators) in the Ratios: the corresponding formula must be of the type \code{num = ~num.1 + \ldots + num.k} (\code{den = ~den.1 + \ldots + den.l}). The \code{design} variables referenced by \code{num} (\code{den}) must be \code{numeric}.

If \code{cross=TRUE}, the function computes estimates for \emph{all} the Ratios between pairs of variables coming from \code{num} and \code{den} (that is \code{k*l} estimates for the formulae above). If, on the contrary, \code{cross=FALSE} (the default), Ratios get estimated parallel-wise and \R recycling rule is applied whenever \code{k!=l}: for the formulae above, this generates \code{r} Ratios, where \code{r=max(k,l)}.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatR} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatR} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates.  By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the variance of the estimator under the actual sampling design and the variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}.\cr
Being Ratios nonlinear estimators, the design effect is estimated on the linearized version of the estimator (that is: for the estimator of the total of the linearized variable, aka "Woodruff transform").\cr
When dealing with domain estimation, the design effects referring to a given subpopulation are currently computed by taking the ratios between the actual variance estimates and those that would have been obtained if a simple random sampling were carried out \emph{within} that subpopulation. This is the same as the \code{srssubpop} option for Stata's function \code{estat}.

Missing values (\code{NA}) in interest variables should be avoided. If \code{na.rm=FALSE} (the default) they generate NAs in estimates (or even an error, if \code{design} is calibrated). If \code{na.rm=TRUE}, observations containing NAs are dropped, and estimates get computed on non missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random}: should this not be the case, computed estimates would be \emph{biased}. Notice that the \code{na.rm=TRUE} option is only allowed for a single Ratio, i.e. if \code{num} and \code{den} reference a \emph{single} interest variable.
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values.
}

\section{Warning}{
It can happen that, in some subpopulations, the estimate of the Total of some \code{den} variables turns out to be zero. In such cases \code{svystatR} estimates are either \code{NaN} or \code{Inf}, and \code{NaN} is returned for the corresponding SE estimates.
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
Estimators of Totals and Means \code{\link{svystatTM}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
# Creation of a design object:
data(sbs)
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
     fpc=~fpc)

# Estimation of the average value added per employee
# at the nation level:
svystatR(des,~va.imp2,~emp.num)

# The same as above by economic activity macro-sector:
svystatR(des,~va.imp2,~emp.num,~nace.macro,vartype="cvpct")


# Another design object:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Estimation of the ratios y1/x1, y1/x2, y2/x1 and y2/x2 by region,
# notice the use of argument cross:
svystatR(des,~y1+y2,~x1+x2,by=~regcod,cross=TRUE)

# ... compare the latter with the default (i.e. cross=FALSE)
svystatR(des,~y1+y2,~x1+x2,by=~regcod)


# Estimation of the ratios z/x1, z/x2 e z/x3
# for the whole population (notice the recycling rule):
svystatR(des,~z,~x1+x2+x3,conf.int=TRUE)

# Estimators of means can be thought as 
# estimators of ratios:
svystatTM(des,~income,estimator="Mean")  
svystatR(des.addvars(des,ones=1),num=~income,den=~ones)


##################################################
# Household-level averages in household surveys. #
##################################################

# For an introduction on this topic, see ?svystatTM examples.

  # Load survey data:
  data(data.examples)

  # Define the survey design (variable famcod identifies households) 
  exdes<-e.svydesign(data=example,ids=~towcod+famcod,strata=~stratum,
         weights=~weight)

  # Collapse strata to eliminate lonely PSUs
  exdes<-collapse.strata(design=exdes,block.vars=~sr:procod)

  # Now add new convenience variables to the design object:
    ## 'ones':       to estimate individuals counts
    ## 'housize':    to classify individuals by household size
    ## 'houdensity': to estimate households counts
  exdes<-des.addvars(exdes,
                     ones=1,
                     housize=factor(ave(famcod,famcod,FUN = length)),
                     houdensity=ave(famcod,famcod,FUN = function(x) 1/length(x))
                    )

  # Estimate the average number of household components by region:
  svystatR(exdes,num=~ones,den=~houdensity,by=~regcod,
           vartype="cvpct",conf.int=TRUE)

  # Estimate the average household income for the whole population:
  svystatR(exdes,num=~income,den=~houdensity,vartype="cvpct",
           conf.int=TRUE)

  # ...for household size categories:
  svystatR(exdes,num=~income,den=~houdensity,by=~housize,
           vartype="cvpct",conf.int=TRUE)

  # ...and for province and household size:
  svystatR(exdes,num=~income,den=~houdensity,by=~housize:procod,
           vartype="cvpct")
}
\keyword{survey}