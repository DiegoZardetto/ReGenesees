\name{svystatTM}
\alias{svystatTM}
\alias{coef.svystatTM}
\alias{SE.svystatTM}
\alias{VAR.svystatTM}
\alias{cv.svystatTM}
\alias{deff.svystatTM}
\alias{confint.svystatTM}

\title{Estimation of Totals and Means in Subpopulations}

\description{
Computes estimates, standard errors and confidence intervals for Totals and Means in subpopulations.
}

\usage{
svystatTM(design, y, by = NULL, estimator = c("Total", "Mean"),
          vartype = c("se", "cv", "cvpct", "var"),
          conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
          na.rm = FALSE)

\method{coef}{svystatTM}(object, ...)
\method{SE}{svystatTM}(object, ...)
\method{VAR}{svystatTM}(object, ...)
\method{cv}{svystatTM}(object, ...)
\method{deff}{svystatTM}(object, ...)
\method{confint}{svystatTM}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{y}{Formula defining the variables of interest. Only \code{numeric} or \code{factor} variables are allowed.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{estimator}{\code{character} specifying the desired estimator: it may be \code{'Total'} (the default) or
                  \code{'Mean'}.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.int}{Compute confidence intervals for the estimates? The default is
                 \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{deff}{Should the design effect be computed? The default is \code{FALSE} (see \sQuote{Details}).}
  \item{na.rm}{Should missing values (if any) be removed from the variables of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
  \item{object}{An object of class \code{svystatTM}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Totals and Means using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors for nonlinear estimators (e.g. calibration estimators) are calculated using the Taylor linearization technique.

The mandatory argument \code{y} identifies the variables of interest, that is the variables for which estimates are to be calculated. The corresponding formula should be of the type \code{y=~var1+\ldots+varn}. The \code{design} variables referenced by \code{y} should be \code{numeric} or \code{factor} (variables of other types - e.g. \code{character} - will generate and error). It is admissible to specify for \code{y} "mixed" formulae that simultaneously contain quantitative (\code{numeric}) variables and qualitative (\code{factor}) variables.

To override the restriction to formulae of the type \code{y=~var1+\ldots+varn}, the AsIs operator \code{I()} can be used (see \sQuote{Examples}). Though the latter opportunity could appear quite useful in some occasion, actually it should be almost always possible to find a work-around by using other functions of the \pkg{ReGenesees} package.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatTM} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatTM} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

The optional argument \code{estimator} makes it possible to select the desired estimator. If\cr \code{estimator="Total"} (the default option), \code{svystatTM} calculates, for a given variable of interest \code{vark}, the estimate of the total (when \code{vark} is \code{numeric}) or the estimate of the absolute frequency distribution (when \code{vark} is \code{factor}). Similarly, if \code{estimator="Mean"}, the function calculates the estimate of the mean (when \code{vark} is \code{numeric}) or the the estimate of the relative frequency distribution (when \code{vark} is \code{factor}).

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates. By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the variance of the estimator under the actual sampling design and the variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}.\cr
Understanding what 'equivalent' estimator actually means is straightforward when dealing with Horvitz-Thompson estimators of Totals and Means. This is not the case when, on the contrary, the estimator to which the deff refers is a nonlinear estimator (e.g. for Calibration estimators of Totals and Means). In such cases, the standard approach is to use as 'equivalent' estimator the linearized version of the original estimator (that is: the estimator of the total of the linearized variable, aka "Woodruff transform").\cr
When dealing with domain estimation, the design effects referring to a given subpopulation are currently computed by taking the ratios between the actual variance estimates and those that would have been obtained if a simple random sampling were carried out \emph{within} that subpopulation. This is the same as the \code{srssubpop} option for Stata's function \code{estat}.

Missing values (\code{NA}) in interest variables should be avoided. If \code{na.rm=FALSE} (the default) they generate NAs in estimates (or even an error, if \code{design} is calibrated). If \code{na.rm=TRUE} and \code{estimator="Mean"}, observations containing NAs are dropped, and estimates get computed on non-missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random} (MCAR): should this not be the case, computed estimates would be \emph{biased}. Since, even under the MCAR assumption, estimates of totals and counts solely based on complete cases would be biased, function \code{svystatTM} adopts a model-based ratio estimator when \code{na.rm=TRUE} and \code{estimator="Total"}. This is obtained by (i) first estimating the mean or proportion based on complete cases only, and (ii) then multiplying the result by the estimated population size based on all observations (i.e. both missing and non-missing). This model-based ratio estimator is asymptotically design-unbiased under the MCAR assumption. Notice that the \code{na.rm=TRUE} option is only allowed if \code{y} references a \emph{single} interest variable.
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
Estimators of Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Quantiles \code{\link{svystatQ}}, Multiple Regression Coefficients \code{\link{svystatB}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
# Load survey data:
data(data.examples)

# Creation of a design object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)


# Estimation of the total of 3 quantitative variables for the whole
# population:
svystatTM(des,~y1+y2+y3)


# Estimation of the total of the same 3 variables by region, with SE
# and CV%: 
svystatTM(des,~y1+y2+y3,~regcod,vartype=c("se","cvpct"))


# Estimation of the mean of the same 3 variables by marstat and sex:
svystatTM(des,~y1+y2+y3,~marstat:sex,estimator="Mean")


# Estimation of the absolute frequency distribution of the qualitative
# variable age5c for the whole population, with the design effect:
svystatTM(des,~age5c,deff=TRUE)

# MARGINAL relative frequency distributions
  # Estimation of the relative frequency distribution of the qualitative
  # variable age5c for the whole population:
  svystatTM(des,~age5c,estimator="Mean")

# CONDITIONAL relative frequency distributions
  # Estimation of the relative frequency distribution of the qualitative
  # variable marstat by sex:
  svystatTM(des,~marstat,~sex,estimator="Mean")

# JOINT relative frequency distributions
  # Estimation of the relative frequency of the joint distribution of sex
  # and marstat:
    # *First Solution* (using the AsIs operator I()):
    svystatTM(des,~I(sex:marstat),estimator="Mean")
    # *Second Solution* (adding a new variable to des):
    des2 <- des.addvars(des, sex.marstat=sex:marstat)
    svystatTM(des2,~sex.marstat,estimator="Mean")
    # *Third Solution* (exploiting estimators of Shares, see also ?svystatS):
    # Add new variable 'ones' to estimate counts of final units (individuals)
    # and estimate the share of people for classes of sex and marstat
    des2 <- des.addvars(des, ones=1)
    svystatS(des2,~ones,classes=~sex:marstat)


# Estimation of the mean income inside provinces, with confidence intervals
# at a confidence level of 0.9:
svystatTM(des,~income,~procod,estimator="Mean",conf.int=TRUE,conf.lev=0.9)


# Quantitative and qualitative variables together: estimation of the
# total of income and of the absolute frequency distribution of sex,
# by marstat:
svystatTM(des,~income+sex,~marstat)


# Estimating totals in domains for "incomplete" partitions: more on
# the AsIs operator I()

  # Estimation of the total income (plus cvpct) ONLY in region 7:
  svystatTM(des,~I(income*(regcod=="7")),vartype="cvpct")
  # Alternative solution (adding a new variable to des):
  des2 <- des.addvars(des, inc_reg7=I(income*(regcod=="7")))
  svystatTM(des2,~inc_reg7,vartype="cvpct")

  # Estimation of the total income (plus cvpct) ONLY in regions 6 and 10:
  svystatTM(des,~I(income*as.numeric(regcod \%in\% c("6","10"))),vartype="cvpct")
  # Alternative solution (adding a new variable to des):
  des2 <- des.addvars(des, inc_reg6.10=I(income*(regcod \%in\% c("6","10"))))
  svystatTM(des2,~inc_reg6.10,vartype="cvpct")

  # Compare with the corresponding estimates for the "complete" partition,
  # i.e. for regions:
  svystatTM(des,~income,~regcod,vartype="cvpct")


# Under default settings lonely PSUs produce errors in standard
# errors estimation (notice we didn't pass the fpcs):
data(fpcdat)
des.lpsu<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,
          weights=~w)
\dontrun{
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
}

# This can be circumvented in different ways, namely:
old.op <- options("RG.lonely.psu"="adjust")
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
options(old.op)

# or otherwise:
old.op <- options("RG.lonely.psu"="average")
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
options(old.op)

# but see also ?collapse.strata for a better alternative.


####################################################
# Household-level estimation in household surveys. #
####################################################

# Large scale household surveys typically adopt a 2-stage sampling
# design with municipalities as PSUs and households as SSUs, in order
# to eventually collect information on each individual belonging to
# sampled SSUs. In such a framework (up to possible total nonresponse
# effects), each individual inside a sampled household shares the
# same direct weight, which, in turn, equals the household weight.
# This implies that it is very easy to build estimates referred to
# SSU-level (households) information, despite estimators actually
# involve only individual values. Some examples are given below.

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

  # Estimate the total number of households:
  nhou<-svystatTM(exdes,~houdensity,vartype="cvpct")
  nhou

  # Estimate the total number of individuals:
  nind<-svystatTM(exdes,~ones,vartype="cvpct")
  nind

  # Thus the average number of individuals per household is:
  coef(nind)/coef(nhou)

  # ...which can be obtained also as a ratio (along with
  # its estimated sampling variability):
  svystatR(exdes,~ones,~houdensity,vartype="cvpct")

  # Estimate the number and proportion of individuals living in households
  # of given sizes:
  nind.by.housize<-svystatTM(exdes,~housize,vartype="cvpct")
  nind.by.housize
  
  pind.by.housize<-svystatTM(exdes,~housize,estimator="Mean",var="cvpct")
  pind.by.housize

  # Estimate the number of households by household size:
  nhou.by.housize<-svystatTM(exdes,~houdensity,~housize,vartype="cvpct")
  nhou.by.housize

  # Notice that estimates of individuals and household counts are consistent,
  # indeed:
  coef(nind.by.housize)/coef(nhou.by.housize)
}
\keyword{survey}