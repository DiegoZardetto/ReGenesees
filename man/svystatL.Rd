\name{svystatL}
\alias{svystatL}
\alias{coef.svystatL}
\alias{SE.svystatL}
\alias{VAR.svystatL}
\alias{cv.svystatL}
\alias{deff.svystatL}
\alias{confint.svystatL}

\title{Estimation of Complex Estimators in Subpopulations}

\description{
Computes estimates, standard errors and confidence intervals for Complex Estimators in subpopulations. A Complex Estimator can be any analytic function of (Horvitz-Thompson or Calibration) estimators of Totals and Means.
}

\usage{
svystatL(design, expr, by = NULL,
         vartype = c("se", "cv", "cvpct", "var"),
         conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
         na.rm = FALSE)

\method{coef}{svystatL}(object, ...)
\method{SE}{svystatL}(object, ...)
\method{VAR}{svystatL}(object, ...)
\method{cv}{svystatL}(object, ...)
\method{deff}{svystatL}(object, ...)
\method{confint}{svystatL}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{expr}{R expression defining the Complex Estimator (see \sQuote{Details}).}
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
  \item{object}{An object of class \code{svystatL}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Complex Estimators using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors are calculated using the Taylor linearization technique.

Function \code{svystatL} can handle any user-defined estimator that can be expressed as an analytic function of Horvitz-Thompson or Calibration estimators of Totals or Means, by automatically linearizing them (see [Zardetto, 15] for details).

The mandatory argument \code{expr}, which identifies the Complex Estimator, must be an object of class \code{expression}. It can be specified just a single Complex Estimator at a time, i.e. \code{length(expr)} must be equal to \code{1}. Any analytic function of estimators of Totals and Means is allowed.\cr
Inside \code{expr} the estimator of the Total of a variable is simply represented by the \emph{name} of the variable itself. To represent the estimator of the Mean of a variable \code{y}, the expression \code{y/ones} has to be used (\code{ones} being the convenience name of an artificial variable - created on-the-fly - whose value is \code{1} for each elementary sampling unit, so that its Total estimator actually estimates the size of the population). Variables referenced inside \code{expr} must obviously belong to \code{design} and must be \code{numeric}.\cr
At a minimal level, \code{svystatL} can be used to estimate Totals, Means, Ratios, etc., thus reproducing the same results achieved by using the corresponding dedicated functions \code{svystatTM}, \code{svystatR}, etc. For instance, calling \code{svystatL(design, expression(y/x))} is equivalent to invoking \code{svystatR(design, ~y, ~x)}, while using \code{svystatL(design, expression(y/ones))} or\cr
\code{svystatTM(design, ~y, estimator = "Mean")} achieves an identical result.

The mathematical expression of a Complex Estimator, as specified by argument \code{expr}, can involve \emph{'parameters'}, that is \emph{symbols} representing given, non-random, scalar, numeric \emph{values}. For each parameter appearing in \code{expr}, the value corresponding to its symbol will be searched following \R standard scoping rules, see e.g. the first example in Section \sQuote{Examples} for a practical illustration.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatL} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatL} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates. By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the variance of the estimator under the actual sampling design and the variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}.\cr
For nonlinear estimators, the design effect is estimated on the linearized version of the estimator (that is for the estimator of the total of the linearized variable, aka "Woodruff transform").\cr
When dealing with domain estimation, the design effects referring to a given subpopulation are currently computed by taking the ratios between the actual variance estimates and those that would have been obtained if a simple random sampling were carried out \emph{within} that subpopulation. This is the same as the \code{srssubpop} option for Stata's function \code{estat}.

Missing values (\code{NA}) in interest variables should be avoided. If \code{na.rm=FALSE} (the default) they generate NAs in estimates (or even an error, if \code{design} is calibrated). If \code{na.rm=TRUE}, observations containing NAs are dropped, and estimates get computed on non missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random}: should this not be the case, computed estimates would be \emph{biased}.
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values.
}

\section{Warning}{
When the linearized variable corresponding to a Complex Estimator is ill defined (because the estimator gradient is singular at the Taylor series expansion point), SE estimates returned by \code{svystatL} are \code{NaN}.
}

\author{
Diego Zardetto
}

\references{
Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.

Kish, L. (1995). \emph{\dQuote{Methods for design effects}}. Journal of Official Statistics, Vol. 11, pp. 55-77.

European Commission, Eurostat, (2013). \emph{\dQuote{Handbook on precision requirements and variance estimation for ESS households surveys: 2013 edition}}, Publications Office. doi: 10.2785/13579
}

\seealso{
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, and all of the above \code{\link{svystat}}.
}

\examples{
#####################################################
# A first example: the Ratio Estimator of a Total.  #
#####################################################

# Creation of a design object:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Recall that ratio estimators of Totals rely on auxiliary
# information. Thus, suppose you want to estimate the total
# of income and suppose you know from an external source that
# the population size is, say, 1E6:
pop <- 1E6

# To obtain the ratio estimator of total income, you can do as follows:
  ## A) Directly plug the numeric value of pop into expr
svystatL(des, expression(1E6 * (income/ones)), vartype = "cvpct")

  ## B) Treat pop as a parameter and let R find its actual value (1E6) inside
  ##    the calling environment of svystatL (the .GlobalEnv)
svystatL(des, expression(pop * (income/ones)), vartype = "cvpct")

# NOTE: Method B) can be very useful for simulation purposes, as it avoids
#       having to directly type in numbers when invoking svystatL (something
#       that only a human in an interactive R session could do).

# By comparing the latter result with the ordinary estimator of the mean...
svystatTM(des,~income,vartype="cvpct")

# ...one can appreciate the variance reduction stemming from the correlation
# between numerator and denominator:
Corr(des, expression(income), expression(ones))


############################################################
# A complex example: estimation of the Population Standard #
# Deviation of a variable.                                 #
############################################################

# Creation of another design object:
data(sbs)
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
        fpc=~fpc)
 
# Suppose you want to estimate the standard deviation of the
# population distribution of value added (va.imp2):
sbsdes<-des.addvars(sbsdes,va.imp2.sq=va.imp2^2)
svystatL(sbsdes,expression( sqrt( (ones/(ones-1))*
                                  ( (va.imp2.sq/ones)-(va.imp2/ones)^2 )
                                )
                          ), conf.int=TRUE)

# The estimate above and the associated confidence interval (which
# involves the estimate of the sampling variance of the complex
# estimator) turn out to be very sound: indeed the TRUE value of the
# parameter is:
sd(sbs.frame$va.imp2)


###############################################
# Estimation of Geometric and Harmonic Means. #
###############################################

## 1. Harmonic Mean
 # Recall that the the harmonic mean of a positive variable,
 # say z, can be computed as 1/mean(1/z). Thus, for instance,
 # to get a survey estimate of the harmonic mean of emp.num,
 # you can do as follows:
sbsdes<-des.addvars(sbsdes,emp.num.m1=1/emp.num)
h<-svystatL(sbsdes,expression( ones/emp.num.m1 ),
            conf.int=TRUE)
h

 # You can easily verify that the obtained estimate is close
 # to the true value (as computed from the sampling frame) and
 # covered by the 95% confidence interval:
1/mean(1/sbs.frame$emp.num)

## 2. Geometric Mean
 # Recall that the the geometric mean of a non negative variable,
 # say z, can be computed as exp(mean(log(z))). Thus, for instance,
 # to get a survey estimate of the geometric mean of emp.num,
 # you can do as follows:
sbsdes<-des.addvars(sbsdes,log.emp.num=log(emp.num))
g<-svystatL(sbsdes,expression( exp(log.emp.num/ones) ),
            conf.int=TRUE)
g

 # You can easily verify that the obtained estimate is close
 # to the true value (as computed from the sampling frame) and
 # covered by the 95% confidence interval:
exp(mean(log(sbs.frame$emp.num)))

## 3. Comparison with the arithmetic mean
 # If you compute the arithmetic mean estimate:
a<-svystatTM(sbsdes,~emp.num,estimator="Mean")
a

 #...you easily verify the expected hierachy,
 # i.e. harmonic <= geometric <= arithmetic:
H<-coef(h)
G<-coef(g)
A<-coef(a)
stopifnot(H <= G && G <= A)


#################################################################
# Further complex examples: estimation of Population Regression #
# Coefficients (for a model with a single predictor).           #
#################################################################

# Suppose you want to estimate of the slope of the population
# regression y vs. emp.num. You can do as follows:

  ## 1. No intercept model: y ~ emp.num - 1
   # Get survey estimate:
   sbsdes<-des.addvars(sbsdes,y4emp.num=y*emp.num,
                              emp.num.sq=emp.num^2)
   svystatL(sbsdes,expression(y4emp.num/emp.num.sq),
            conf.int=TRUE)

   # Compare with the actual slope from the population fit:
   pop.fit<-lm(y~emp.num-1,data=sbs.frame)
   coef(pop.fit)

   # ...a very good agreement.

  ## 2. The model with intercept: y ~ emp.num
   # Get survey estimate:
   svystatL(sbsdes,expression(  (ones*y4emp.num  - y*emp.num)/
                                (ones*emp.num.sq - emp.num^2)
                             ),
            conf.int=TRUE)
   
   # Compare with the actual slope from the population fit:
   pop.fit<-lm(y~emp.num,data=sbs.frame)
   coef(pop.fit)

   # ...again a very good agreement.

# Notice that both results above could be obtained also
# by using ReGenesees specialized function svystatB:

  ## 1.
  svystatB(sbsdes,y~emp.num-1,conf.int=TRUE)

  ## 2.
  svystatB(sbsdes,y~emp.num,conf.int=TRUE)

# Notice also - incidentally - that the estimate of the intercept
# turns out to be less accurate than the one we obtained for the slope,
# with about a 6% overestimation.
}
\keyword{survey}