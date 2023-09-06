\name{svystatB}
\alias{svystatB}
\alias{coef.svystatB}
\alias{coef.svystatB.by.list}
\alias{SE.svystatB}
\alias{SE.svystatB.by.list}
\alias{VAR.svystatB}
\alias{VAR.svystatB.by.list}
\alias{cv.svystatB}
\alias{cv.svystatB.by.list}
\alias{deff.svystatB}
\alias{deff.svystatB.by.list}
\alias{confint.svystatB}
\alias{confint.svystatB.by.list}
\alias{summary.svystatB}
\alias{summary.svystatB.by}
\alias{summary.svystatB.by.list}
\alias{print.svystatB.by.list}

\title{Estimation of Population Regression Coefficients in Subpopulations}

\description{
Computes estimates, standard errors and confidence intervals for Multiple Regression Coefficients in subpopulations.
}

\usage{
svystatB(design, model, by = NULL,
         vartype = c("se", "cv", "cvpct", "var"),
         conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
         na.rm = FALSE)

\method{coef}{svystatB}(object, ...)
\method{SE}{svystatB}(object, ...)
\method{VAR}{svystatB}(object, ...)
\method{cv}{svystatB}(object, ...)
\method{deff}{svystatB}(object, ...)
\method{confint}{svystatB}(object, ...)
\method{summary}{svystatB}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{model}{Formula specifying the linear model whose coefficients have to be estimated.}
  \item{by}{Formula specifying the variables that define the "estimation domains" (see \sQuote{Details}). If \code{NULL}
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
  \item{object}{An object of class \code{svystatB}.}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any).}
}

\details{
This function computes weighted estimates for Multiple Regression Coefficients using suitable weights depending on the class of \code{design}: calibrated weights for class \code{cal.analytic} and direct weights otherwise. Standard errors are calculated using the Taylor linearization technique.

The mandatory argument \code{model} identifies the regression model whose population coefficients have to be estimated (for details on model specification, see e.g. \code{\link{lm}}). The \code{design} variables referenced by \code{model} should be \code{numeric} or \code{factor} (variables of other types - e.g. \code{character} - will need to be converted in advance, e.g. using function \code{\link{des.addvars}}).

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svystatB} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{svystatQ} twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. Note that, to prevent obvious collinearity issues, the variables referenced by argument \code{by} must \emph{not} appear in the input \code{model} formula: otherwise, the program will stop and print an error message.

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

The optional argument \code{deff} allows to request the design effect [Kish 1995] for the estimates. By default \code{deff=FALSE}, that is the design effect is not provided. The design effect of an estimator is defined as the ratio between the variance of the estimator under the actual sampling design and the variance that would be obtained for an 'equivalent' estimator under a hypothetical simple random sampling without replacement of the same size. To obtain an estimate of the design effect comparing to simple random sampling \emph{\dQuote{with replacement}}, one must use \code{deff="replace"}.\cr
For nonlinear estimators, the design effect is estimated on the linearized version of the estimator (that is for the estimator of the total of the linearized variable, aka "Woodruff transform").\cr
When dealing with domain estimation, the design effects referring to a given subpopulation are currently computed by taking the ratios between the actual variance estimates and those that would have been obtained if a simple random sampling were carried out \emph{within} that subpopulation. This is the same as the \code{srssubpop} option for Stata's function \code{estat}.

Missing values (\code{NA}) in model variables should be avoided. If \code{na.rm=FALSE} (the default) they generate an error. If \code{na.rm=TRUE}, observations containing NAs in model variables are dropped, and estimates get computed on non missing values only. This implicitly assumes that missing values hit interest variables \emph{completely at random}: should this not be the case, computed estimates would be \emph{biased}.

The \code{summary} method invoked on regression coefficients (say \code{b}) estimated via \code{svystatB}, gives p-values and significance codes for the component-wise test \code{b = 0}. Such values are computed assuming that the distribution of the regression coefficients estimators is normal (which is asymptotically true for large scale surveys). This assumption has the advantage of overcoming the problem of choosing the "right" statistic and assessing its "right" number of degrees of freedom when using data from a complex survey (see e.g. [Korn, Graubard 1990]).
}

\section{Collinearity, Aliasing and Impacts in Domain Estimation}{
Function \code{svystatB} overcomes problems arising from exact collinearity between model variables via \emph{\sQuote{aliasing}} (see the \sQuote{Examples} Section). Put simply, aliasing discards redundant (i.e. collinear) regressors, yielding exact estimates and standard errors for non-aliased regression coefficients (namely, the same results that would be obtained with a reduced - no collinearity - model). Note that for the way aliasing works, the order of the terms in the linear model formula definitely matters.

Collinearity between variables may manifest itself in subsets of the sample, and with different patterns across subsets. In domain estimation, this phenomenon can have an impact on the structure of \code{svystatB}'s output. In fact, owing to aliasing, the estimable regression coefficients - for the same input linear model - can be different across domains. In such cases, the domain estimates produced by \code{svystatB} can no longer be stored in a \code{data.frame}, and the output object will instead be a \emph{list} (see the \sQuote{Examples} Section). Note that for these \code{svystatB}'s return objects (whose \code{class} is \code{svystatB.by.list}) no variance \code{\link{extractors}} are currently available.

Note also that, for the reasons above, the usage of \code{svystatB} via function \code{\link{svystat}} is restricted: it is not allowed to specify \code{\link{svystat}}'s arguments \code{by} and \code{group} when \code{kind == "B"}. 
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values. In special cases (see Section \sQuote{Collinearity, Aliasing and Impacts in Domain Estimation}), a \emph{list} object.
}

\author{
Diego Zardetto
}

\references{
Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.

Kish, L. (1995). \emph{\dQuote{Methods for design effects}}. Journal of Official Statistics, Vol. 11, pp. 55-77.

European Commission, Eurostat, (2013). \emph{\dQuote{Handbook on precision requirements and variance estimation for ESS households surveys: 2013 edition}}, Publications Office. doi: 10.2785/13579

Korn, E.L., Graubard, B.I. (1990) \emph{\dQuote{Simultaneous testing of regression coefficients with complex survey
data: Use of Bonferroni t statistics}}. The American Statistician, 44, 270-276.
}

\seealso{
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Quantiles \code{\link{svystatQ}}, Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}, and all of the above \code{\link{svystat}}.
}

\examples{
######################################################
# A simple regression model with a single predictor. #
# Let's compare the estimated regression coefficient #
# to its true value computed on the sampling frame.  #
######################################################

# Load sbs data:
data(sbs)

# Create a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
        fpc=~fpc)
 
# The population scatterplot of y vs emp.num reveals a linear
# behaviour:
plot(sbs.frame$emp.num,sbs.frame$y,
     col=rgb(50,205,50,100,maxColorValue=255),pch=16)

# Compute the population fit of the linear regression
# model y~emp.num-1 (no intercept):
pop.fit<-lm(y~emp.num-1,data=sbs.frame)
abline(pop.fit,col="red",lwd=2,lty=2)

# The obtained population R-squared is quite significant
# (greater than 0.7):
pop.R2<-summary(pop.fit)$r.squared
pop.R2

# The population regression coefficient is:
B<-coef(pop.fit)
B

# Now let's estimate B on the basis of the sbs sample and
# let's build a 95% confidence interval for the obtained estimate:
svystatB(sbsdes,y~emp.num-1,conf.int=TRUE)

# Thus, the confidence interval covers the true value of B.

# Notice that using ReGenesees Complex Estimators function
# svystatL, you would have obtained exactly the same results:
sbsdes<-des.addvars(sbsdes,y4emp.num=y*emp.num,
                           emp.num.sq=emp.num^2)
svystatL(sbsdes,expression(y4emp.num/emp.num.sq),
         conf.int=TRUE)


##################################
# A multiple regression example. #
##################################

# Let's estimate the coefficients of a model describing
# value added (variable va.imp2) as a linear function
# of number of employees by region and of nace.macro:
b <- svystatB(sbsdes,va.imp2~emp.num:region+nace.macro,vartype="cvpct")
b

# To obtain p-values and significance codes for the
# component-wise test t=0, you can exploit the
# summary method:
summary(b)

# Notice that estimators normality is assumed.


##########################################
# Obtaining domain means via regression. #
##########################################

# The domain mean of a numeric variable can be thought
# as a regression coefficient. Suppose you need the
# average number of employees by macro-sector, you can
# do as follows:
svystatB(sbsdes,emp.num~nace.macro-1)
 
# ...which, indeed, gives exactly the same results of:
svystatTM(sbsdes,y=~emp.num,by=~nace.macro,estimator="Mean")


##########################
# Handling collinearity. #
##########################

# Function svystatB overcomes problems arising from exact
# collinearity between model variables via 'aliasing'.
# To understand how aliasing works, let's build a manifestly
# redundant linear model:
svystatB(sbsdes,y~emp.num+I(2*emp.num)+I(3*va.imp2)+va.imp2-1)

# The obtained warning message shows that order definitely matters
# in aliasing, indeed:
svystatB(sbsdes,y~emp.num+I(2*emp.num)+va.imp2+I(3*va.imp2)-1)

# Notice also that aliasing gives exact estimates and standard errors
# for non-aliased regression coefficients (i.e. the same results that
# would be obtained with a reduced - no collinearity - model):
svystatB(sbsdes,y~emp.num+va.imp2-1)


###############################################
# Handling missing values in model variables. #
###############################################

# Load fpcdat:
data(fpcdat)

# Now, let's introduce some NAs in survey data:
fpcdat$y[c(1,3)]<-NA
fpcdat$x[c(3,5)]<-NA

# Create a design object:
fpcdes<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
        fpc=~fpc1+fpc2)

# Let's estimate regression coefficients of model z~y+x
  # na.rm=FALSE (the default) leads to an error:
  \dontrun{
  svystatB(fpcdes,z~y+x)
  }

  # whereas na.rm=TRUE simply drops all the cases
  # with missing data in model variables:  
  svystatB(fpcdes,z~y+x,na.rm=TRUE)


##################################
# Handling non positive weights. #
##################################

# Non positive direct weights are not allowed, anyway some
# calibrated weights can sometimes turn out to be <= 0. The
# corrisponding observations would be dropped by svystatB.

  # Prepare a template for population totals:
  pop<-pop.template(fpcdes,~z+pl.domain-1)

  # Fill it with silly values in order to obtain some negative g-weights:
  pop[1,]<-c(20000,90,10,90)

  # Calibrate:
  fpccal<-e.calibrate(fpcdes,pop)
  
  # We got 2 negative calibrated weights:
  g.range(fpccal)
  sum(weights(fpccal)<=0)

  # Now, let's estimate regression coefficients of model z~y+x
  # and pay attantion to the warnings:
  svystatB(fpccal,z~y+x,na.rm=TRUE)


####################################################################
# Domain estimates of simple and multiple regression coefficients. #
####################################################################

# Estimate the coefficients of the simple regression y ~ emp.num by domains
# obtained crossing region and nace.macro:
bb <- svystatB(sbsdes, model= va.imp2 ~emp.num, by= ~region:nace.macro)
bb

# Obtain p-values and significance codes via the summary method:
summary(bb)

# You have yet another method to estimate domain means of numeric variables.
# Suppose you need the average number of employees by macro-sector, you can
# do as follows:
svystatB(sbsdes, model= emp.num ~1, by= ~nace.macro)

# ...which gives exactly the same results of:
svystatB(sbsdes, model= emp.num ~nace.macro -1)

# ...and, of course, of:
svystatTM(sbsdes, y= ~emp.num, by= ~nace.macro, estimator= "Mean")


# One multiple regression example:
svystatB(sbsdes, model = y ~ va.imp2 + emp.num, by = ~region:nace.macro)


# A case of differential aliasing across domain (note the warning messages).
# A list-like object is returned:
svystatB(sbsdes, va.imp2 ~emp.num:emp.cl + nace.macro, by= ~region:public)

}
\keyword{survey}