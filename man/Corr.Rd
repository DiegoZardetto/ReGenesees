\name{Corr}
\alias{Corr}
\alias{CoV}

\title{Design Covariance and Correlation of Complex Estimators in Subpopulations}

\description{
Estimates the covariance and the correlation of Complex Estimators in subpopulations. A Complex Estimator can be any analytic function of (Horvitz-Thompson or Calibration) estimators of Totals and Means.
}

\usage{
CoV(design, expr1, expr2,
     by = NULL, na.rm = FALSE)
Corr(design, expr1, expr2,
     by = NULL, na.rm = FALSE)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{expr1}{R expression defining the first Complex Estimator (see \sQuote{Details}).}
  \item{expr2}{R expression defining the second Complex Estimator (see \sQuote{Details}).}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{na.rm}{Should missing values (if any) be removed from the variables of interest? The default is
              \code{FALSE} (see \sQuote{Details}).}
}

\details{
This function allows to estimate the covariance and the correlation of two arbitrary Complex Estimators. Estimates are calculated using the Taylor linearization technique.

The mandatory arguments \code{expr1} and \code{expr2} identify the Complex Estimators: both must be of class \code{expression}. For further details on the syntax and the semantics of such expressions, see \code{\link{svystatL}}.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{CoV} (\code{Corr}) refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}: if you need to compute estimates for domains \code{B1} and \code{B2} \emph{separately}, you have to call \code{CoV} (\code{Corr}) twice. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. 

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
}

\seealso{
Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, and Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}.
}

\examples{
########################################
# Some checks and some simple examples #
# to illustrate the syntax.            # 
########################################
# Load survey data:
data(data.examples)

# Creation of a design object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Let's start with some natural checks:
  ## The covariance of any estimator with itself is its variance
  ## (use mean income as an example):
     CoV(des,expression(income/ones),expression(income/ones))
     VAR(svystatL(des,expression(income/ones)))
     VAR(svystatTM(des,~income,estimator="Mean"))

  ## The correlation of any estimator with itself is 1
  ## (use mean income as an example):
     Corr(des,expression(income/ones),expression(income/ones))

# Switch to non trivial examples:
  ## Correlation of mean income with population size:
     Corr(des,expression(income/ones),expression(ones))

  ## Correlation of mean income with total income:
     # at population level:
     Corr(des,expression(income/ones),expression(income))
     # for regions:
     Corr(des,expression(income/ones),expression(income),by=~regcod)


  ## Correlation of a product of two totals and a ratio of two totals:
     # at population level:
     Corr(des,expression(y1*y2),expression(x1/x2))
     # for provinces:
     Corr(des,expression(income/ones),expression(income),by=~procod)


######################################################
# A more meaningful and complex example: correlation #
# between Geometric, Harmonic and Arithmetic Means.  #
######################################################
# Creation of another design object:
data(sbs)
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
        fpc=~fpc)

# Let's use variable emp.num, which is ok as it is always strictly positive:
  ## Add a convenience variable for estimating the harmonic mean (see ?svystatL
  ## for details) and prepare the formal estimator expression:
     sbsdes<-des.addvars(sbsdes,emp.num.m1=1/emp.num)
     h<-expression(ones/emp.num.m1)  
  
  ## Add a convenience variable for estimating the geometric mean (see ?svystatL
  ## for details) and prepare the formal estimator expression:
     sbsdes<-des.addvars(sbsdes,log.emp.num=log(emp.num))
     g<-expression(exp(log.emp.num/ones))

  ## prepare the formal estimator expression for the arithmetic mean:
     m<-expression(emp.num/ones)

# Now compute correlations:
  ## Harmonic with Arithmetic:
     Corr(sbsdes,h,m)

  ## Geometric with Arithmetic:
     Corr(sbsdes,g,m)

  ## Harmonic with Geometric:
     Corr(sbsdes,h,g)

  ## Hence, while correlations g-m and g-h are high, correlation h-m is low.


####################################################
# Another example: is a ratio estimator of a total #
# expected to be more efficient than an HT one?    #
####################################################

# Let's recall that the ratio estimator of a total is
# expected to be more efficient than HT, if the
# correlation of numerator and denominator exceeds
# half of the ratio between the CVs of denominator
# and numerator.

# Compute the HT estimate of the total of value added (variable va.imp2):
VA<-svystatTM(sbsdes,~va.imp2)
VA

# Compute the HT estimate of the total of emp.num:
EMP<-svystatTM(sbsdes,~emp.num)
EMP

# Now estimate the correlation of the numerator
# and denominator totals:
corr <- Corr(sbsdes,expression(va.imp2),expression(emp.num))
corr

# and compare it with (1/2)*( CV(den)/CV(num) )
stopifnot( corr > 0.5*cv(EMP)/cv(VA) )

# As the comparison holds TRUE, we expect an efficiency gain
# of the ratio estimator of the total compared to HT.
# Let's check...:

  # Compute the ratio estimate of the total of value added using
  # as auxiliary variable the number of emloyees, whose total
  # is 984394:
  TOT.emp.num <- sum(sbs.frame$emp.num)
  TOT.emp.num

  VA.ratio<-svystatL(sbsdes, expression(TOT.emp.num * (va.imp2/emp.num)))
  VA.ratio

  # Compare standard errors sizes:
  SE(VA.ratio)
  SE(VA)
  stopifnot( SE(VA.ratio) < SE(VA) )

# ...as expected.
}
\keyword{survey}