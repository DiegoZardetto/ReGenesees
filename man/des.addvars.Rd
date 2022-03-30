\name{des.addvars}
\alias{des.addvars}

\title{Add Variables to Design Objects}

\description{
Modifies an analytic object by adding new variables to it. 
}

\usage{
des.addvars(design, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{\ldots}{\code{tag = expr} arguments defining columns to be added to \code{design}.}
}

\details{
This function adds to the data frame contained in \code{design} the \emph{new} variables defined by the \code{tag = expr} arguments. A \code{tag} can be specified either by means of an identifier or by a character string; \code{expr} can be any expression that it makes sense to evaluate in the \code{design} environment.

For each argument \code{tag = expr} bound to the formal argument \code{\ldots} the added column will have \emph{name} given by the \code{tag} value and \emph{values} obtained by evaluating the \code{expr} expression on \code{design}. Any input expression not supplied with a \code{tag} will be ignored and will therefore have no effect on the \code{des.addvars} return value.

Variables to be added to the input object have to be \emph{new}: namely it is not possible to use \code{des.addvars} to modify the values in a pre-existing \code{design} column. This an intentional feature meant to safeguard the integrity of the relations between survey data and sampling design metadata stored in \code{design}.
}

\value{
An object of the same class of \code{design}, containing new variables but supplied with exactly the same metadata.
}

\author{Diego Zardetto}

\references{
Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata, \code{\link{e.calibrate}} for calibrating weights.
}

\examples{
data(data.examples)

# Creation of an analytic object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Adding the new 'ones' variable to estimate the number
# of final units in the population: 
des<-des.addvars(des,ones=1)
svystatTM(des,~ones)

# Recoding a qualitative variable:
des<-des.addvars(des,agerange=factor(ifelse(age5c==1,
                                     "young","not-young")))
svystatTM(des,~agerange,estimator="Mean")
svystatTM(des,~income,~agerange,estimator="Mean",conf.int=TRUE)

# Algebraic operations on numeric variables:
des<-des.addvars(des,z2=z^2)
svystatTM(des,~z2,estimator="Mean")

# A more interesting example: estimating the
# percentage of population with income below
# the poverty threshold (defined as 0.6 times
# the median income for the whole population):
Median.Income <- coef(svystatQ(des, ~income,probs=0.5))
Median.Income
des <- des.addvars(des,
                   status = factor(
                          ifelse(income < (0.6 * Median.Income),
                          "poor", "non-poor")
                                  )
                   )
svystatTM(des,~status,estimator="Mean")
# Mean income for poor and non-poor:
svystatTM(des,~income,~status,estimator="Mean")

### NOTE: Procedure above yields *correct point estimates* of the share of poor
###       population and their average income, while *variance estimation is
###       approximated* since we neglected the sampling variability of the
###       estimated poverty threshold.
}
\keyword{survey}