\name{population.check}
\alias{population.check}

\title{Compliance Test for Known Totals Data Frames}

\description{
Checks whether a known population totals data frame conforms to the standard required by \code{e.calibrate} for a specific calibration problem. 
}

\usage{
population.check(df.population, data, calmodel, partition = FALSE)
}

\arguments{
  \item{df.population}{Data frame of known population totals.}
  \item{data}{Data frame of survey data (or an object inheriting from class \code{analytic}).}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model. \code{FALSE} (the default) implies no calibration domains.}
}

\details{
The behaviour of this function depends on the outcome of the test. If \code{df.population} is found to conform to the standard, the function first converts it into an object of class \code{pop.totals} and then invisibly returns it. Failing this, the function stops and prints an error message: the meaning of the message should help the user diagnose the cause of the problem.
 
The mandatory argument \code{df.population} identifies the known totals data frame for which compliance with the standard is to be checked.

The mandatory argument \code{data} identifies the survey data frame on which the calibration problem is defined (or, as an alternative, an \code{analytic} object built upon that data frame).

The mandatory argument \code{calmodel} symbolically defines the calibration model you intend to use: it identifies the auxiliary variables and the constraints for the calibration problem. The \code{data} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}).

The optional argument \code{partition} specifies the variables that define the calibration domains for the model. The default value (\code{FALSE}) means either that there are not calibration domains or that you want to solve the problem globally (even though it could be factorized). If a formula is passed through the \code{partition} argument the program checks that \code{calmodel} actually describes a "reduced model", that is it does not reference any of the partition variables; if this is not the case, the program stops and prints an error message. Notice that a formula like \code{by=~D1+D2} will be automatically translated into the factor-crossing formula \code{by=~D1:D2}. The \code{data} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}). Note that, if the \code{partition} formula involves two or more factors, their crossed levels will be ordered according to operator \code{:} (that is, those from the \emph{rightmost} variable will vary fastest).
}

\value{
An invisible object of class \code{pop.totals}. The \code{pop.totals} class is a specialization of the \code{data.frame} class; this means that an object built by \code{pop.template} inherits from the \code{data.frame} class and you can use on it every method defined on that class.
}

\note{
The \code{population.check} function can be used to convert a known totals data frame that conforms to the standard required by \code{e.calibrate} into an object of class \code{pop.totals}. The usefulness of this conversion lies in the fact that, once you have known totals with this "certified format", you can invoke \code{e.calibrate} without specifying the values for the \code{calmodel} and \code{partition} arguments (this means that the function is able to extract them directly from the attributes of the \code{pop.totals} object). 
}

\author{Diego Zardetto}

\references{
Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\seealso{
\code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for the definition of the class \code{pop.totals} and to build a "template" data frame for known population totals, \code{\link{fill.template}} to automatically fill the template when a sampling frame is available.
}

\examples{
data(data.examples)

# Suppose you have to calibrate the example survey data frame
# on the totals of x1 by sex and you want the partitioned solution.
# Start creating a design object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Then build a template data frame for the known totals:
pop<-pop.template(data=example,calmodel=~x1-1,partition=~sex)
pop
class(pop)

# Now fill NAs with the actual values for the population
# totals (suppose 123 for sex="f" and 456 for sex="m"): 
pop[,"x1"]<-c(123,456)
pop
class(pop)

# Finally check if pop complies with the e.calibrate standard:
population.check(df.population=pop,data=example,calmodel=~x1-1,
                 partition=~sex)


# If, despite keeping the content unchanged, we altered the
# structure of the data frame (for example, by changing the 
# order of its rows)...
pop.mod<-pop ; pop.mod[1,]<-pop[2,] ; pop.mod[2,]<-pop[1,]
pop
pop.mod

# ...we would obtain an error:
\dontrun{
population.check(df.population=pop.mod,data=example,calmodel=~x1-1,
                 partition=~sex)
}


# Remember that, if the known totals have been converted
# into the pop.totals "format" by means of population.check,
# it is possible to invoke e.calibrate without specifying
# calmodel and partition: 

class(pop04p)
pop04p
descal04p<-e.calibrate(design=des,df.population=pop04p,
           calfun="logit",bounds=bounds,aggregate.stage=2)


# ...this option is not allowed if the known totals
# are not of class 'pop.totals' even if they conform to the
# standard:

pop04p.mod<-data.frame(pop04p)
class(pop04p.mod)
pop04p.mod
\dontrun{
e.calibrate(design=des,df.population=pop04p.mod,calfun="logit",
             bounds=bounds,aggregate.stage=2)
}
}
\keyword{survey}