\name{bounds.hint}
\alias{bounds.hint}

\title{A Hint for Range Restricted Calibration}

\description{
Suggests a sound \code{bounds} value for which \code{e.calibrate} is likely to converge.
}

\usage{
bounds.hint(design, df.population, 
    calmodel = if (inherits(df.population, "pop.totals"))
                   attr(df.population, "calmodel"), 
    partition = if (inherits(df.population, "pop.totals")) 
                    attr(df.population, "partition") else FALSE,
    msg = TRUE)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{df.population}{Data frame containing the known population totals for the auxiliary variables.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model; \code{FALSE} (the default) implies no calibration domains.}
  \item{msg}{Enables printing of a summary description of the result (the default is \code{TRUE}).}
}

\details{
The function \code{bounds.hint} returns a \code{bounds} value for which \code{e.calibtrate} is \emph{likely} to converge. This interval is just a sound hint, \emph{not} an exact result (see \sQuote{Note}).

The mandatory argument \code{design} identifies the \code{analytic} object on which the calibration problem is defined.

The mandatory argument \code{df.population} identifies the known totals data frame.

The argument \code{calmodel} symbolically defines the calibration model you want to use: it identifies the auxiliary variables and the constraints for the calibration problem. The \code{design} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}). The argument can be omitted provided \code{df.population} is an object of class \code{pop.totals} (see \code{\link{population.check}}). 

The optional argument \code{partition} specifies the variables that define the calibration domains for the model. The default value (\code{FALSE}) means either that there are not calibration domains or that you want to solve the problem globally (even though it could be factorized). The \code{design} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}). The argument can be omitted provided \code{df.population} is an object of class \code{pop.totals} (see \code{\link{population.check}}).

The optional argument \code{msg} enables/disables printing of a summary description of the achieved result.
}

\value{
A numeric vector of length 2, representing the \emph{suggested} value for the \code{bounds} argument of \code{e.calibrate}. The attributes of that vector store additional information, which can lead to better understand why a given calibration problem is (un)feasible (see \sQuote{Examples}).  
}

\note{
Assessing the feasibility of an arbitrary calibration problem is not an easy task. The problem is even more difficult whenever additional \emph{"range restrictions"} are imposed. Indeed, even if one assumes that the calibration constraints define a consistent system, one also has to choose the \code{bounds} such that the feasible region is non-empty.

One can argue that there must exist a minimun-length interval \eqn{I=[L,U]} such that, if it is covered by \code{bounds}, the specified calibration problem is feasible. Unfortunately in order to compute exactly that minimun-length interval \eqn{I} one should solve a big linear programming problem [Vanderhoeft 01]. As an alternative, a trial and error procedure has been frequently proposed [Deville et al 1993; Sautory 1993]: (i) start with a very large interval \code{bounds.0}; (ii) if convergence is achieved, shrink it so as to obtain a new interval \code{bounds.1}; (iii) repeat until you get a sufficiently tight feasible interval \code{bounds.n}. The drawback is that this procedure can cost a lot of computer time since, for each choice of the \code{bounds}, the full calibration problem has to be solved.

A rather easy task is, on the contrary, the one of finding at least a given specific interval \eqn{I^*=[L^*,U^*]} such that, if it is \emph{not} covered by \code{bounds}, the current calibration problem is \emph{surely unfeasible}. This means that any feasible \code{bounds} value must necessarily contain the \eqn{I^*} interval. The function \code{bounds.hint}: (i) first identifies such an \eqn{I^*} interval (by computing the range of the ratios between known population totals and corresponding direct Horvitz-Thompson estimates), (ii) then builds a new interval \eqn{I^{sugg}} with same midpoint and double length. The latter is the \emph{suggested} value for the \code{bounds} argument of \code{e.calibrate}. The return value of \code{bounds.hint} should be understood as a useful starting guess for \code{bounds}, even though there is definitely no warranty that the calibration algorithm will actually converge.
}

\author{Diego Zardetto}

\references{
Vanderhoeft, C. (2001) \emph{"Generalized Calibration at Statistic Belgium"}, Statistics Belgium Working Paper n. 3.

Deville, J.C., Sarndal, C.E. and Sautory, O. (1993) \emph{"Generalized Raking Procedures in Survey Sampling"}, Journal of the American Statistical Association, Vol. 88, No. 423, pp.1013-1020.

Sautory, O. (1993) \emph{"La macro CALMAR: Redressement d'un Echantillon par Calage sur Marges"}, Document de travail de la Direction des Statistiques Demographiques et Sociales, no. F9310.
}

\seealso{
\code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for constructing known totals data frames in compliance with the standard required by \code{e.calibrate}, \code{\link{population.check}} to check that the known totals data frame satisfies that standard, \code{\link{g.range}} to compute the range of the obtained g-weights, and \code{\link{check.cal}} to check if calibration constraints have been fulfilled. 
}

\examples{
# Creation of the object to be calibrated:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Calibration (partitioned solution) on the marginal distribution
# of age in 5 classes (age5c) inside provinces (procod)
# (totals in pop06p). Get a hint for feasible bounds:
hint<-bounds.hint(des,pop06p,~age5c-1,~procod)

# Let's verify if calibration converges with the suggested 
# value for the bounds argument (i.e. c(0.219, 1.786) ):
descal06p<-e.calibrate(design=des,df.population=pop06p,
           calmodel=~age5c-1,partition=~procod,calfun="logit",
           bounds=hint,aggregate.stage=2)

# Now let's verify that calibration fails, if bounds don't cover
# the interval [0.611, 1.394]:
\dontrun{
descal06p<-e.calibrate(design=des,df.population=pop06p,
           calmodel=~age5c-1,partition=~procod,calfun="logit",
           bounds=c(0.62,1.50),aggregate.stage=2,force=FALSE)
}
# The warning message raised by e.calibrate tells that
# the population total of variable age5c5 (i.e. the fifth
# age class frequency) was not matched.

# By analysing ecal.status one understands that calibration
# failed due to the sub-task identified by procod 30:
ecal.status

# this is easily explained by inspecting the "bounds"
# attribute of the bounds.hint output object:
attr(hint,"bounds")

# indeed the specified lower bound (0.62) was too high
# for procod 30, where instead a value ~0.61 was required.

# Recall that you can always "force" a calibration task that
# would not converge:
descal06p.forced<-e.calibrate(design=des,df.population=pop06p,
                  calmodel=~age5c-1,partition=~procod,calfun="logit",
                  bounds=c(0.62,1.50),aggregate.stage=2,force=TRUE)

# Notice, also, that forced sub-tasks can be tracked down by
# directly looking at ecal.status...
ecal.status

# ...or by using function check.cal:
check.cal(descal06p.forced)
}
\keyword{survey}