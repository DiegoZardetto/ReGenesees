\name{aux.estimates}
\alias{aux.estimates}

\title{Quick Estimates of Auxiliary Variables Totals}

\description{
Quickly estimates the totals of the auxiliary variables of a calibration model.
}

\usage{
aux.estimates(design,
              calmodel = if (inherits(template, "pop.totals"))
                         attr(template, "calmodel"), 
              partition = if (inherits(template, "pop.totals")) 
                          attr(template, "partition") else FALSE,
              template = NULL)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model (see \sQuote{Details}); \code{FALSE} (the default) implies no calibration domains.}
  \item{template}{An object of class \code{pop.totals}, be it a template or the actual known totals data frame for the calibration task.}
}

\details{
The main purpose of function \code{aux.estimates} is to make easy the task of estimating the totals of \emph{all} the auxiliary variables involved in a calibration model (separately inside distinct calibration domains, if specified). Even if such totals can be estimated also by repeatedly invoking function \code{\link{svystatTM}}, this may reveal very tricky in practice, because real-world calibration tasks (e.g. in the field of Official Statistics) can simultaneously involve hundreds of auxiliary variables. Moreover, total estimates provided by function \code{\link{svystatTM}} are always complemented by sampling errors, whose estimation is very computationally demanding.

Function \code{aux.estimates}, on the contrary, \emph{only} provides estimates of totals (i.e. without associated sampling errors), thus being very quick to be executed. Moreover, \code{aux.estimates} is able to compute, \emph{in just a single shot}, all the totals of the auxiliary variables of a calibration model, no matter how complex the model is. Lastly, as a third strong point, the totals estimated by \code{aux.estimates} will be returned exactly in the same \emph{standard format} in which the known population totals for the related calibration task need to be represented (see \code{\link{pop.template}}, \code{\link{population.check}}, \code{\link{fill.template}}).

It may be useful to point out that, besides having been designed to handle auxiliary variables involved in calibration models, function \code{aux.estimates} could be also used for computing \emph{general} estimates of totals inside subpopulations in a very effective way (see \sQuote{Examples}).
}

\value{
An object of class \code{pop.totals}, thus inheriting from class \code{data.frame} storing the estimated totals in a standard format.
}

\author{Diego Zardetto}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata, \code{\link{svystatTM}} for calculating estimates and standard errors of totals, \code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for constructing known totals data frames in compliance with the standard required by \code{e.calibrate}, \code{\link{population.check}} to check that the known totals data frame satisfies that standard, \code{\link{fill.template}} to automatically fill the template when a sampling frame is available.
}

\examples{
# Load sbs data:
data(sbs)

# Build a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

# Now suppose you have to perform a calibration process which
# exploits as auxiliary information:
#  i) the total number of employees (emp.num)
#     by class of number of employees (emp.cl) crossed with nace.macro; 
# ii) the total number of enterprises (ent)
#     by region crossed with nace.macro;

# Build a template for the known totals:
pop<-pop.template(sbsdes,
     calmodel=~emp.num:emp.cl + region -1,
     partition=~nace.macro)

# Use the fill.template function to automatically compute
# the totals from the universe (sbs.frame) and safely fill
# the template:
pop<-fill.template(sbs.frame,template=pop)
pop

# You can now use aux.estimates to verify how much difference
# exists between the target totals and the initial HT estimates:
aux.HT<-aux.estimates(sbsdes,template=pop)
aux.HT 

# If you calibrate, ...
sbscal<-e.calibrate(sbsdes,pop)

# ... you can verify that CAL estimates exactly match the known totals:
aux.CAL<-aux.estimates(sbscal,template=pop)
aux.CAL

# Recall that you can also use aux.estimates for computing
# general estimates of totals inside subpopulations (even
# not related to any calibration task).
# E.g. estimate the total of value added inside areas:
aux.estimates(sbsdes,~va.imp2-1,~area)

# ...and compare to svystatTM (notice also
# the increased execution time):
svystatTM(sbsdes,~va.imp2,~area) 
}
\keyword{survey}