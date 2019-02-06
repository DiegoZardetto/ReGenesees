\name{ext.calibrated}
\alias{ext.calibrated}

\title{Make ReGenesees Digest Externally Calibrated Weights}

\description{
Enables \pkg{ReGenesees} to provide correct variance estimates of (functions of) calibration estimators, even if the survey weights have not been calibrated by \pkg{ReGenesees}.
}

\usage{
ext.calibrated(data, ids, strata = NULL, weights,
               fpc = NULL, self.rep.str = NULL, check.data = TRUE,
               weights.cal, calmodel, partition = FALSE, sigma2 = NULL)
}

\arguments{
  \item{data}{The same as in function \code{\link{e.svydesign}}.}
  \item{ids}{The same as in function \code{\link{e.svydesign}}.}
  \item{strata}{The same as in function \code{\link{e.svydesign}}.}  
  \item{weights}{The same as in function \code{\link{e.svydesign}}.}
  \item{fpc}{The same as in function \code{\link{e.svydesign}}.} 
  \item{self.rep.str}{The same as in function \code{\link{e.svydesign}}.}
  \item{check.data}{The same as in function \code{\link{e.svydesign}}.}
  \item{weights.cal}{Formula identifying the externally calibrated weights.}
  \item{calmodel}{The same as in function \code{\link{e.calibrate}}.}
  \item{partition}{The same as in function \code{\link{e.calibrate}}.}
  \item{sigma2}{The same as in function \code{\link{e.calibrate}}.}
}

\details{
Owing to \pkg{ReGenesees}'s ability to provide proper variance estimates for (complex functions of) calibration estimators, some users may be tempted to exploit \pkg{ReGenesees} in the estimation phase \emph{even if} they did \emph{not} use \pkg{ReGenesees} for calibration.

This result \emph{cannot} be achieved naively, by simply passing to \pkg{ReGenesees} function \code{\link{e.svydesign}} the survey data and supplying the externally calibrated weights through its \code{weights} argument.

Indeed, variance estimation methods of \pkg{ReGenesees}'s summary statistics functions (\code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}} and \code{\link{svystat}}) are dispatched according to the class of the input \code{design} object:
\enumerate{
\item If the design object is un-calibrated (i.e. its class is \sQuote{analytic}), variance formulas are appropriate to Horvitz-Thompson estimators (and functions of them).

\item If the design object is calibrated (i.e. its class is \sQuote{cal.analytic}), variance formulas are appropriate to Calibration estimators (and functions of them).
}

Therefore, the naive approach of passing the externally calibrated weights \code{weights.cal} to \code{\link{e.svydesign}} as if they were initial or design weights cannot succeed, since it would result in HT-like variance estimates, leading generally to \emph{variance overestimation} (with bigger upward bias for variables that are better explained by the calibration model).

Function \code{ext.calibrated} has been designed exactly to avoid the aforementioned pitfalls and to allow \pkg{ReGenesees} provide correct variance estimates of (functions of) calibration estimators, even if the survey weights have been calibrated externally by other software.

Argument \code{weights.cal} identifies the externally calibrated weights of the units included in the sample. The \code{data} variable referenced by \code{weights.cal} must be numeric. Currently, only \emph{positive} externally calibrated weights can be handled (see the dedicated section below).

Other arguments to \code{ext.calibrated} derive either from function \code{\link{e.svydesign}} or from function \code{\link{e.calibrate}}. The former serve the purpose of passing the survey data and the corresponding sampling design metadata, the latter are meant to tell \code{ext.calibrated} how the externally calibrated weights have been obtained.
}

\section{What if externally calibrated weights happen to be negative?}{
From a methodological perspective, negative calibration weights are legitimate. However, owing to software implementation details whose modification would not be trivial, function \code{ext.calibrated} is \emph{not} yet able to cope with this case. Note that the problem is actually due to the \emph{external} origin of the negative calibration weights. In fact, \pkg{ReGenesees} calibration and estimation facilities are entirely able to cope with possibly negative calibration weights, provided they were computed \emph{internally}.
}

\value{
An object of class \code{cal.analytic}, storing the original survey data \emph{plus} all the sampling design and calibration metadata needed for proper variance estimation.
}

\note{
Exactly as \pkg{ReGenesees}'s base functions \code{\link{e.svydesign}} and \code{\link{e.calibrate}} would do, \code{ext.calibrated} too will wrap inside its return value a local copy of \code{data}. As usual, this copy will be stored inside the \code{variables} slot of the output list. As usual, again, the calibrated weights will be accessible by using the \code{\link{weights}} function.
}

\author{Diego Zardetto.}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata, and \code{\link{e.calibrate}} for calibrating survey weights within \pkg{ReGenesees}.
}

\examples{

# Load data sbs data
data(sbs)

#########################################################################
# Simulate an external calibration procedure and compute some benchmark #
# estimates and errors to test function ext.calibrated                  #
#########################################################################
# Define a survey design
sbsdes <- e.svydesign(data= sbs, ids= ~id, strata= ~strata, weights= ~weight,
                      fpc= ~fpc)

# Build a template for population totals
pop <- pop.template(data= sbsdes, calmodel= ~y:nace.macro + emp.cl + emp.num - 1,
                    partition= ~dom3)

# Have a look at the template structure
pop.desc(pop)

# Fill the template
pop <- fill.template(universe= sbs.frame, template= pop)

# Calibrate
sbscal <- e.calibrate(design= sbsdes, df.population= pop, calfun= "logit",
                      bounds= c(0.8, 1.3), sigma2= ~ emp.num)

# Compute benchmark estimates and errors (average value added per employee by
# region) to be later compared with those obtained by using ext.calibrated 
benchmark <- svystatR(design= sbscal, num= ~va.imp2, den= ~emp.num, by= ~region)
benchmark

# Extract the 'externally' calibrated weights...
w <- weights(sbscal)

#...and add these 'externally' calibrated weights to the original survey data
sbs.ext <- data.frame(sbs, w.ext = w)

# NOTE: Now sbs.ext is just a data frame, without any knowledge of the
#       calibration metadata formerly stored inside sbscal (i.e. the object
#       calibrated by ReGenesees)


##############################################################
# Let ReGenesees digest the 'externally' calibrated weights, #
# then re-compute benchmark estimates and errors for testing #
##############################################################
# Simply pass survey data along with sampling design and calibration model
# metadata
sbscal.ext <- ext.calibrated(data= sbs.ext, ids= ~id, strata= ~strata,
                             weights= ~weight, fpc = ~fpc,
                             weights.cal= ~w.ext,
                             calmodel= ~y:nace.macro + emp.cl + emp.num - 1,
                             partition= ~dom3, sigma2= ~emp.num)

# Have a look at the output
sbscal.ext

# Now re-compute benchmark estimates and errors by means of new object
# ext.sbscal
test <- svystatR(design= sbscal.ext, num= ~va.imp2, den= ~emp.num, by= ~region)
test

################################################################
# Compare benchmark estimates and errors to those derived from #
# ext.calibrated return object                                 #
################################################################
benchmark
test

# ...and they are identical, as it must be.

# NOTE: All utility tools yield exactly the same results, e.g.
identical(weights(sbscal), weights(sbscal.ext))
identical(g.range(sbscal), g.range(sbscal.ext))


##########################################################################
# Show that the naive idea of directly passing the externally calibrated #
# weights to e.svydesign does NOT work properly for variance estimation  #
##########################################################################
naive <- e.svydesign(data= sbs.ext, ids= ~id, strata= ~strata,
                     weights= ~w.ext, fpc = ~fpc)

# Estimated sampling errors derived by this naive design object...
svystatR(design= naive, num= ~va.imp2, den= ~emp.num, by= ~region)

#...do NOT match benchmark values, overestimating them:
benchmark
}
\keyword{survey}