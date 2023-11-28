\name{gvf.input}
\alias{gvf.input}
\alias{plot.gvf.input}

\title{
Prepare Input Data to Fit GVF Models
}

\description{
Transforms a set of computed survey statistics into a suitable (\code{data.frame}-like) data structure, in order to fit a Generalized Variance Function model.
}

\usage{
gvf.input(design, ..., stats = list(...))

\method{plot}{gvf.input}(x, ...)
}

\arguments{
  \item{design}{The design object (of class \code{analytic} or inheriting from it) from which the input survey statics are supposed to have been derived.}
  \item{\dots}{For function \code{gvf.input}, objects containing survey statistics. For \code{plot}, further arguments passed to or from other methods.}
  \item{stats}{A \code{list} storing survey statistic objects (see \sQuote{Details}).}
  \item{x}{The object of class \code{gvf.input} to plot.}
}

\details{
Given a set of survey statistic objects (via arguments \sQuote{\code{\ldots}} or \code{stats}) and a design object (\code{design}) from which those statics are supposed to have been derived, function \code{gvf.input} builds a data structure that can be fed to \pkg{ReGenesees} GVF model fitting function \code{\link{fit.gvf}}.

Argument \sQuote{\code{\ldots}} can be bound to an arbitrary number of objects. These objects must be output of survey statistics functions, i.e. \code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}}, \code{\link{svySigma}}, and \code{\link{svySigma2}}.\cr
All input objects passed to \sQuote{\code{\ldots}} must derive from estimators of the \emph{same} kind (as returned by function \code{\link{estimator.kind}}). For the same reason, objects of \emph{mixed} kind (see \code{\link{estimator.kind}}) are not allowed. Since function \code{\link{svystatL}} can actually handle estimators of \emph{different} kinds, objects of kind \sQuote{\code{Complex Estimator}} are the only exception to the rule.

Argument \code{stats} can be used as an alternative to argument \sQuote{\code{\ldots}}: one has only to store the survey statistic objects into a \code{list} and bind such list to \code{stats}. Note that, if both are passed, argument \code{stats} will prevail on \sQuote{\code{\ldots}} (see \sQuote{Examples}).

Should \emph{any} input object be a survey statistic derived from a design object \emph{other than} \code{design}, the function would raise an error.

The \code{plot} method for \code{gvf.input} objects produces a matrix of scatterplots with polynomial smoothers.
}

\value{
An object of class \code{gvf.input}, inheriting from class \code{data.frame}: basically a data frame supplied with appropriate attributes.\cr
Each row of the data frame contains an \emph{estimate} along with its estimated sampling error, expressed in terms of \emph{standard error}, \emph{coefficient of variation}, \emph{variance}, and - whenever available - \emph{design effect}.\cr
The data frame has the following structure:
  \describe{
    \item{\code{name}}{The name of the original estimate, \code{factor}.}
    \item{\code{Y}}{The value of the original estimate, \code{numeric}.}
    \item{\code{SE}}{The standard error of the original estimate, \code{numeric}.}
    \item{\code{CV}}{The coefficient of variation of the original estimate, \code{numeric}.}
    \item{\code{VAR}}{The variance of the original estimate, \code{numeric}.}
    \item{\code{DEFF}}{The design effect of the original estimate (if available), \code{numeric}.}
  }

Note that by inspecting the \emph{attributes} of a \code{gvf.input} object, one can always identify which design object and which kind of estimator generated that object (see \sQuote{Examples}). 
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{estimator.kind}} to assess what kind of estimates are stored inside a survey statistic object, \code{\link{svystat}} as a useful alternative to prepare the input for GVF model fitting, \code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load sbs data:
data(sbs)

# Create a design object...
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

# ...and use it to compute some survey statistics:
va<-svystatTM(sbsdes,~va.imp2)
va.reg<-svystatTM(sbsdes,~va.imp2,~region)
va.area<-svystatTM(sbsdes,~va.imp2,~area)

# Now suppose you want to fit a GVF model on the estimates and errors computed
# above: you must prepare your input as follows:
ee<-gvf.input(sbsdes,va,va.reg,va.area)

# Inspect the obtained data structure:
ee
str(ee)   # Note the "design" and "stats.kind" attributes
plot(ee)

# Note that, instead of argument '...', you could have used argument 'stats'
# as follows:
va.list<-list(va,va.reg,va.area)
ee2<-gvf.input(sbsdes,stats=va.list)

# ...obtaining exactly the same result:
identical(ee,ee2)

# Note also that, if both are passed, argument 'stats' prevails on '...':
  # indeed, while:
  gvf.input(sbsdes,va.reg)
  # we would get again:
  gvf.input(sbsdes,va.reg,stats=va.list)

}
\keyword{survey}
