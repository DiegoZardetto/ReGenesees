\name{getBest}
\alias{getBest}
\alias{getBest.default}
\alias{getBest.gvf.fit}
\alias{getBest.gvf.fits}
\alias{getBest.gvf.fit.gr}
\alias{getBest.gvf.fits.gr}

\title{Identify the Best Fit GVF Model}

\description{
Given a set of competing fitted GVF models, this function selects the best model according to a given criterion. 
}

\usage{
getBest(object,
        criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
}

\arguments{
  \item{object}{Typically, an object containing \emph{many} fitted GVF models (i.e. of class \code{gvf.fits} or \code{gvf.fits.gr}).}
  \item{criterion}{The quality criterion to be used for model selection. Default is \code{R^2}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{
Given a set of competing fitted GVF models, this function selects the best model according to a given criterion.

Four goodness-of-fit criteria are available: \code{R^2}, \code{adjusted R^2}, \code{AIC}, and \code{BIC} (see \code{\link{getR2}}).

If \code{object} is a \emph{set} of GVF models fitted to \emph{grouped} data (i.e. of class \code{gvf.fits.gr}), the function will return the fitted GVF model with best \emph{average score} in the given criterion \emph{over the groups}.
}

\value{
A single GVF fitted model.
}

\section{Methodological Warning}{
Each one of the available criteria has its own specificities and limitations (e.g. it is senseless to use AIC to compare two GVF models with different response variables). It is up to the user to select the measure which is appropriate to his goals.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load example data:
data(AF.gvf)

# Inspect available estimates and errors of counts:
str(ee.AF)

# List available registered GVF models:
GVF.db


## (A) A *a set* of GVF models fitted to the same data ##
# Fit example data to all registered GVF models:
mm <- fit.gvf(ee.AF)
summary(mm)

# Get the best model according to adjusted R^2:
mm.best <- getBest(mm, criterion = "adj.R2")
mm.best

# NOTE: The *first* model has been selected. A thorough model comparison
#       by means of diagnostic plots would have led to the same result:
plot(mm, 1:3)


## (B) a *set of* GVF models fitted to *grouped* data ##
# We have at our disposal the following survey design object on household data:
exdes

# Use function svystat to prepare *grouped* estimates and errors of counts
# to be fitted separately (here groups are regions):
ee.g <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
str(ee.g)

# Fit all registered GVF model number separately inside groups:
mm.g <- fit.gvf(ee.g)
summary(mm.g)

# Get the best model according to R^2:
mm.g.best <- getBest(mm.g)
mm.g.best

# NOTE: Again, the *first* model has been selected. A thorough model comparison
#       by means of diagnostic plots would have led to the same result:
plot(mm.g, 1:3)
}

\keyword{survey}