\name{getR2}
\alias{getR2}
\alias{getR2.default}
\alias{getR2.gvf.fit}
\alias{getR2.gvf.fits}
\alias{getR2.gvf.fit.gr}
\alias{getR2.gvf.fits.gr}
\alias{AIC}
\alias{AIC.gvf.fits}
\alias{AIC.gvf.fit.gr}
\alias{AIC.gvf.fits.gr}
\alias{BIC}
\alias{BIC.gvf.fits}
\alias{BIC.gvf.fit.gr}
\alias{BIC.gvf.fits.gr}

\title{Quality Measures on Fitted GVF Models}

\description{
These functions extract goodness-of-fit measures from fitted GVF models.
}

\usage{
getR2(object, adjusted = FALSE, ...)

\method{AIC}{gvf.fits}(object, ...)

\method{BIC}{gvf.fits}(object, ...)
}

\arguments{
  \item{object}{An object containing one or more fitted GVF models.}
  \item{adjusted}{Should the \emph{adjusted} \code{R^2} be computed? The default is \code{FALSE}}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{
These functions compute three goodness-of-fit measures on fitted GVF models: \code{R^2}, \link[stats]{AIC}, and \link[stats]{BIC}. Such measures can help compare the relative quality of competing GVF models, hence facilitating model selection (see also function \code{\link{getBest}}).

Though \code{object} can also be a \emph{single} fitted GVF model, these functions are principally meant to compare \emph{different} GVF models fitted to the same data (i.e. the same estimates and errors).

To request the \emph{adjusted} \code{R^2}, use function \code{getR2} and specify \code{adjusted = TRUE}.
}

\value{
If \code{object} is a single GVF model (class \code{gvf.fit}), the requested quality measure.

If \code{object} is a set of GVF models fitted to the same data (class \code{gvf.fits}), a vector whose elements store the requested quality measure for each GVF model.

If \code{object} is a single GVF model fitted to "grouped" data (class \code{gvf.fit.gr}), a list whose components store the requested quality measure for the corresponding groups.

If \code{object} is a set of GVF models fitted to "grouped" data (class \code{gvf.fits.gr}), a list whose components store vectors whose elements report the requested quality measure for each GVF model of each group.
}

\section{Methodological Warning}{
Each one of the provided quality measures has its own specificities and limitations (e.g. it is senseless to use AIC to compare two GVF models with different response variables). It is up to the user to select the measure which is appropriate to his goals.
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

## (A) A *single* fitted GVF model ##
# Fit example data to registered GVF model number one:
m <- fit.gvf(ee.AF, 1)

# Compute some goodness-of-fit measures:
getR2(m)
AIC(m)


## (B) A *a set* of GVF models fitted to the same data ##
# Fit example data to all registered GVF models:
mm <- fit.gvf(ee.AF)

# Compute some goodness-of-fit measures:
getR2(mm, adjusted = TRUE)
BIC(mm)


## (C) a *single* GVF model fitted to *grouped* data ##
# We have at our disposal the following survey design object on household data:
exdes

# Use function svystat to prepare *grouped* estimates and errors of counts
# to be fitted separately (here groups are regions):
ee.g <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
str(ee.g)

# Fit registered GVF model number one separately inside groups:
m.g <- fit.gvf(ee.g, 1)

# Compute some goodness-of-fit measures:
getR2(mm)
AIC(mm)


## (D) a *set of* GVF models fitted to *grouped* data ##
# Fit all registered GVF model number separately inside groups:
mm.g <- fit.gvf(ee.g)

# Compute some goodness-of-fit measures:
getR2(mm.g, adjusted = TRUE)
BIC(mm.g)
}

\keyword{survey}