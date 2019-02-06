\name{gvf.misc}
\alias{gvf.misc}
\alias{coef.gvf.fit}
\alias{coef.gvf.fits}
\alias{coef.gvf.fit.gr}
\alias{coef.gvf.fits.gr}
\alias{residuals.gvf.fit}
\alias{residuals.gvf.fits}
\alias{residuals.gvf.fit.gr}
\alias{residuals.gvf.fits.gr}
\alias{fitted.gvf.fit}
\alias{fitted.gvf.fits}
\alias{fitted.gvf.fit.gr}
\alias{fitted.gvf.fits.gr}
\alias{predict.gvf.fit}
\alias{predict.gvf.fits}
\alias{predict.gvf.fit.gr}
\alias{predict.gvf.fits.gr}
\alias{effects.gvf.fit}
\alias{effects.gvf.fits}
\alias{effects.gvf.fit.gr}
\alias{effects.gvf.fits.gr}
\alias{rstandard.gvf.fit}
\alias{rstandard.gvf.fits}
\alias{rstandard.gvf.fit.gr}
\alias{rstandard.gvf.fits.gr}
\alias{rstudent.gvf.fit}
\alias{rstudent.gvf.fits}
\alias{rstudent.gvf.fit.gr}
\alias{rstudent.gvf.fits.gr}
\alias{anova.gvf.fit}
\alias{anova.gvf.fits}
\alias{anova.gvf.fit.gr}
\alias{anova.gvf.fits.gr}
\alias{vcov.gvf.fit}
\alias{vcov.gvf.fits}
\alias{vcov.gvf.fit.gr}
\alias{vcov.gvf.fits.gr}

\title{Miscellanea: Methods for Fitted GVF Models}

\description{These methods extract information from fitted GVF model(s).}

\usage{
\method{coef}{gvf.fit}(object, ...)
\method{coef}{gvf.fits}(object, ...)
\method{coef}{gvf.fit.gr}(object, ...)
\method{coef}{gvf.fits.gr}(object, ...)

\method{residuals}{gvf.fit}(object, ...)
\method{residuals}{gvf.fits}(object, ...)
\method{residuals}{gvf.fit.gr}(object, ...)
\method{residuals}{gvf.fits.gr}(object, ...)

\method{fitted}{gvf.fit}(object, ...)
\method{fitted}{gvf.fits}(object, ...)
\method{fitted}{gvf.fit.gr}(object, ...)
\method{fitted}{gvf.fits.gr}(object, ...)

\method{predict}{gvf.fit}(object, ...)
\method{predict}{gvf.fits}(object, ...)
\method{predict}{gvf.fit.gr}(object, ...)
\method{predict}{gvf.fits.gr}(object, ...)

\method{effects}{gvf.fit}(object, ...)
\method{effects}{gvf.fits}(object, ...)
\method{effects}{gvf.fit.gr}(object, ...)
\method{effects}{gvf.fits.gr}(object, ...)

\method{rstandard}{gvf.fit}(model, ...)
\method{rstandard}{gvf.fits}(model, ...)
\method{rstandard}{gvf.fit.gr}(model, ...)
\method{rstandard}{gvf.fits.gr}(model, ...)

\method{rstudent}{gvf.fit}(model, ...)
\method{rstudent}{gvf.fits}(model, ...)
\method{rstudent}{gvf.fit.gr}(model, ...)
\method{rstudent}{gvf.fits.gr}(model, ...)

\method{anova}{gvf.fit}(object, ...)
\method{anova}{gvf.fits}(object, ...)
\method{anova}{gvf.fit.gr}(object, ...)
\method{anova}{gvf.fits.gr}(object, ...)

\method{vcov}{gvf.fit}(object, ...)
\method{vcov}{gvf.fits}(object, ...)
\method{vcov}{gvf.fit.gr}(object, ...)
\method{vcov}{gvf.fits.gr}(object, ...)
}

\arguments{
  \item{object}{An object containing one or more fitted GVF models (see \sQuote{Usage} for the allowed classes).}
  \item{model}{An object containing one or more fitted GVF models (see \sQuote{Usage} for the allowed classes).}
  \item{\dots}{Further arguments passed to or from other methods (see corresponding .lm methods).}
}

\details{
These methods can be used to extract information from fitted GVF model(s).

For more details on their usage, please read the help pages of the methods with same name defined on class .lm by package \pkg{stats} (e.g. \code{\link[stats]{coef}}, \code{\link[stats]{fitted}}, etc.).
}

\value{
The requested information, wrapped into an \R object whose structure depends on the class of the input fitted GVF model(s) (i.e. arguments \code{object} and/or \code{model}).
}

\author{
Diego Zardetto
}

\seealso{\code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load example data:
data(AF.gvf)

# Inspect available estimates and errors of counts:
head(ee.AF)
summary(ee.AF)

# List available registered GVF models:
GVF.db

## (A) A *single* fitted GVF model ##
# Fit example data to registered GVF model number one:
m <- fit.gvf(ee.AF, 1)

# Extract some information:
coef(m)
fitted(m)


## (B) A *a set* of GVF models fitted to the same data ##
# Fit example data to registered GVF models for frequencies (i.e. number 1:3):
mm <- fit.gvf(ee.AF, 1:3)

# Extract some information:
r.mod <- residuals(mm)
lapply(r.mod, head)

r.sta <- rstandard(mm)
lapply(r.sta, head)

r.stu <- rstudent(mm)
lapply(r.stu, head)


## (C) a *single* GVF model fitted to *grouped* data ##
# We have at our disposal the following survey design object on household data:
exdes

# Use function svystat to prepare *grouped* estimates and errors of counts
# to be fitted separately (here groups are regions):
ee.g <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
str(ee.g)

# Fit registered GVF model number one separately inside groups:
m.g <- fit.gvf(ee.g, 1)

# Extract some information:
coef(m.g)
fitted(m.g)


## (D) a *set of* GVF models fitted to *grouped* data ##
# Fit all registered GVF models for frequencies (i.e. number 1:3) separately
# inside groups:
mm.g <- fit.gvf(ee.g, 1:3)
# Extract some information:
coef(mm.g)
fitted(mm.g)
}

\keyword{survey}