\name{plot.gvf.fit}
\alias{plot.gvf.fit}
\alias{plot.gvf.fits}
\alias{plot.gvf.fit.gr}
\alias{plot.gvf.fits.gr}

\title{
Diagnostic Plots for Fitted GVF Models
}

\description{
This function provides basic diagnostic plots for fitted GVF model(s).
}

\usage{
\method{plot}{gvf.fit}(x, which.more = 1:3, id.n = 3, labels.id = names(residuals(x)),
     cex.id = 0.75, label.pos = c(4, 2), cex.caption = 1, Main = NULL, ...)

\method{plot}{gvf.fits}(x, which.more = NULL, id.n = 3, labels.id = names(residuals(x)),
     cex.id = 0.75, label.pos = c(4, 2), cex.caption = 1, Main = NULL, ...)

\method{plot}{gvf.fit.gr}(x, which.more = 1:3, id.n = 3, labels.id = NULL,
     cex.id = 0.75, label.pos = c(4, 2), cex.caption = 1, ...)

\method{plot}{gvf.fits.gr}(x, which.more = NULL, id.n = 3, labels.id = NULL,
     cex.id = 0.75, label.pos = c(4, 2), cex.caption = 1, ...)
}

\arguments{
  \item{x}{An object containing one or more fitted GVF models (see \sQuote{Usage} for the allowed classes).}
  \item{which.more}{Select additional plots beyond the default one (\sQuote{Observed vs Fitted}). 
                    Can be any subset of vector \code{1:6} with up to three elements.}
  \item{id.n}{Number of points to be initially labelled in each plot, starting with the most extreme.}
  \item{labels.id}{Vector of labels, from which the labels for extreme points will be chosen. \code{NULL} uses observation numbers.}
  \item{cex.id}{Magnification of point labels.}
  \item{label.pos}{Positioning of labels, for the left half and right half of the graph(s) respectively.}
  \item{cex.caption}{Controls the size of \code{caption}.}
  \item{Main}{Optional string to be added to automatic plot titles.}
  \item{\dots}{Other parameters to be passed through to plotting functions.}
}

\details{
Diagnostic plots can be useful both for assessing the goodness of a GVF model fit qualitatively, and for selecting the \dQuote{best} GVF model among different alternatives.

This function can provide any of the following \code{7} plots:
\preformatted{
 (0)  'Observed vs Fitted'
 (1)  'Residuals vs Fitted'
 (2)  'Normal Q-Q'
 (3)  'Scale-Location'
 (4)  'Cook's distance'
 (5)  'Residuals vs Leverage'
 (6)  'Cook's distances vs Leverage/(1-Leverage)'
 }

The \emph{\sQuote{Residuals vs Fitted}} plot is special in that it will be always provided: this explains its zero-th order in the list above. The rest of the list, namely plots \code{1:6}, exactly matches the numbering convention of function \code{\link{plot.lm}}.

Additional plots - beyond \emph{\sQuote{Residuals vs Fitted}} - can be requested through argument \code{which.more}. Any subset of \code{1:6} is allowed, provided its \code{length} does not exceed \code{3}. Therefore, at most \code{4} plots will be generated simultaneously.

Note that the \emph{default} behaviour of this function do depend on whether input object \code{x} stores \emph{one} or \emph{more than one} fitted GVF models. In the first case, plots \code{0:3} will be returned in a multiple plot with a 2x2 layout. In the second case, only the default plot number \code{0} will be returned, opening a new graphics frame for each different GVF model.

Argument \code{id.n} specifies how many points have to be labelled, starting with the most extreme in terms of residuals: this applies to all plots.

Argument \code{Main} is expected to be seldom (if ever) useful: its main purpose is programming consistency at a deeper level.

All the other arguments have the same meaning as in function \code{\link{plot.lm}}.
}

\references{
See \code{\link{plot.lm}} and references therein.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load example data:
data(AF.gvf)

# Inspect available estimates and errors of counts:
str(ee.AF)

# List available registered GVF models:
GVF.db


############################################
# Diagnostic plots for fitted GVF model(s) #
############################################

## (A) Plots of a *single* fitted GVF model ##
# Fit example data to registered GVF model number one:
m <- fit.gvf(ee.AF, 1)

## Default call yields 4 plots:
plot(m)

# Play with argument 'which.more':
## which.more = NULL yields the "Observed vs Fitted" plot only:
plot(m, which.more = NULL)

## which.more = 1 adds the "Residuals vs Fitted" plot:
plot(m, which.more = 1)

## subsets of 1:6 with length <= 3 are allowed:
plot(m, which.more = c(1:2,4))

# Just for illustration, play with other parameters:
plot(m, id.n = 6, col = "blue", pch = 20)


## (B) Plots of *many* fitted GVF models ##
# Fit example data to registered GVF models for frequencies (i.e. number 1:3):
mm <- fit.gvf(ee.AF, 1:3)

## Default call yields the "Observed vs Fitted" plot reported separately for
## each model in subsequent graphics frames:
plot(mm)

# Play with argument 'which.more':
## which.more = 1:3 yields subsequent 2x2 plots:
plot(mm, which.more = 1:3)

## again, subsets of 1:6 are allowed:
plot(mm, which.more = 1)


######################################################
# Diagnostic plots for "grouped" fitted GVF model(s) #
######################################################
# We have at our disposal the following survey design object on household data:
exdes

# Use function svystat to prepare "grouped" estimates and errors of counts
# to be fitted separately (here groups are regions):
ee <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
lapply(ee, head)

## (C) Plots of a *single* GVF model fitted to different groups ##
# Fit registered GVF model number one separately inside groups:
m.g <- fit.gvf(ee, 1)

## Default call yields 4 plots reported separately for each group
## in subsequent graphics frames:
plot(m.g)

## Play with argument 'which.more' to select different plots:
plot(m.g, which.more = c(1:2,4))


## (D) Plots of *many* GVF model fitted to different groups ##
# Fit all registered GVF models for frequencies separately inside groups:
mm.g <- fit.gvf(ee, 1:3)

## Default call yields the "Residuals vs Fitted" plot reported separately
## for each group in subsequent graphics frames:
plot(mm.g)

## Play with argument 'which.more' to add more plots:
plot(mm.g, which.more = 1:3, id.n = 6, col = "blue", pch = 20)
}

\keyword{survey}
