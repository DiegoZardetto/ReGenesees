\name{fit.gvf}
\alias{fit.gvf}
\alias{fit.gvf.default}
\alias{fit.gvf.gvf.input}
\alias{fit.gvf.gvf.input.gr}
\alias{gvf.fit}
\alias{gvf.fits}
\alias{gvf.fit.gr}
\alias{gvf.fits.gr}
\alias{print.gvf.fit}
\alias{print.gvf.fits}
\alias{print.gvf.fit.gr}
\alias{print.gvf.fits.gr}
\alias{[.gvf.fits}
\alias{[[.gvf.fits}
\alias{summary.gvf.fit}
\alias{summary.gvf.fits}
\alias{summary.gvf.fit.gr}
\alias{summary.gvf.fits.gr}
\alias{GVF}

\title{
Fit GVF Models
}

\description{
This function fits one or more GVF models to a set of survey statistics.
}

\usage{
fit.gvf(gvf.input, model = NULL, weights = NULL)

\method{print}{gvf.fit}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{gvf.fits}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{gvf.fit.gr}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{gvf.fits.gr}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{[}{gvf.fits}(x, ...)
\method{[[}{gvf.fits}(x, ...)

\method{summary}{gvf.fit}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
\method{summary}{gvf.fits}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
\method{summary}{gvf.fit.gr}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
\method{summary}{gvf.fits.gr}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}

\arguments{
  \item{gvf.input}{An object of class \code{gvf.input} (or \code{gvf.input.gr}), containing the data to fit.}
  \item{model}{The GVF model(s) to be fitted (see \sQuote{Details}). \code{NULL} (the default) requires to fit \emph{all} the registered GVF models currently available in \code{\link{GVF.db}}.}
  \item{weights}{Formula specifying the weights to be used for fitting (via weighted least squares), if any. \code{NULL} (the default) means that ordinary least squares will be used. See also \sQuote{Details}.}

  \item{x}{An object of class \code{gvf.fits}, storing fitted GVF models.}
  \item{digits}{Minimal number of significant digits, see \code{\link{print.default}}.}

  \item{object}{Any output of \code{fit.gvf}, storing one (more than one) fitted GVF model(s).}
  \item{correlation}{Should the correlation matrix of the estimated parameters be returned and printed? Logical, with default \code{FALSE}.}
  \item{symbolic.cor}{Should the correlations be printed in symbolic form (see \code{\link{symnum}}) rather than as numbers. Logical, with default \code{FALSE}.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{
Function \code{fit.gvf} fits one or more GVF models to a set of survey statistics. The rationale for fitting multiple models to the same data is primarily for comparison purposes: the user is expected to eventually choose his preferred model, in order to obtain sampling errors predictions.

Argument \code{gvf.input} specifies the set of (pre computed) estimates and errors to which GVF models are to be fitted, as prepared by functions \code{\link{gvf.input}} and/or \code{\link{svystat}}.

One or more GVF models can be fitted \emph{simultaneously} to the same data, depending on the way argument \code{model} is passed.

Argument \code{model} can be either:

(1) \code{NULL} (the default) meaning \emph{all} the registered models currently available in \code{\link{GVF.db}};

(2) any sub-vector of \code{GVF.db$Model.id}, i.e. an integer vector identifying an arbitrary selection of registered models;

(3) an arbitrary (single) formula, i.e. any custom, user-defined GVF model.\cr

When \code{model} is passed via options (1) or (2), function \code{fit.gvf} can take advantage of any additional information available inside \code{GVF.db}, e.g. to warn the user in case a GVF model is not deemed to be appropriate for the kind of estimates contained into \code{gvf.input} (see \sQuote{Examples}).

Argument \code{weights} enables fitting the specified GVF model(s) via weighted least squares. By default \code{weights = NULL} and ordinary least squares are used. The weights must be passed by a formula referencing variables belonging to \code{gvf.input}. For instance, to weight observations according to reciprocals of squared CVs, one can use \code{weights} \code{=} \code{~I(CV^-2)}.
}

\value{
An object containing one or more fitted GVF models, depending on the way argument \code{model} was passed.

Let's first focus on input objects of class \code{gvf.input}.\cr
If \code{model} specifies a single GVF model, the output object will be of class \code{gvf.fit} and inherit from class \code{\link{lm}}.\cr
If \code{model} specifies many GVF models, the output object will be of class \code{gvf.fits} and inherit from class \code{\link{list}}. Hence, it will be possible to subset \code{gvf.fits} objects via methods \code{[} and \code{[[}. Note, moreover, that each component (in the sense of class \code{\link{list}}) of a \code{gvf.fits} object will be of class \code{gvf.fit}.

When, instead, the input object has class \code{\link{gvf.input.gr}}, i.e. it stores \emph{\dQuote{grouped}} estimates and errors, model fitting is performed \emph{separately} for different groups. Therefore, applying \code{fit.gvf} \emph{always} results in \emph{many} fitted GVF models.\cr
If \code{model} specifies a single GVF model, the output object will be of class \code{gvf.fit.gr} and inherit from class \code{\link{list}}. Each slot of the list will contain the same GVF model fitted to a specific group.\cr
If \code{model} specifies many GVF models, the output object will be of class \code{gvf.fits.gr} and again inherit from class \code{\link{list}}. Each slot of the list will now contain \emph{a second list} storing different GVF models fitted to a specific group.\cr
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{estimator.kind}} to assess what kind of estimates are stored inside a survey statistic object, \code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load example data:
data(AF.gvf)

# Now we have at our disposal a set of estimates and errors
# of Absolute Frequencies:
str(ee.AF)

# And the available registered GVF models are listed below:
GVF.db


###########################################
# How to specify the GVF model(s) to fit? #
###########################################

## (A) How to specify a *single* GVF model ##

#### (A.1) Select one registered model using its 'Model.id' as reported in
####       the GVF.db archive
# Let's fit, for instance, the GVF model with Model.id = 1:
m <- fit.gvf(ee.AF, model = 1)

# Inspect the result:
class(m) 
m
summary(m)

# Now let's fit GVF model with Model.id = 4
m <- fit.gvf(ee.AF, model = 4)
# Beware of the NOTE reported when printing or summarizing this fitted model:
m
summary(m)


#### (A.2) Specify the GVF model to fit by providing its formula directly, e.g.
####       because it is not available in GVF.db (yet):
m <- fit.gvf(ee.AF, model = CV ~ I(1/Y^2) + I(1/Y) + Y + I(Y^2))
m
summary(m)


## (B) How to specify a *many* GVF models simultaneously ##

#### (B.1) Use a subset of column 'Model.id' of GVF.db
# Let's, for instance, fit all the available GVF models which are appropriate
# to Frequencies, as reported in column 'Estimator.kind' of GVF.db
mm <- fit.gvf(ee.AF, model = 1:3)

# Inspect the result:
class(mm)
length(mm)
mm
summary(mm)

# Note that you can subset the output fitted models as a list:
mm.31 <- mm[c(3,1)]
class(mm.31)
mm.31

# and:
mm.2 <- mm[[2]]
class(mm.2)
mm.2


#### (B.2) Not specifying any GVF model, or specifying model = NULL, causes
####       *all* the available models in GVF.db to be fitted simultaneously:
mm <- fit.gvf(ee.AF)

# Inspect the result:
class(mm)
length(mm)
mm
summary(mm)


#########################################################
# How to fit GVF model(s) via *weighted* least squares? #
#########################################################
# Weights can be specified by a formula. Of course, the 'weights' formula must
# reference variables belonging to gvf.input.

# Let's use the built-in GVF model with Model.id = 1 and weight observations
# according to reciprocals of squared CVs:
mw <- fit.gvf(ee.AF, model = 1, weights = ~I(CV^-2))
mw

# Compute ordinary least squares fit:
m <- fit.gvf(ee.AF, model = 1)
m

# Compare the results:
summary(mw)
summary(m)


#########################################################################
# Fitting GVF model(s) to "grouped" estimates and errors: a quick ride. #
#########################################################################
# Recall we have at our disposal the following survey design object 
# defined on household data:
exdes

# Now use function svystat to prepare "grouped" estimates and errors
# of counts to be fitted separately (here groups are regions):
ee.g <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
class(ee.g)
ee.g

## Fit a *single* registered GVF model separately inside groups ##
m.g <- fit.gvf(ee.g, model = 1)

# Inspect the result:
class(m.g)
length(m.g)
m.g
summary(m.g)

# Can subset the result as a list, e.g. to get the fitted model of region '7':
m.g7 <- m.g[["7"]]
class(m.g7)
summary(m.g7)


## Fit *many* registered GVF models separately inside groups ##
mm.g <- fit.gvf(ee.g, model = 1:3)

# Inspect the result:
class(mm.g)
length(mm.g)
mm.g
summary(mm.g)

# Still can subset the result as a list, but now each component is a list
# itself. To get the fitted models of region '7', simply:
mm.g7 <- mm.g[["7"]]
class(mm.g7)
summary(mm.g7)

# And to isolate GVF fitted model number 2 for region '7', simply:
mm.g7.2 <- mm.g7[[2]]
class(mm.g7.2)
summary(mm.g7.2)

}

\keyword{survey}
