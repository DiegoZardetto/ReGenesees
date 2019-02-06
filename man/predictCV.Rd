\name{predictCV}
\alias{predictCV}
\alias{predictCV.default}
\alias{predictCV.gvf.fit}
\alias{predictCV.gvf.fits}
\alias{predictCV.gvf.fit.gr}
\alias{predictCV.gvf.fits.gr}

\title{
Predict CV Values via Fitted GVF Models
}

\description{
This function predicts the CV values associated to given estimates, based on fitted GVF model(s).
}

\usage{
predictCV(object, new.Y = NULL, scale = NULL, df = Inf,
          interval = c("none", "confidence", "prediction"), level = 0.95,
          na.action = na.pass, pred.var = NULL, weights = 1)
}

\arguments{
  \item{object}{An object containing one or more fitted GVF models.}
  \item{new.Y}{A data frame storing new estimates whose CVs have to be predicted. If omitted or \code{NULL}, CVs arising from the fitted GVF model(s) will be returned.}
  \item{scale}{Scale parameter for standard error calculation. 
               See also \code{\link{predict.lm}}.}
  \item{df}{Degrees of freedom for scale. 
            See also \code{\link{predict.lm}}.}
  \item{interval}{Type of interval calculation. Can be abbreviated. 
                  See also \code{\link{predict.lm}}.}
  \item{level}{Confidence (or tolerance) level for intervals. 
               See also \code{\link{predict.lm}}.}
  \item{na.action}{Function determining what should be done with missing values in \code{new.Y}. The default is to predict \code{NA}. 
                   See also \code{\link{predict.lm}}.}
  \item{pred.var}{The variance(s) for future observations to be assumed for prediction intervals. 
                  See also \code{\link{predict.lm}}.}
  \item{weights}{Variance weights for prediction. This can be a numeric vector or a one-sided model formula. In the latter case, it is interpreted as an expression evaluated in \code{new.Y}. See also \code{\link{predict.lm}}.}
}

\details{
The main motivation for building and fitting a GVF model is to exploit the fitted model to \emph{predict} the sampling error associated to a given estimate, instead of having to \emph{compute} directly an estimate of such sampling error. Function \code{predictCV} is relevant to that scope.

Despite different GVF models can specify as response term (call it \code{'resp'}) different functions of variables \code{'SE'}, \code{'CV'}, and \code{'VAR'} (see e.g. [Wolter 07]), function \code{predictCV} adopts variable \code{'CV'} as a universal \emph{pivot}. This means that \code{predictCV} can handle \emph{only} fitted GVF models which are \emph{registered} (that is already stored inside \code{\link{GVF.db}}), and for which variable \code{Resp.to.CV} is \emph{not} \code{NA}. Indeed, it is variable \code{Resp.to.CV} of data frame \code{GVF.db} which tells \code{predictCV} how to transform the response of an arbitrary GVF model (\code{'resp'}) into the pivot measure of variability (\code{'CV'}).

By default \code{new.Y = NULL} and CVs (and intervals, if any) obtained by transforming fitted response values will be returned. If passed, argument \code{new.Y} must be a data frame storing new estimates for which CVs have to be predicted. Such input estimates have to be stored in column \code{Y} of data frame \code{new.Y}. Moreover, if \code{object} stores GVF model(s) fitted to \emph{grouped} data (namely, it has class \code{gvf.fit.gr} or \code{gvf.fits.gr}), then \code{new.Y} should also have columns identifying the \emph{groups} to which input estimates are referred (see \sQuote{Examples}). The \emph{only exception} is the following: if the \code{new.Y} data frame has just the \code{Y} column, then CVs will be predicted for \emph{all groups} (see \sQuote{Examples}). The function will check for consistency between groups available in \code{object} and in \code{new.Y} (see \sQuote{Examples}).

If \code{interval = "none"} (the default), the function will return predicted CVs only. Otherwise, lower and upper bounds of confidence (or prediction) intervals around predicted CVs will be also provided. Use argument \code{level} to specify the desired confidence (or tolerance) level for those intervals. 

All the other arguments have the same meaning as in function \code{\link{predict.lm}}.
}

\value{

If \code{object} is a single GVF model (classes \code{gvf.fit} and \code{gvf.fit.gr}), a data frame.

If \code{object} is a set of GVF models fitted to the same data (classes \code{gvf.fits} and \code{gvf.fits.gr}), a list of data frames, one for each input GVF model.

The output data frame(s) will store input estimates \code{new.Y} plus additional columns:
\describe{
    \item{\code{CV.fit}}{Predicted CV value, \code{numeric}.}
    \item{\code{CV.lwr}}{Lower bound of requested interval (if any), \code{numeric}.}
    \item{\code{CV.upr}}{Upper bound of requested interval (if any), \code{numeric}.}
  }

Of course, lower and upper bounds for CVs will be reported only when \code{interval != "none"}.
}

\references{
Wolter, K.M. (2007) \emph{"Introduction to Variance Estimation"}, Second Edition, Springer-Verlag, New York.
}

\author{
Diego Zardetto
}

\note{
Please read the \sQuote{Note} section of \code{\link{predict.lm}}.
}

\seealso{
\code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, and \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it.
}

\examples{
############################################
# Simple examples to illustrate the syntax #
############################################
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

# Not passing 'new.Y' yield CVs from fitted response values:
p <- predictCV(m)

# Take a look:
head(p)
with(p, plot(CV, CV.fit, col = "blue", pch = 20))

# Now let's predict CV values for new estimates of counts
# e.g. Y = c(1000, 5000, 10000, 50000, 100000)
# First, put these values into a data frame:
new.Y <- data.frame(Y = c(1000, 5000, 10000, 50000, 100000))
new.Y

# Then, compute predicted values and confidence intervals:
p <- predictCV(m, new.Y, interval = "confidence")
p

# NOTE: Should we ever need it, we could also use function predict
#       to predict *response* values instead of CVs:
predict(m, new.Y, interval = "confidence")


## (B) A *a set* of GVF models fitted to the same data ##
# Fit example data to registered GVF models for frequencies (i.e. number 1:3):
mm <- fit.gvf(ee.AF, 1:3)

# Let's predict CV values for the same new estimates of counts used above,
# i.e. Y = c(1000, 5000, 10000, 50000, 100000).

# Separate predictions will be obtained from the three fitted GVF models
pp <- predictCV(mm, new.Y, interval = "confidence")
pp

# NOTE: The WARNING above arises from the third fitted GVF model and explains
#       the appearance of NaN at the lower bound of the CV confidence interval
#       for input Y = 100000. Indeed, the response of the third model is the
#       squared CV (which ought to be *positive*), but the prediction for the
#       lower endpoint of the confidence interval happens to be *negative*:
predict(mm, new.Y, interval = "confidence")


## (C) a *single* GVF model fitted to *grouped* data ##
# We have at our disposal the following survey design object on household data:
exdes

# Use function svystat to prepare *grouped* estimates and errors of counts
# to be fitted separately (here groups are regions):
ee.g <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)

# Inspect these grouped estimates and errors of counts:
lapply(ee.g, head)
lapply(ee.g, summary)

# Fit registered GVF model number one, separately inside regions '6', '7',
# and '10':
m.g <- fit.gvf(ee.g, 1)

# Suppose we want to predict CV values for the same new estimates of counts used
# above, i.e. Y = c(1000, 5000, 10000, 50000, 100000).
# Obviously, we need tell to what groups (i.e. regions) should these Y values
# be referred. Therefore, input data frame new.Y should have columns identifying
# the groups (i.e. regions).

# Case 1: all known regions (i.e. regions '6', '7', and '10')
#         Here we can exploit the exception described in section 'Details' and
#         avoid building group variables explicitly:

 # Predict:
p.g <- predictCV(m.g, new.Y, interval = "confidence")
p.g


# Case 2: a subset of known regions (e.g. region '7' only)
#         Here we must build group variables explicitly:
new.Y.g <- data.frame(Y = c(1000, 5000, 10000, 50000, 100000),
                      regcod = 7)
new.Y.g

 # Predict:
p.g <- predictCV(m.g, new.Y.g, interval = "confidence")
p.g


# Case 3: a subset of known regions (e.g. region '7') *plus* some *unknown*
#         region (e.g. region '11').
#         Unknown groups will be tacitly *discarded*:
new.Y.g <- data.frame(Y = c(1000, 5000, 10000, 50000, 100000),
                      regcod = rep(c(7, 11), each = 5))
new.Y.g

 # Predict:
p.g <- predictCV(m.g, new.Y.g, interval = "confidence")
p.g


# Case 4: only *unknown* regions (e.g. regions '11' and '12).
#         This will raise an *error*:
new.Y.g <- data.frame(Y = c(1000, 5000, 10000, 50000, 100000),
                      regcod = rep(c(11, 12), each = 5))
new.Y.g

\dontrun{
 # Predict:
p.g <- predictCV(m.g, new.Y.g, interval = "confidence")
}


# Case 5: *unknown* group variables (e.g. 'region' instead of 'regcod').
#         This will raise an *error*:
new.Y.g <- data.frame(Y = c(1000, 5000, 10000, 50000, 100000),
                      region = rep(c(11, 12), each = 5))
new.Y.g

\dontrun{
 # Predict:
p.g <- predictCV(m.g, new.Y.g, interval = "confidence")
}


## (D) a *set of* GVF models fitted to *grouped* data ##
# Fit all registered GVF models for frequencies (i.e. number 1:3) separately
# inside groups:
mm.g <- fit.gvf(ee.g, 1:3)

# Predict CV values for the same new estimates of counts used above,
# i.e. Y = c(1000, 5000, 10000, 50000, 100000), for all the regions:
# Again, here we can exploit the exception described in section 'Details'
# and avoid building group variables explicitly:

 # Predict:
pp.g <- predictCV(mm.g, new.Y)
pp.g

# NOTE: The WARNING above explains the appearance of NaN for some predicted
#       CV values stemming from the third GVF model. The reason causing this
#       behaviour is exactly the same as discussed in previous example (B).       



####################################################
# Estimating CVs: Prediction vs Direct Calculation #
####################################################
# Load example data:
data(AF.gvf)

# Fit available registered GVF models for frequencies:
mm <- fit.gvf(ee.AF, model=1:3)

# Get the best fitted model:
mbest <- getBest(mm, criterion="adj.R2")
mbest
# Note: adjusted R^2 used as a 'quick and dirty' criterion, as a thorough model
#       comparison via diagnostic plots would have given the same result.

# Compute directly the estimates and errors of a set of absolute frequencies
# which did not belong to the previously fitted data ee.AF, e.g. the joint
# distribution of marstat and regcod:
marstat.regcod <- svystatTM(exdes, ~I(marstat:regcod))
marstat.regcod

# Predict CVs of the joint distribution of marstat and regcod 
# by means of the selected GVF model:
  # First, prepare data with which to predict:
    newdata <- gvf.input(exdes, marstat.regcod)
  # Then, compute CV predictions:
    p.marstat.regcod <- predictCV(mbest, new.Y = newdata, interval="prediction")

# Inspect the results:
p.marstat.regcod

# Plot of computed and predicted CVs with prediction error bars:
  # plot starts #
plot(p.marstat.regcod$Y, p.marstat.regcod$CV.fit, pch=19, col="red",
     ylim=range(p.marstat.regcod$CV.lwr, p.marstat.regcod$CV.upr, p.marstat.regcod$CV),
     xlab="Absolute Frequency Estimate", ylab="Coefficient of Variation",
     main="Estimated and GVF Predicted CVs\n(joint distribution of marstat and regcod)")
segments(x0=p.marstat.regcod$Y, y0=p.marstat.regcod$CV.lwr, y1=p.marstat.regcod$CV.upr,
         col="red")
points(p.marstat.regcod$Y, p.marstat.regcod$CV, pch=0)
legend("topright", title="CV Estimation Method",
       legend=c("Direct Estimate", "GVF Predicted Value", "GVF Prediction Interval"),
       pch=c(0,19,124), col=c("black", "red", "red"), inset=rep(0.05, 2))
  # plot ends #
}

\keyword{survey}
