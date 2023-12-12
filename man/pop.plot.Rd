\name{pop.plot}
\alias{pop.plot}
\alias{pop.plot.pop.totals}
\alias{pop.plot.spc.pop}

\title{Plot Calibration Control Totals vs Current Estimates}

\description{
Draw a scatter plot of calibration control totals vs current estimates.
}

\usage{
pop.plot(pop.totals, ...)

\method{pop.plot}{pop.totals}(pop.totals, design,
                              xlab = "Current Estimates",
                              ylab = "Calibration Control Totals",
                              lcol = c("red", "green", "blue"),
                              lwd = c(1, 1, 1),
                              lty = c(2, 1, 2),
                              verbose = TRUE, ...)

\method{pop.plot}{spc.pop}(pop.totals, design,
                           xlab = "Current Estimates",
                           ylab = "Calibration Control Totals",
                           lcol = c("red", "green", "blue"),
                           lwd = c(1, 1, 1),
                           lty = c(2, 1, 2),
                           verbose = TRUE, ...)
}

\arguments{
  \item{pop.totals}{A known totals data frame for a calibration task. Must be of class \code{pop.totals} (for ordinary calibration tasks) or \code{spc.pop} (for special purpose calibration tasks).}
  \item{design}{A design object to compute estimates of the calibration control totals. Must be of class of class \code{analytic} (or inheriting from it).}
  \item{xlab}{A suggested label for the x axis. See also \code{\link{plot}}.}
  \item{ylab}{A suggested label for the y axis. See also \code{\link{plot}}.}
  \item{lcol}{Colors of reference lines, see \sQuote{Details}. See also \code{\link{par}}.}
  \item{lwd}{Width of reference lines, see \sQuote{Details}. See also \code{\link{par}}.}
  \item{lty}{Type of reference lines, see \sQuote{Details}. See also \code{\link{par}}.}
  \item{verbose}{Print on screen minimum and maximum slopes?}
  \item{\dots}{Other parameters to be passed through to plotting functions.}
}

\details{
Function \code{\link{pop.plot}} draws a scatter plot of calibration control totals vs current estimates derived from object \code{design}. This plot is sometimes very telling and may provide a first-level, rough assessment of how hard the calibration problem at hand will turn out to be in terms of constrained optimization.

A trivial calibration problem would correspond to points lying exactly on the \code{y = x} bisector line. A moderate spread around the bisector line without evident and anomalous patterns is usually a sign of a well behaved calibration task. Larger spreads signal increasing calibration complexity. The maximum and minimum slopes of lines connecting the origin to the points in the scatter plot provide often a reasonable clue about the minimum possible range for the calibration bounds \code{[L, U]}, see \code{\link{bounds.hint}}. Evident and anomalous patterns emerging from the scatter plot not only generally imply harder calibration tasks, but can even reveal signs of underlying biases in the auxiliary variables. For instance, an extremely imbalanced pattern, with the overwhelming majority of points lying below (or above) the \code{y = x} line would be a striking symptom of upward (or downward) bias of the current estimates with respect to the true population totals. 

In addition to the scatter plot of calibration control totals vs current estimates, function \code{\link{pop.plot}} draws three useful reference lines: (1) the maximum slope line connecting the origin to the points, (2) the \code{y = x} bisector line, and (3) the minimum slope line connecting the origin to the points.

If \code{verbose = TRUE}, the minimum and maximum slopes get printed on screen.

When the calibration control totals in \code{pop.totals} span different orders of magnitude (as it often happens in enterprise surveys) it is beneficial to look at the plot in logarithmic scale, which can be achieved by setting \code{log = "xy"} via argument \code{\dots}.

If \code{pop.totals} is of class \code{spc.pop}, i.e. it stores known totals for a special purpose calibration tasks, which are all zero (see \code{\link{prep.calBeta}}), then just one reference vertical line at \code{x = 0} is added to the scatter plot.

If \code{pop.totals} is of class \code{spc.pop} but it is a \emph{fused} data frame, encompassing simultaneous control totals for a special purpose calibration task and an ordinary calibration task (see \code{\link{pop.fuse}}), than two separate scatter plots are drawn.

Note that you can also pass to argument \code{design} the calibrated object obtained using \code{pop.totals} as known totals data frame: if calibration perfectly converged, points will lie on the \code{y = x} bisector line, otherwise the spread around the bisector line will provide a visual assessment of how significantly the calibration task failed.
}

\value{
The main purpose of the function is to draw a plot, anyway it returns invisibly the minimum and maximum slopes whenever it is meaningful to do so (thus excluding non-fused special purpose calibration task).
}

\author{Diego Zardetto}

\seealso{
\code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for the definition of the class \code{pop.totals} and to build a template data frame for known population totals, \code{\link{fill.template}} to automatically fill the template when a sampling frame is available, \code{\link{pop.desc}} to obtain a natural language description of control totals data frames.
}

\examples{
## First prepare some design objects to work with:
# Load household data:
data(data.examples)

# Build a design object:
exdes<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
       weights=~weight)

# Now suppose the known population totals for the calibration task
# are stored in object pop07p:
class(pop07p)
dim(pop07p)

# These totals refer to the 
# 1) Joint distribution of sex and age10c (age in 10 classes)
#    at the region level;
# 2) Joint distribution of sex and age5c (age in 5 classes)
#    at the province level;

# For a natural language description, run:
# pop.desc(pop07p)	   

## Here is the totals vs estimates scatter plot:
pop.plot(pop07p, exdes, pch = 20)

# ...or with a log-log scale:
pop.plot(pop07p, exdes, pch = 20, log = "xy")

# Now calibrate (e.g. with the unbounded linear function):
excal07p <- e.calibrate(exdes, pop07p)

# As calibration converged...
check.cal(excal07p)

# ...the points now lie on the y = x bisector:
pop.plot(pop07p, excal07p, pch = 20)

## You are encouraged to try the function on control totals for special
## purpose calibration tasks (both simple and fused, see ?prep.calBeta,
## ?pop.calBeta, and ?pop.fuse). It might be interesting.

}
\keyword{survey}