\name{drop.gvf.points}
\alias{drop.gvf.points}
\alias{drop.gvf.points.default}
\alias{drop.gvf.points.gvf.fit}
\alias{drop.gvf.points.gvf.fit.gr}

\title{
Drop Outliers and Refit a GVF Model
}

\description{
This function drops observations (alleged outliers) from a fitted GVF model and simultaneously re-fits the model.
}

\usage{
drop.gvf.points(x, method = c("pick", "cut"), which.plot = 1:2,
                res.type = c("standard", "student"), res.cut = 3,
                id.n = 3, labels.id = NULL,
                cex.id = 0.75, label.pos = c(4, 2),
                cex.caption = 1, col = NULL, drop.col = "red",
                ...)
}

\arguments{
  \item{x}{An object containing a \emph{single} fitted GVF model (i.e. of class \code{gvf.fit} or \code{gvf.fit.gr}).}
  \item{method}{\code{character} specifying the method for identifying observations to be dropped
                (see \sQuote{Details}); it may be either \code{'pick'} (the default) or \code{'cut'}.}
  \item{which.plot}{\code{integer} controlling the nature of the plot(s) that are used to identify and/or
                    visualize the observations to be dropped: \code{1} means \sQuote{Observed vs Fitted},
                    \code{2} means \sQuote{Residuals vs Fitted} (see \sQuote{Details}).}
  \item{res.type}{\code{character} specifying what kind of residuals must be used.}
  \item{res.cut}{A positive value: observations to be dropped will be those with residuals whose absolute
                 value exceeds \code{'res.cut'}. Only meaningful if \code{method} is \code{'cut'}.}
  \item{id.n}{Number of points to be initially labelled in each plot, starting with the most extreme.
              Only meaningful if \code{method} is \code{'pick'}.}
  \item{labels.id}{Vector of labels, from which the labels for extreme points will be chosen. \code{NULL} uses observation numbers.}
  \item{cex.id}{Magnification of point labels.}
  \item{label.pos}{Positioning of labels, for the left half and right half of the graph(s) respectively.}
  \item{cex.caption}{Controls the size of \code{caption}.}
  \item{col}{Color to be used for the points in the plot(s).}
  \item{drop.col}{Color to be used to visualize and annotate the points to be dropped in the plot(s).}
  \item{\dots}{Other parameters to be passed through to plotting functions.}
}

\details{
This function drops observations (alleged outliers) from a \emph{single} fitted GVF model and simultaneously re-fits the model.
As a side effect, the function prints on screen the induced change for selected quality measures (see, e.g.,
\code{\link{getR2}}).

If \code{method = "pick"}, observations to be dropped are identified interactively by clicking on points of a plot (see \sQuote{Note}).
Argumemt \code{which.plot} determines the nature of the plot: value \code{1} is for \emph{\sQuote{Observed vs Fitted}},
value \code{2} is for \emph{\sQuote{Residuals vs Fitted}}. In the latter case, argument \code{res.type} specifies what
kind of residuals have to be plotted. Argument \code{id.n} specifies how many points have to be labelled
initially, starting with the most extreme in terms of the selected residuals: this applies to both kinds of plots.

If \code{method = "cut"}, observations to be dropped are those with residuals whose absolute value exceeds the
value of argument \code{res.cut}. Again, argument \code{res.type} specifies what kind of residuals have to be used
(and plotted). The points which have been cut will be highlighted on a plot, whose nature is again specified by
argument \code{which.plot}. If \code{which.plot = 1:2}, dropped points will be visualized on both the
\emph{\sQuote{Observed vs Fitted}} and the \emph{\sQuote{Residuals vs Fitted}} graphs simultaneously.

Argument \code{drop.col} controls the color to be used to visualize and annotate in the plot(s) the points to be
dropped. All the other arguments have the same meaning as in function \code{\link{plot.lm}}.
}

\value{
An object of the same class as \code{x} (i.e. either \code{gvf.fit} or \code{gvf.fit.gr}), containing the original GVF model re-fitted after dropping (alleged) outliers.
}

\author{
Diego Zardetto
}

\note{
For \code{method = "pick"}, function \code{drop.gvf.points} is only supported on those screen devices for which
function \code{\link{identify}} is supported. The identification process can be terminated either by right-clicking the mouse 
and selecting 'Stop' from the menu, or from the 'Stop' menu on the graphics window.
}

\seealso{
\code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Load example data:
data(AF.gvf)

# Inspect available estimates and errors of counts:
str(ee.AF)

# List available registered GVF models:
GVF.db

# Fit example data to registered GVF model number one:
m <- fit.gvf(ee.AF, model=1)
m
summary(m)

##############################################################
# Method 'pick': identify outlier observations to be dropped #
# interactively by clicking on points of a plot.             #
##############################################################
  # Using the 'Observed vs Fitted' plot (the default):
\dontrun{
m1 <- drop.gvf.points(m)
m1
summary(m1)
}

 # Using the 'Residuals vs Fitted' plot with standardized
 # residuals (the default) and increasing id.n to get more
 # labelled points to guide your choices:
\dontrun{
m1 <- drop.gvf.points(m, which.plot = 2, id.n = 10)
m1
summary(m1)
}

 # The same as above, but with studentized residuals and
 # playing with colors:
\dontrun{
m1 <- drop.gvf.points(m, which.plot = 2, id.n = 10, res.type = "student",
                      col = "blue", drop.col = "green", pch = 20)
m1
summary(m1)
}


#############################################################
# Method 'cut': identify outlier observations to be dropped #
# by specifying a threshold for the absolute values of the  #
# residuals.                                                #
#############################################################
 # Using default threshold on standardized residuals and visualizing
 # dropped observations on both 'Observed vs Fitted' and 'Residuals
 # vs Fitted' plots:
m1 <- drop.gvf.points(m, method ="cut")
m1
summary(m1)

 # Using a custom threshold on studentized residuals and visualizing
 # dropped observations on the 'Observed vs Fitted' plot:
m1 <- drop.gvf.points(m, method ="cut",  res.type = "student",
                      res.cut = 2.5, which.plot = 1)
m1
summary(m1)

 # The same as above, but visualizing dropped observations on the
 # 'Residuals vs Fitted' plot:
m1 <- drop.gvf.points(m, method ="cut",  res.type = "student",
                      res.cut = 2.5, which.plot = 2)
m1
summary(m1)

 # You can obviously "cut"/"pick" alleged outliers again from an already
 # "cut"/"picked" fitted GVF model:
m2 <- drop.gvf.points(m1, method ="cut",  res.type = "student",
                      res.cut = 2.5, col = "blue", pch = 20) 
m2
summary(m2)


#################################################################
# Identifying outlier observations to be dropped from "grouped" #
# GVF fitted models (i.e. x has class 'gvf.fit.gr').            #
#################################################################
 # Recall we have at our disposal the following survey design object 
 # defined on household data:
exdes

 # Now use function svystat to prepare "grouped" estimates and errors
 # of counts to be fitted separately (here groups are regions):
ee <- svystat(exdes, y=~ind, by=~age5c:marstat:sex, combo=3, group=~regcod)
ee
plot(ee)

 # Fit registered GVF model number one separately inside groups:
m <- fit.gvf(ee, model=1)
m
summary(m)

 # Now drop alleged outliers separately inside groups:

   #####################################################
   # Method 'pick': work interactively group by group. #
   #####################################################
\dontrun{
   m1 <- drop.gvf.points(m, which.plot = 2, res.type = "student", col = "blue",
                         pch = 20)
   m1
   summary(m1)
}

   #########################################################
   # Method 'cut': apply the same threshold to all groups. #
   #########################################################
   m1 <- drop.gvf.points(m, method ="cut",  res.type = "student", res.cut = 2)
   m1
   summary(m1)

}

\keyword{survey}
