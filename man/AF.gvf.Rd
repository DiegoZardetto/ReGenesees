\name{AF.gvf}
\alias{AF.gvf}
%- AF is a list storing estimates and errors of the following statistics:
\alias{AF}
%-- Whole Population Count
%---- ind 
%-- Marginal Counts
%---- sex
%---- marstat
%---- age5c
%---- age10c
%---- regcod
%---- procod
%-- Two-way Counts
%---- sex.marstat
%---- sex.age5c
%---- marstat.age5c
%---- sex.regcod
%---- sex.procod
%-- Three-way Counts
%---- sex.marstat.age5c
%---- sex.marstat.age10c
%---- sex.marstat.regcod
%---- sex.marstat.procod
%-- Four-way Counts
%---- sex.marstat.age5c.regcod
%- The design object used to compute estimates and errors
\alias{exdes}
%- The GVF input built upon estimates and errors above
\alias{ee.AF}
\docType{data}

\title{Example Data for GVF Model Fitting}
\description{
File \code{AF.gvf} containes a set of summary statistics that can be used to illustrate \pkg{ReGenesees} facilities for fitting Generalized Variance Functions models. These summary statistics have kind \sQuote{Absolute Frequency} (see function \code{\link{estimator.kind}}), i.e. involve estimates and errors of counts.
}

\usage{data(AF.gvf)}

\format{
Each row of the \code{ee.AF} data frame represents an estimated absolute frequency along with its estimated sampling error (expressed in terms of standard error, coefficient of variation and variance). The data frame has \code{349} rows, and the following \code{5} columns:
  \describe{
    \item{\code{name}}{The name of the original estimate, \code{factor} with \code{349} levels.}
    \item{\code{Y}}{The value of the original estimate (an absolute frequency), \code{numeric}.}
    \item{\code{SE}}{The standard error of the original estimate, \code{numeric}.}
    \item{\code{CV}}{The coefficient of variation of the original estimate, \code{numeric}.}
    \item{\code{VAR}}{The variance of the original estimate, \code{numeric}.}
  }
}

\details{
Object \code{AF} is a list storing estimates and errors of counts (namely, summary statistics of kind \sQuote{Absolute Frequency}) computed on survey design object \code{exdes}. The names of the slots of list \code{AF} indicate the nature of the corresponding estimates, e.g. element \code{AF[["sex.marstat"]]} stores estimates and errors of the joint absolute frequency distribution of variables \code{sex} and \code{marstat} (see \sQuote{Examples}).

Object \code{ee.AF} is the \code{gvf.input} object built upon all such summary statistics, via function \code{\link{gvf.input}} (see \sQuote{Examples}).
}

\seealso{
\code{\link{estimator.kind}} to assess what kind of estimates are stored inside a survey statistic object, \code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
data(AF.gvf)

# Inspect object AF
class(AF)
length(AF)
names(AF)
AF$sex.marstat
class(AF$sex.marstat)

# Inspect gvf.input object ee.AF
head(ee.AF)
str(ee.AF)
plot(ee.AF)

# The design object used to compute ee.AF is the following:
exdes

# How has object ee.AF been built?
foo <- gvf.input(exdes, stats = AF)
identical(ee.AF, foo)
}

\keyword{datasets}