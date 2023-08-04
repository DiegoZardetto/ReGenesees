\name{e.svydesign}
\alias{e.svydesign}
\alias{analytic}
\alias{summary.analytic}

\title{Specification of a Complex Survey Design}

\description{
Binds survey data and sampling design metadata.
}

\usage{
e.svydesign(data, ids, strata = NULL, weights,
            fpc = NULL, self.rep.str = NULL, check.data = TRUE)

\method{summary}{analytic}(object, ...)
}

\arguments{
  \item{data}{Data frame of survey data.}
  \item{ids}{Formula identifying clusters selected at subsequent sampling stages (PSUs, SSUs, \ldots).}
  \item{strata}{Formula identifying the stratification variable; \code{NULL} (the default) implies no stratification.}  
  \item{weights}{Formula identifying the initial weights for the sampling units.}
  \item{fpc}{Formula identifying finite population corrections at subsequent sampling stages (see \sQuote{Details}).} 
  \item{self.rep.str}{Triggers an approximate variance estimation method for multistage designs (see \sQuote{Details}). If not \code{NULL} (the default), must be a formula identifying self-representing strata (SR), if any.}
  \item{check.data}{Check out the correct nesting of \code{data} clusters? The default is \code{TRUE}.}

  \item{object}{An object of class \code{analytic}, as returned by \code{e.svydesign}.}
  \item{\dots}{Arguments for future extensions.}
}

\details{
This function has the purpose of binding in an effective and persistent way the survey data to the metadata describing the adopted sampling design. Both kinds of information are stored in a complex object of class \code{analytic}, which extends the \code{survey.design2} class from the \pkg{survey} package. The sampling design metadata are then used to enable and guide processing and analyses provided by other functions in the \pkg{ReGenesees} package (such as \code{\link{e.calibrate}}, \code{\link{svystatTM}}, \ldots).

The \code{data}, \code{ids} and \code{weights} arguments are mandatory, while \code{strata}, \code{fpc}, \code{self.rep.str} and \code{check.data} arguments are optional. The \code{data} variables that are referenced by \code{ids}, \code{weights} and, if specified, by \code{strata}, \code{fpc}, \code{self.rep.str} must not contain any missing value (\code{NA}). Should empty levels be present in any factor variable belonging to \code{data}, they would be dropped.

The \code{ids} argument specifies the cluster identifiers. It is possible to specify a multistage sampling design by simply using a formula which involves the identifiers of clusters selected at subsequent sampling stages. For example, \code{ids=~id.PSU + id.SSU} declares a two-stage sampling in which the first stage units are identified by the \code{id.PSU} variable and second stage ones by the \code{id.SSU} variable.

The \code{strata} argument identifies the stratification variable. The \code{data} variable referenced by \code{strata} (if specified) must be a \code{factor}. By default the sample is assumed to be non-stratified.

The \code{weights} argument identifies the initial (or direct) weights for the units included in the sample. The \code{data} variable referenced by \code{weights} must be \code{numeric}. Direct weights must be strictly positive.

The \code{fpc} formula serves the purpose of specifying the finite population corrections at subsequent sampling stages. By default \code{fpc=NULL}, which implies with-replacement sampling.\cr
If the survey has only one stage, then the \code{fpc}s can be given either as the total population size in each stratum or as the fraction of the total population that has been sampled. In either case the relevant population size must be expressed in terms of sampling units (be they elementary units or clusters). That is, sampling 100 units from a population stratum of size 500 can be specified as 500 or as 100/500=0.2. Thus, passing to \code{fpc} a column of zeros, means again with-replacement sampling.\cr
For multistage sampling the population size (or the sampling fraction) for each sampling stage should also be specified in \code{fpc}. For instance, when \code{ids=~id.PSU + id.SSU} the \code{fpc} formula should look like \code{fpc=~fpc.PSU + fpc.SSU}, with  variable \code{fpc.PSU} giving the population sizes (or sampling fractions) in each stratum for the first stage units, while variable \code{fpc.SSU} gives population sizes (or sampling fractions) for the second stage units in each sampled PSU. Notice that if you choose to pass to \code{fpc} population totals (rather than sampling rates) at a given stage, then you must do the same for all stages (and vice versa).\cr
If \code{fpc} is specified but for fewer stages than \code{ids}, sampling is assumed to be \emph{complete} for subsequent stages. The function will check that \code{fpc}s values at each sampling stage do not vary within strata.

When dealing with a two-stage (multistage) stratified sampling design that includes \emph{self-representing (SR) strata} (i.e. strata containing only one PSU selected with probability 1), the only (leading) contribution to the variance of SR strata arises from the second stage units (\emph{\dQuote{variance PSUs}}).\cr
When \code{options("\link{RG.ultimate.cluster}")} is \code{FALSE} (which is the default for \pkg{ReGenesees}), variance estimation for SR strata is correctly handled provided the survey \code{fpc}s have been properly specified. In particular, if \code{fpc=~fpc.PSU + fpc.SSU} and one specifies \code{fpc}s in terms of sampling fractions, then, inside SR strata, \code{fpc.PSU} must be always equal to one. When, on the contrary, the \emph{\dQuote{Ultimate Cluster Approximation}} holds (i.e. \code{options("\link{RG.ultimate.cluster}")} has been set to \code{TRUE}) the SR strata give no contribution at all to the sampling variance.

A compromise solution (adopted by former existing survey software) is the one of retaining, for both SR and not-SR strata, only the leading contribution to the sampling variance. This means that only the SSUs are relevant for SR strata, whereby only the PSUs matter in not-SR strata. This compromise solution can be achieved by using the \code{self.rep.str} argument. If this argument is actually specified (as a formula referencing the \code{data} variable that identifies the SR strata), a warning is generated in order to remind the user that an approximate variance estimation method will be adopted on that design. Notice that, when choosing the \code{self.rep.str} option, the user must ensure that the variable referenced by \code{self.rep.str} is \code{logical} (with value \code{TRUE} for SR strata and \code{FALSE} otherwise) or \code{numeric} (with value \code{1} for SR strata and \code{0} otherwise) or \code{factor} (with levels \code{"1"} for SR strata and \code{"0"} otherwise).

The optional argument \code{check.data} allows to check out the correct nesting of \code{data} clusters (PSUs, SSUs, \ldots). If \code{check.data=TRUE} the function checks that every unit selected at stage \code{k+1} is associated to one and only one unit selected at stage \code{k}. For a stratified design the function checks also the correct nesting of clusters within strata.
}

\section{PPS Sampling Designs}{
Probability proportional to size sampling \emph{with replacement} does not pose any problem: one must simply specify 
\code{fpc=NULL} and pass the right weights. This holds also for multistage designs, where PSUs are selected with replacement with PPS inside strata. Moreover, when the PSUs are sampled with replacement, the only contribution to the variance arises from the estimated PSU totals, and one can simply ignore any available information about subsequent sampling stages.\cr

For unequal probability sampling \emph{without replacement}, on the contrary, in order to get correct variance estimates, one should know the second-order inclusion probabilities under the sampling design at hand. Unluckily, these probabilities cannot generally be computed, thus one has to resort to some viable approximation. The easier one rests on pretending that PSUs were sampled with replacement, even if this is not actually the case. It is worth stressing that this approach will result in conservative estimates. Moreover, the variance overestimation is expected to be negligible as long as the actual sampling fractions of PSUs are close to zero. Notice that this "with replacement" approximation can be achieved by either not specifying \code{fpc}, or by passing to the PSUs term of \code{fpc} a column of zeros.
}

\value{
An object of class \code{analytic}. The \code{print} method for that class gives a concise description of the sampling design. The \code{summary} method provides further details. Objects of class \code{analytic} persistently store input survey data inside their \code{variables} component. Weights can be accessed by using the \code{\link{weights}} function.
}

\note{
The \code{analytic} class is a specialization of the \code{survey.design2} class from the \pkg{survey} package [Lumley 06]; this means that an object created by \code{e.svydesign} inherits from the \code{survey.design2} class and you can use on it every method defined on the latter class.
}

\author{Diego Zardetto.}

\references{
Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.

Lumley, T. (2006) \emph{\dQuote{survey: analysis of complex survey samples}}, \url{https://CRAN.R-project.org/package=survey}.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\seealso{
\code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}} for calculating estimates and standard errors, \code{\link{e.calibrate}} for calibrating weights, \code{\link{ReGenesees.options}} for setting/changing variance estimation options, \code{\link{collapse.strata}} for the suggested way of handling lonely PSUs, \code{\link{weights}} to extract weights.
}

\examples{
##############################################################
# The following examples illustrate how to create objects    #
# (of class 'analytic') defining different sampling designs. #
# Note: sometimes the same survey data will be used to       #
# define more than one design: this serves only the purpose  #
# of illustrating e.svydesign syntax.                        #
##############################################################

data(data.examples)
# Two-stage stratified cluster sampling design (notice that
# the design contains lonely PSUs):
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~stratum,
     weights=~weight)
des

# Use the summary() function if you need some additional details, e.g.:
summary(des)

# Use the 'variables' slot to extract survey data, e.g.:
head(des$variables)

# Use the weights() function to extract weights, e.g.:
summary(weights(des))

# Again the same design, but using collapsed strata (SUPERSTRATUM variable)
# to remove lonely PSUs:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)
des

# Two stage cluster sampling (no stratification):
des<-e.svydesign(data=example,ids=~towcod+famcod,weights=~weight)
des

# Stratified unit sampling design:
des<-e.svydesign(data=example,ids=~key,strata=~SUPERSTRATUM,
     weights=~weight)
des


data(sbs)
# One-stage stratified unit sampling without replacement
# (notice the presence of the fpc argument):
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
     fpc=~fpc)
des

# Same design as above but ignoring the finite population corrections:
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight)
des


data(fpcdat)
# Two-stage stratified cluster sampling without replacement
# (notice that the fpcs are specified for both stages):
des<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
     fpc=~fpc1+fpc2)
des

# Same design as above but assuming complete sampling for the
# second stage units (notice fpcs have been passed only for the
# first stage):
des<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
     fpc=~fpc1)
des

# Again a two-stage stratified cluster sampling without replacement but
# specified in such a way as to retain, in the estimation phase, only
# the leading contribution to the sampling variance (i.e. the one arising
# from SSUs in SR strata and PSUs in not-SR strata). Notice that the
# self.rep.str argument is used:
des<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
     fpc=~fpc1+fpc2, self.rep.str=~sr)
des
}
\keyword{survey}