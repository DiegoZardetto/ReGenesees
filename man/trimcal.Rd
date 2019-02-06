\name{trimcal}
\alias{trimcal}

\title{Trim Calibration Weights while Preserving Calibration Constraints}

\description{
This function trims calibration weights to a bounded interval, while preserving all the calibration constraints.
}

\usage{
trimcal(cal.design, w.range = c(-Inf, Inf),
        maxit = 50, epsilon = 1e-07, force = TRUE)
}

\arguments{
  \item{cal.design}{Object of class \code{cal.analytic} containing the calibration weights to be trimmed.}
  \item{w.range}{The interval to which trimmed calibration weights must be bound (see \sQuote{Details}).
                 The default is \code{c(-Inf,Inf)}, which would leave the calibration weights \emph{unchanged}.}
  \item{maxit}{The same as in function \code{\link{e.calibrate}}.}
  \item{epsilon}{The same as in function \code{\link{e.calibrate}}.}
  \item{force}{The same as in function \code{\link{e.calibrate}}.}
}

\details{
Extreme calibration weights might determine unstable estimates and inflate sampling error estimates. To avoid this risk, extreme weights may be \emph{trimmed} by using some suitable procedure (see e.g. [Potter 90], [Valliant, Dever, Kreuter 13]). Despite no rigorous justifications exist for any proposed trimming method, sometimes practitioners are (or feel) compelled to apply a trimming step before estimation. This happens more frequently when interest variables are highly skewed at the population level, like in business surveys or in social surveys with a focus on economic variables (e.g. income, see [Verma, Betti, Ghellini 07], [EUROSTAT 16]).

Unfortunately, the most common trimming techniques do \emph{not} preserve the calibration constraints: if the input weights to the trimming algorithm are calibrated, typically the trimmed weights will \emph{not} reproduce the calibration totals. As a consequence, users have to calibrate again the weights after trimming and iterate through the trimming and calibration steps, until a set of weights is obtained that respects both the trimming bounds and the calibration controls.

Function \code{trimcal} overcomes this limitation: it allows to trim calibration weights to a specific interval while \emph{simultaneously} preserving all the calibration totals. To achieve this result, a \emph{constrained trimming algorithm} is used, which is similar in spirit to the GEM (Generalized Exponential Method) of [Folsom, Singh 2000], but adopts the range restricted \emph{euclidean} distance - instead of the \emph{logit} - for numerical stability considerations.

When \code{w.range} is passed, \emph{both} the trimming limits it defines must be \emph{positive}. In other words, all calibration weights have to be positive after trimming. The purpose of this condition is to enable sound variance estimation on the trimmed object (see also below).

Note that \code{trimcal} is allowed to trim only \emph{already calibrated} weights, i.e. the input object \code{cal.design} must be of class \code{cal.analytic}. This is a deliberate design choice, as trimming design (or initial) weights is methodologically unsound.

Note also that, in case the original calibration weights were asked to be constant within clusters selected at a given sampling stage (via argument \code{aggregate.stage} of \code{\link{e.calibrate}}), \code{trimcal} will \emph{preserve} that property (see \sQuote{Examples}).

Note lastly that \code{trimcal} will not trim further weights that have \emph{just} been trimmed. This is again a deliberate design choice, devised to discourage over-trimming and cosmetic adjustements of the survey weights. Of course, it is instead entirely legitimate to calibrate again a trimmed object: after that, a further trimming step will be allowed.

From a variance estimation perspective, the trimmed object returned by function \code{trimcal} is treated as an ordinary calibrated object. More precisely, the trimming step is regarded as a \emph{\dQuote{finalization}} of the weight adjustment procedure which generated \code{cal.design}, i.e. as a completion of the previous calibration step. Call \code{w}, \code{w.cal} and \code{w.cal.trim} the starting weights, the calibrated weights of \code{cal.design} and the trimmed calibration weights as computed by function \code{trimcal}, respectively. Variance estimates computed on the trimmed object will pretend that one passed from \code{w} to \code{w.cal.trim} \emph{directly} (\code{w -> w.cal.trim}), rather then in two steps (\code{w -> w.cal -> w.cal.trim}). Note incidentally that function \code{\link{get.residuals}}, when invoked on a trimmed object, will behave consistently with the variance estimation approach documented here.
}

\section{Trimming Process Diagnostics}{
Function \code{trimcal} exploits a \emph{constrained} trimming algorithm to adjust the calibration weights so that (i) they fall within a bounded interval but (ii) still preserve all the calibration totals. \emph{When this task is unfeasible, the algorithm will fail}. As a consequence, the adjusted weights returned in the output object will respect the range restrictions set by \code{w.range}, but some of the calibration constraints will be broken. Exactly as for function \code{\link{e.calibrate}}, in order to asses the degree of violation of the calibration constraints introduced by trimming, the user can exploit function \code{\link{check.cal}} (or, equivalently, the diagnostic data structure \code{\link{ecal.status}} available in the \code{.GlobalEnv}).
}   

\value{
A calibrated object of class \code{cal.analytic}, storing \emph{trimmed} calibration weights.
}

\section{Methodological Warning}{
Trimming the calibration weights can result in introducing a \emph{bias} in calibration estimates. Of course, one must hope that this \emph{unknown} bias will turn out to be small compared to the \emph{unknown} gain in precision obtained by trimming. In any case - since the actual effect of trimming weights on the MSE of the estimators is unclear - function \code{trimcal} should be used sparingly and carefully.
}

\author{Diego Zardetto}

\references{
Potter, F.J. (1990) \emph{"A study of procedures to identify and trim extreme sampling weights"}. Proceedings of the Survey Research Methods Section, American Statistical Association, pp. 225-230.

Folsom, R.E., Singh, A.C. (2000) \emph{"The generalized exponential model for sampling weight calibration for extreme values, nonresponse, and poststratification"}. Proceedings of the Section on Survey Research Methods, American Statistical Association, pp. 598-603.

Verma, V., Betti, G., Ghellini, G. (2007) \emph{Cross-sectional and longitudinal weighting in a rotational household panel: applications to EU-SILC}, Statistics in Transition, 8(1), pp. 5-50.

Valliant, R., Dever, J., Kreuter, F. (2013) \emph{"Practical Tools for Designing and Weighting Survey Samples"}. Springer-Verlag, New York.

EUROSTAT (2016) \emph{"EU statistics on income and living conditions (EU-SILC) methodology - data quality"}, \url{http://ec.europa.eu/eurostat/statistics-explained/index.php/EU_statistics_on_income_and_living_conditions_(EU-SILC)_methodology_%E2%80%93_data_quality
 }.
}

\seealso{
\code{\link{e.calibrate}} for calibrating survey weights within \pkg{ReGenesees} and \code{\link{check.cal}} to check if calibration constraints have been fulfilled.
}

\examples{
######################
## Data preparation ##
######################
# Load sbs data:
data(sbs)
# Build a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)
# Build a population totals template and fill it with actual known totals:
pop<-pop.template(sbsdes, calmodel=~(emp.num+ent):emp.cl:nace.macro-1, ~region)
pop <- fill.template(sbs.frame, pop)
# Calibrate:
sbscal <- e.calibrate(sbsdes, pop)
# Have a look at the calibration weights distribution:
summary(weights(sbscal))


#######################
## Trimming examples ##
#######################
## Example 1
# Now suppose we want to trim these calibration weights to, say, the bounded
# interval [0.5, 20]. Let's use our trimcal() function:
sbstrim <- trimcal(sbscal, c(0.5, 20))

# Have a look at the trimmed object:
sbstrim

# Let's first verify that the trimmed calibration weights actually obey the
# imposed range restrictions...
summary(weights(sbstrim))
# ...ok, as it must be.

# Second, let's verify that the trimmed object still preserves all the
# calibration constraints:
check.cal(sbstrim)

# or, more explicitly:
all.equal(aux.estimates(sbscal,  template=pop),
          aux.estimates(sbstrim, template=pop))
# ...ok, as it must be.

# Let's have a look at the scatterplots of calibrated and trimmed weights:
plot(weights(sbsdes), weights(sbscal), pch = 20, col = "red",
     xlab = "Direct weights", ylab = "Calibration (red) and Trimmed (blue) weights")
points(weights(sbsdes), weights(sbstrim), pch = 20, col = "blue")
abline(h = c(0.5, 20), col = "green")

# Last, compute estimates and estimated sampling errors on the trimmed object
# as you would do on ordinary calibrated objects, e.g.
  # before trimming:
svystatTM(sbscal, y = ~va.imp2, by = ~nace.macro, estimator = "Mean")
  # after trimming:
svystatTM(sbstrim, y = ~va.imp2, by = ~nace.macro, estimator = "Mean")


## Example 2
# If w.range is too tight, constrained trimming can fail:
sbstrim2 <- trimcal(sbscal, c(1, 20))

# As a consequence, the trimmed weights will respect the range restrictions...
summary(weights(sbstrim2))

# ...but some of the calibration constraints will be broken:
check.cal(sbstrim2)


## Example 3
# If calibration weights were asked to be constant within clusters, the same
# will hold true for the trimmed calibration weights.

# Load household data:
data(data.examples)

# Define a survey design object:
des <- e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
                   weights=~weight)

# Calibrate asking that all individuals within any household share the same
# calibration weight (re-use an example from ?e.calibrate):
descal <- e.calibrate(design=des,df.population=pop04p,
                      calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
                      bounds=bounds,aggregate.stage=2)

# Have a look at the calibration weights distribution:
summary(weights(descal))

# Trim the calibration weights to, say, the bounded interval [150, 850]
destrim <- trimcal(descal, c(150, 850))

# Do trimmed calibration weights obey the imposed range restrictions?
summary(weights(destrim))

# Verify that trimmed weights are still equal within households:
any( tapply( weights(destrim), destrim$variables$famcod,
             function(x) {length(unique(x)) > 1} ) )

# FALSE, as it must be.


#############################################
## Allowed and forbidden trimming policies ##
#############################################
# Let's illustrate some design restrictions on function trimcal():
# 1) Trimming limits must be both positive:
\dontrun{
trimcal(sbscal, c(-0.05, 18))
}

# 2) Trimming design (or direct, or initial) weights is not allowed, you can
#    only trim calibration weights:
\dontrun{
trimcal(sbsdes, c(1, 18))
}

# 3) You cannot trim further weights that have just been trimmed:
\dontrun{
trimcal(sbstrim, c(1, 18))
}

# 4) You can calibrate again trimmed weights...
pop2<-pop.template(sbsdes, calmodel=~(emp.num+ent):area-1)
pop2<-fill.template(sbs.frame,pop2)
sbscal2<-e.calibrate(sbstrim, pop2) 
summary(weights(sbscal2))

# ...after that, a further trimming step is allowed:
sbstrim2 <- trimcal(sbscal2, c(0.6, 19))
sbstrim2
summary(weights(sbstrim2))
}
\keyword{survey}