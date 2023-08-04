\name{svyDelta}
\alias{svyDelta}
\alias{details}
\alias{coef.svyDelta}
\alias{SE.svyDelta}
\alias{VAR.svyDelta}
\alias{cv.svyDelta}
\alias{confint.svyDelta}

\title{Estimation of a Measure of Change from Two Not Necessarily Independent Samples}

\description{
Computes estimates and sampling errors of a Measure of Change from two not necessarily independent samples. The Measure of Change can be any analytic function of Horvitz-Thompson or Calibration estimators derived from the two samples.
}

\usage{
svyDelta(expr, design1, design2, by = NULL,
         des.INDEP = FALSE, rho.STRAT = c("Full", "noJump", "noStrat"),
         vartype = c("se", "cv", "cvpct", "var"),
         conf.int = FALSE, conf.lev = 0.95)

details(object, print.call = TRUE, ...)

\method{coef}{svyDelta}(object, ...)
\method{SE}{svyDelta}(object, ...)
\method{VAR}{svyDelta}(object, ...)
\method{cv}{svyDelta}(object, ...)
\method{confint}{svyDelta}(object, ...)
}

\arguments{
  \item{expr}{R expression defining the Measure of Change (see \code{\link{svystatL}} for the basic syntax; see also section \sQuote{Details}).}
  \item{design1}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata related to the first sample. See \sQuote{Details}.}
  \item{design2}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata related to the second sample. See \sQuote{Details}.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population. See \sQuote{Details}.}
  \item{des.INDEP}{Have the samples in design1 and design2 been selected independently? The default is
                 \code{FALSE}. See \sQuote{Details}.}
  \item{rho.STRAT}{To which extent should stratification be considered when estimating correlations? The default is
                 \code{"Full"}, which accommodates even \emph{dynamic stratification} but can sometimes be too computationally demanding. See \sQuote{Details} and \sQuote{Methodological and Computational Remarks} for alternatives.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose
                one or more of: standard error (\code{'se'}, the default), coefficient of variation (\code{'cv'}),
                percent coefficient of variation (\code{'cvpct'}), or variance (\code{'var'}).}
  \item{conf.int}{Compute confidence intervals for the estimates? The default is
                 \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
  \item{object}{An object created by invoking function \code{svyDelta}, about which details are sought.}
  \item{print.call}{Print on screen the call of \code{svyDelta} that generated the input object?}
  \item{\dots}{Additional arguments to \code{coef}, \ldots, \code{confint} methods (if any). For function \code{details}, arguments for future extensions (currently unused).}
}

\details{
Function \code{svyDelta} computes estimates and sampling errors of a Measure of Change from two not necessarily independent samples.

For \sQuote{Measure of Change} we mean here any function that can be used to compare cross-sectional estimates of two population parameters derived from two samples. The easiest Measure of Change is the simple difference between the estimators of the same population parameter referred to two distinct survey occasions. Estimating this difference and its sampling variance is central to hypothesis testing (e.g. for impact evaluation). More complex Measures of Change can of course be considered, for instance nonlinear ones, like the relative difference between two estimators. In addition, the cross-sectional estimators being compared in the Measure of Change can themselves be complex (i.e. expressed as nonlinear functions of Horvitz-Thompson or Calibration estimators). 

Function \code{svyDelta} can accommodate all the cases sketched above.

When the two samples being analyzed are independent, estimating the sampling variance of a simple Measure of Change is fairly straightforward: one just has to add variance estimates derived from the two samples. The problem becomes non-trivial when the two samples are not independent, because sampling covariance terms come into play. Yet, the analysis of non-independent samples is often required in concrete applications, including panel studies (with or without sample rotation), repeated cross-sectional surveys making use of coordinated samples, quasi-experimental designs in which control and treatment groups share some level of sampling information (e.g. observations are taken within same PSUs), etc. 

If the samples in \code{design1} and \code{design2} are not independent, function \code{svyDelta} estimates the needed sampling covariance terms following the method of [Berger, Priam 16]. Note that this method derives design-based covariance estimates from estimates of \emph{correlation} (\code{rho}), which are computed in the first place.

If the Measure of Change is \emph{complex}, \code{svyDelta} automatically linearizes it. Automatic linearization is performed as function \code{\link{svystatL}} would do, along the lines illustrated in [Zardetto, 15], see also section \sQuote{Methodological and Computational Remarks} below.

The mandatory argument \code{expr}, which identifies the Measure of Change, must be an object of class \code{expression}. You can specify just a single Measure of Change at a time, i.e. \code{length(expr)} must be equal to \code{1}. Any analytic function of estimators of Totals and Means derived from the two input design objects \code{design1} and \code{design2} is allowed.

The basic syntax follows the same rules described for function \code{\link{svystatL}}. Inside \code{expr}, the estimator of the Total of a variable is simply represented by the \emph{name} of the variable itself. The reserved name \code{ones} can be used to reference an artificial variable (which will be created on-the-fly, if not already present) whose value is \code{1} for each sampling unit, so that its Total estimator actually estimates the size of the population in terms of elementary units. Therefore, expression \code{y/ones} represents the estimator of the Mean of variable \code{y}. Variables referenced inside \code{expr} must be \code{numeric} and belong to \code{design1} and/or \code{design2}.

To unequivocally associate each variable either to \code{design1} or to \code{design2}, suffixes \code{'.1'} and \code{'.2'} must always be used. This makes it possible to reuse the same variable names across the two samples. For instance, \code{expression(y.1 - y.2)}, represents the difference between the Totals of \code{y} estimated from \code{design1} and \code{design2}, respectively. Similarly, \code{expression((y.1 - y.2) / y.1)} represents the relative difference of the same Totals (using the estimate derived from the first sample as reference value). Along the same vein, the difference between the Means of \code{y} estimated from \code{design1} and \code{design2} can be specified as \code{expression(y.1/ones.1 - y.2/ones.2)}. Likewise, \code{expression(y.1/x.1 - y.2/x.2)} represents the difference between estimates of a Ratio derived from \code{design1} and \code{design2}. 

The mathematical expression of the Measure of Change, as specified by argument \code{expr}, can involve \sQuote{parameters}, that is \emph{symbols} representing given, non-random, scalar, numeric \emph{values}. For each parameter appearing in \code{expr}, the value corresponding to its symbol will be searched following \R standard scoping rules.

Although function \code{svyDelta} is very general and flexible, when the samples are \emph{non-independent} (i.e. \code{des.INDEP = FALSE}, see below), some restrictions exist on arguments \code{design1} and \code{design2}:
\enumerate{
  \item{Objects \code{design1} and \code{design2} must be either both stratified or both unstratified: mixed cases are not allowed.}
  \item{Objects \code{design1} and \code{design2} must be either both element sampling designs or both cluster sampling designs: mixed cases are not allowed.}
  \item{Objects \code{design1} and \code{design2} can be both uncalibrated, both calibrated, or even one uncalibrated and the other calibrated (i.e. mixed cases are allowed).}
}

The method of [Berger, Priam 16] can handle non-independent samples with fixed (i.e. non-random) overlap (including the limiting cases of no overlap and full overlap, for which it reproduces \dQuote{classical} covariance estimates). Therefore, in order for function \code{svyDelta} to work properly, it is a fundamental requirement that the overlapping portion of \code{design1} and \code{design2} can be unequivocally identified. To this end, \strong{the identifiers of sampling units originally used by function} \code{\link{e.svydesign}} \strong{to define objects} \code{design1} \strong{and} \code{design2} \strong{(or, if they are calibrated, their uncalibrated counterparts) must be consistent}. Note that function \code{svyDelta} will use those same identifiers as keys to merge the survey data of \code{design1} and \code{design2} and thereby identify their overlapping subsample. Since, for multi-stage designs, the method of [Berger, Priam 16] resorts to the ultimate cluster approximation, in those cases the fundamental piece of information actually concerns the identifier of the PSUs of \code{design1} and \code{design2}: it is crucial that the latter can be used to identify PSUs that are common to both samples.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{svyDelta} refer to the whole population. If specified, estimation domains must be defined by a formula, following the usual syntactic and semantic rules (see e.g. \code{\link{svystatL}}). Variables referenced in \code{by} must be common to both \code{design1} and \code{design2}, and suffixes \code{'.1'} and \code{'.2'} must not be used for them.

Concerning domain estimation, a methodological caveat is in order. Since the method of [Berger, Priam 16] handles non-independent samples assuming their overlap is fixed, in case \code{design1} and \code{design2} are not independent one should use \code{svyDelta} only for \emph{planned} estimation domains, namely domains that can be obtained by aggregation of sampling strata. In any case, even considering unplanned estimation domains, point estimates produced by \code{svyDelta} will still be correct.

Argument \code{des.INDEP} must be used to tell \code{svyDelta} whether the samples in \code{design1} and \code{design2} were selected independently or not (the default is \code{FALSE}). Note that this argument is indeed necessary, because in \emph{finite population sampling} there is always a chance that two independent samples have a non-null overlap. Therefore, in case of any overlap in \code{design1} and \code{design2}, function \code{svyDelta} will behave differently depending on the value you specify for \code{des.INDEP}:
\itemize{
 \item{If \code{des.INDEP = FALSE}, then any overlap will be regarded as non-accidental and intended by design, thereby triggering the method of [Berger, Priam 16] for calculating the correlation terms needed to estimate the sampling variance of the Measure of Change.}
 \item{If \code{des.INDEP = TRUE}, then any overlap will be regarded as accidental and not intended by design, thereby leading to the sampling variance estimate of the Measure of Change that is appropriate for independent samples.}
}

Note that, when the samples are declared to be \emph{independent} by setting \code{des.INDEP = TRUE}, even the above restrictions 1. and 2. on objects \code{design1} and \code{design2} are obviously lifted, and function \code{svyDelta} can be used to compare estimators coming from arbitrarily mixed sampling designs (e.g. stratified element vs. unstratified cluster).
 
When non-independent samples with sizable overlap are concerned (e.g. for rotating panel surveys), the method of [Berger and Priam 2016] may become \emph{computationally demanding}, especially for highly stratified unit sampling designs (see section \sQuote{Methodological and Computational Remarks}). In such cases, if the invocation of \code{svyDelta} with the default setting \code{rho.STRAT = "Full"} resulted in unaffordable computation time or memory-failure on your machine, you might want to resort to lighter, but approximate, alternatives:
\itemize{
 \item{If \code{rho.STRAT = "noJump"}, then function \code{svyDelta} will disregard the contribution of units that changed stratum from \code{design1} to \code{design2} (if any) when estimating correlation terms (\code{rho}). Note that, in case no such units actually exist, the result will be the \emph{same} that would be obtained with the default setting \code{rho.STRAT = "Full"}, only it will be achieved \emph{using much less computing time and memory}.}
 \item{If \code{rho.STRAT = "noStrat"}, then function \code{svyDelta} will entirely disregard stratification when estimating correlation terms (\code{rho}). This sometimes results in \emph{huge computing time and memory savings}.}
}
In the light of the above, setting \code{rho.STRAT = "noJump"} seems a \emph{very sound option} in many concrete large-scale applications (see section \sQuote{Methodological and Computational Remarks}). 

Note, in any case, that both the approximations \code{rho.STRAT = "noJump"} and \code{rho.STRAT = "noStrat"} are restricted to \emph{correlation} terms only: all the other terms involved in the estimator of the sampling variance of the Measure of Change will still be calculated by properly and fully taking into account the stratification of \code{design1} and \code{design2}.

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.

Function \code{details} provides more details about the results of function \code{svyDelta}. It takes in input an \code{object} returned by \code{svyDelta}, and prints on screen that object together with additional diagnostic information. This additional information (see also section \sQuote{Methodological and Computational Remarks}) is stored in a data frame with one row and the following columns:
\preformatted{
  Column    Meaning
  n1........Sample size of design1 (in terms of elements or PSUs, depending
            on the sampling design)
  n2........Sample size of design2 (in terms of elements or PSUs, depending
            on the sampling design)
  nc........Size of the overlapping sample between design1 and design2
  overlap...Rate of overlap: nc / ( (n1 + n2) / 2 )
  V1........Estimated sampling variance of the linearized Measure of Change,
            arranged in standard difference form, w.r.t. design1: V1(L.1)
  V2........Estimated sampling variance of the linearized Measure of Change,
            arranged in standard difference form, w.r.t. design2: V2(L.2)
  Vind......Estimated sampling variance of the linearized Measure of Change,
            as would be obtained if design1 and design2 were independent: V1 + V2
  rho.......Estimated sampling correlation between L.1 and L.2
  CoV.......Estimated sampling covariance between L.1 and L.2:
            rho * sqrt(V1) * sqrt(V2)
  V.........Estimated sampling variance of the linearized Measure of Change
            w.r.t. design1 and design2 jointly: V1 + V2 - 2 * CoV
}
Note that the same data frame is also invisibly returned by function \code{details}.
}

\value{
An object inheriting from the \code{data.frame} class, whose detailed structure depends on input parameters' values.
}

\section{Methodological and Computational Remarks}{
The method proposed by [Berger and Priam 2016] estimates correlations (\code{rho}) using the estimated residual variance-covariance matrix of a suitable multivariate regression model. The explanatory variables of this multivariate regression model encode information about stratification, and overlap of the two samples. However, the estimation approach remains fully design-based (i.e. not model-based), because the multivariate regression model is exploited only to make calculations easier, whereas its validity is never assumed, nor used, for inferential purposes.

As already mentioned, the method of [Berger and Priam 2016] assumes instead that (i) the size of both samples and, importantly, of their overlapping subsample is fixed, and that (ii) sampling fractions are negligible in both samples. For multi-stage designs, both conditions (i) and (ii) should be fulfilled at the PSU-level.

When stratified designs are concerned, the method is flexible enough to handle \emph{dynamic stratification}. More explicitley, the sample at time 2 can include new strata that were not present at time 1, and sampling units that are common to both samples can belong to different strata at time 1 and 2. However, as already mentioned, this flexibility comes with a computational cost, which can become prohibitive for highly stratified sampling designs.

If no units exist that changed stratum from \code{design1} to \code{design2}, then correlations can be estimated exactly using a more parsimonious multivariate regression model than the one of [Berger and Priam 2016] (i.e. a model with less interaction terms). This is triggered by \code{rho.STRAT = "noJump"}. This option allows for remarkable savings in computing time and memory usage, and typically results in a very good approximation when a moderate number of units actually changed stratum.

In extreme cases, the option of entirely disregarding stratification when estimating correlations (triggered by setting \code{rho.STRAT = "noStrat"}) can sometimes be the only viable solution, at least in ordinary computing environments. Unfortunately, however, the inferential effect of this approximation (e.g. whether it leads to over-estimation or under-estimation) is unclear.

\emph{Complex Measures of Change}, expressed as nonlinear functions of Horvitz-Thompson or Calibration estimators derived from the two samples, are inherently challenging, and would still be so even for independent samples. To estimate the sampling variance of a \emph{complex} Measure of Change, \code{svyDelta} automatically linearizes the estimator specified by argument \code{expression}, call it \strong{Delta}. Automatic linearization is performed as in function \code{\link{svystatL}}, along the lines illustrated in [Zardetto, 15]. Put briefly:
\enumerate{
  \item{Woodruff transforms of the estimator \strong{Delta} are derived with respect to \code{design1} and \code{design2}, respectively.}
  \item{In case any of the two input designs is calibrated, the Woodruff transform associated to it and computed in step 1. includes the appropriate residuals of the involved variables with respect to the calibration model.}
  \item{The linearized expression of the Measure of Change \strong{Delta} is re-arranged in standard difference form, \strong{L.Delta = L.1 - L.2}, where \strong{L.1} and \strong{L.2} are linear functions of the previously computed Woodruff transforms with suitable signs.}
  \item{The method of [Berger and Priam 2016] is applied to \strong{L.Delta} to derive the estimate of the correlation between \strong{L.1} and \strong{L.2} (this is the \code{rho} value that can be extracted using function \code{details}).}
  \item{The other terms of the variance of \strong{L.Delta} are estimated using the respective sampling designs separately}.
  \item{The overall estimate of the variance of \strong{L.Delta} is re-constructed using the terms computed in steps 4. and 5.}
}
}

\author{
Diego Zardetto
}

\references{
Berger, Y. G., Priam, R. (2016). \emph{\dQuote{A simple variance estimator of change for rotating repeated surveys: an application to the European Union Statistics on Income and Living Conditions household surveys}}. Journal of the Royal Statistical Society: Series A (Statistics in Society), 179(1), 251-272.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.

Sarndal, C.E., Swensson, B., Wretman, J. (1992) \emph{\dQuote{Model Assisted Survey Sampling}}, Springer Verlag.
}

\seealso{
Estimators of Complex Analytic Functions of Totals and/or Means \code{\link{svystatL}}. Estimators of Totals and Means \code{\link{svystatTM}}, Ratios between Totals \code{\link{svystatR}}, Shares \code{\link{svystatS}}, Ratios between Shares \code{\link{svystatSR}}, Multiple Regression Coefficients \code{\link{svystatB}}, Quantiles \code{\link{svystatQ}}, and all of the above \code{\link{svystat}}.
}

\examples{
###############################################
# Basic usage and trivial consistency checks. #
###############################################
## Load two non-independent stratified samples of *elementary units*, s1 and s2,
## with 40% overlap (see ?Delta.el)
data(Delta.el)

# Define the design objects
# s1
des1 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = s1)
des1

# s2
des2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = s2)
des2

# Estimate the difference between cross-sectional estimates of the Total of
# variable y
d <- svyDelta(expression(y.2 - y.1), des1, des2)
d

# Now, to know if the difference is statistically significant at alpha = 0.05,
# just see whether the confidence interval covers 0
confint(d)
# it does not, thus the difference is significant

# Have a look at the details, e.g. the estimated correlation rho
details(d)
# ...so rho ~ 0.2

# Can check some of that figures in the details above:
Y1 <- svystatTM(des1, ~y)
Y1
VAR(Y1) # Equal to V1 in details above

Y2 <- svystatTM(des2, ~y)
Y2
VAR(Y2) # Equal to V2 in details above

coef(d)
coef(Y2) - coef(Y1) # Exactly matched, as it must be


## Function svyDelta handles *cluster samples* essentially the same way. The key
## is this time the identifier of PSUs.
## Load two non-independent stratified *cluster samples*, sclus1 and sclus2,
## with 50% overlap at the PSU-level (see ?Delta.clus)
data(Delta.clus)

# Define the design objects
# sclus1
dclus1 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = sclus1)
dclus1

# sclus2
dclus2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = sclus2)
dclus2

# Estimate the difference between cross-sectional estimates of the Total of
# variable y
d <- svyDelta(expression(y.2 - y.1), dclus1, dclus2)
d

# Check statistical significance at alpha = 0.05
confint(d)
# strongly significant (would be so at alpha = 0.01 too)

# Have a look at the details, e.g. the estimated correlation rho
details(d)
# ...so, this time, rho ~ 0.9

# NOTE: If one compared the estimates assuming independence, i.e. disregarding
#       sampling correlations, the difference would be wrongly considered
#       statistically not significant (see the size of Vind above).
#       In fact:
d <- svyDelta(expression(y.2 - y.1), dclus1, dclus2, des.INDEP = TRUE)
details(d)

# Check statistical significance at alpha = 0.05
confint(d)
# non-significant, as anticipated.


###############################
# Complex Measures of Change. #
###############################

## Relative difference (percent) between Totals of y at time 2 and time 1:
d <- svyDelta(expression( 100 * (y.2 - y.1) / y.1 ), des1, des2)
details(d, print.call = FALSE)

## Difference between Means of y at time 2 and time 1:
d <- svyDelta(expression( y.2/ones.2 - y.1/ones.1 ), des1, des2)
details(d, print.call = FALSE)

## Relative difference between Ratios Y/X at time 2 and time 1:
d <- svyDelta(expression( (y.2/x.2 - y.1/x.1) / (y.1/x.1) ), des1, des2)
details(d, print.call = FALSE)

# NOTE: For the cases above, you may want to perform the same checks illustrated
#       for the simple difference between Totals, using suitable ReGenesees
#       functions


#####################################
# Examples with calibrated objects. #
#####################################
# Suppose that:
# (1) Population size was N = 500 at time 1, and N = 501 at time 2
# (2) Total of x was X = 1810 at time 1, and X = 1770 at time 2

# Calibrate on the auxiliary information (1) and (2) above:
# design1
des1 <- des.addvars(des1, ones = 1)
pop1 <- pop.template(des1, calmodel = ~ ones + x - 1)
pop1[,] <- c(500, 1810)
pop1
descal1 <- e.calibrate(des1, pop1)

# design2
des2 <- des.addvars(des2, ones = 1)
pop2 <- pop.template(des2, calmodel = ~ ones + x - 1)
pop2[,] <- c(501, 1770)
pop2
descal2 <- e.calibrate(des2, pop2)

# Compare the estimates of change:
## Calibration estimators
d <- svyDelta(expression(y.2 - y.1), descal1, descal2)
details(d, print.call = FALSE)

## Horvitz-Thompson estimators
d <- svyDelta(expression(y.2 - y.1), des1, des2)
details(d, print.call = FALSE)

# Technically, you can estimate changes in mixed situations, comparing e.g.
# Horvitz-Thompson estimators at time 1 with Calibration estimators at time 2:
## Mix HT and CAL
d <- svyDelta(expression(y.2 - y.1), des1, descal2)
details(d, print.call = FALSE)


########################################################
# Dynamic stratification and computational efficiency. #
########################################################
## Stratification of s1 and s2 is static (see ?Delta.el).
# Therefore the full complexity setting (rho.STRAT = "Full", the default), and
# the alternative and *much more efficient setting* which disregards stratum-
# changer units (if any) in estimating correlations (rho.STRAT = "noJump") will
# yield identical results:

## rho.STRAT = "Full"
d <- svyDelta(expression(y.2 - y.1), des1, des2)
details(d, print.call = FALSE)

## rho.STRAT = "noJump"
d <- svyDelta(expression(y.2 - y.1), des1, des2, rho.STRAT = "noJump")
details(d, print.call = FALSE)

## identical results, as it must be.

## Now check the implications of dynamic stratification.
# Simulate new strata at time 2 (new level "C"):
levels(s2$strata) <- c("A", "B", "C")
# and stratum-changer units from s1 to s2:
s2$strata[4]  <- "C"
s2$strata[14] <- "A"
s2$strata[16] <- "C"
s2$strata[20] <- "C"

# Have a look at the resulting rotation structure of s1 and s2:
s <- merge(s1, s2, by = "id", all = TRUE, suffixes = c("1", "2"))
s <- s[order(s$strata1, s$strata2), ]
s

# Check that strata are dynamic: one unit jumped from "A" to "C", one from "B"
# to "A", and one from "B" to "C":
with(s, table(strata1, strata2, useNA = "ifany"))

# Now recreate object des2 using the updated s2 data
des2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = s2)

# And re-do the comparison:
## rho.STRAT = "Full"
d <- svyDelta(expression(y.2 - y.1), des1, des2)
details(d, print.call = FALSE)

## rho.STRAT = "noJump"
d <- svyDelta(expression(y.2 - y.1), des1, des2, rho.STRAT = "noJump")
details(d, print.call = FALSE)

## The results are slightly different, as expected.

## However, the approximation obtained adopting option rho.STRAT = "noJump"
## is *very good*, and - importantly - *far more efficient* (computation is way
## less memory-hungry and remarkably faster).

## You may want to run the same comparison using real-world stratified samples
## to really appreciate the efficiency gain.
}
\keyword{survey}