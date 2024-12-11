\name{n.prop}
\alias{n.prop}
\alias{prec.prop}
\alias{n.comp2prop}
\alias{pow.comp2prop}
\alias{mde.comp2prop}

\title{
Sample Size Requirements and Power Calculations for Proportions
}

\description{
These functions estimate the minimum sample size required to (i) satisfy specific precision constraints in the estimation of proportions and to (ii) attain specified levels of significance and power in a statistical test that compares two proportions. The inverse problems of finding, given a specified sample size, (iii) the expected precision of the estimator of the proportion, and (iv) the expected power or (v) minimum detectable effect for the test that compares two proportions are also addressed.
}

\usage{
n.prop(prec, prec.ind = c("ME", "RME", "SE", "CV"), P = 0.5,
       DEFF = 1, RR = 1,
       F = 1, hhSize = 1, AVEhh = F * hhSize,
       old.clus.size = NULL, new.clus.size = NULL,
       N = NULL, alpha = 0.05, verbose = TRUE)

prec.prop(n, prec.ind = c("ME", "RME", "SE", "CV"), P = 0.5,
          DEFF = 1, RR = 1,
          F = 1, hhSize = 1, AVEhh = F * hhSize,
          old.clus.size = NULL, new.clus.size = NULL,
          N = NULL, alpha = 0.05, verbose = TRUE)

n.comp2prop(P1, P2 = P1, MDE = abs(P2 - P1), K1 = 1/2,
            alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
            pooled.variance = TRUE,
            DEFF = 1, RR = 1,
            F = 1, hhSize = 1, AVEhh = F * hhSize,
            old.clus.size = NULL, new.clus.size = NULL,
            verbose = TRUE)

pow.comp2prop(n, P1, P2 = P1, MDE = abs(P2 - P1), K1 = 1/2,
              alpha = 0.05, sides = c("two-tailed", "one-tailed"),
              pooled.variance = TRUE,
              DEFF = 1, RR = 1,
              F = 1, hhSize = 1, AVEhh = F * hhSize,
              old.clus.size = NULL, new.clus.size = NULL,
              verbose = TRUE)

mde.comp2prop(n, P1, P2 = P1, K1 = 1/2,
              alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
              pooled.variance = TRUE,
              DEFF = 1, RR = 1,
              F = 1, hhSize = 1, AVEhh = F * hhSize,
              old.clus.size = NULL, new.clus.size = NULL,
              verbose = TRUE)
}

\arguments{
  \item{prec}{The precision level you want to attain in estimation.}
  \item{prec.ind}{The precision indicator for which the specified precison level \code{prec} has to be attained. The following uncertainty measures can be used as precision indicators: Margin of Error (\code{'ME'}, the default), Relative Margin of Error (\code{'RME'}), Standard Error (\code{'SE'}), and Coefficient of Variation (\code{'CV'}, also known as Relative Standard Error).}
  \item{P}{Anticipated estimate of the proportion for your target (sub)population.}
  \item{DEFF}{Anticipated estimate of the design effect of the estimator of the proportion.}
  \item{RR}{Anticipated estimate of the response rate.}
  \item{F}{Anticipated estimate of the proportion of your target subpopulation in the general population that will be sampled.}
  \item{hhSize}{Anticipated estimate of the average household size.}
  \item{AVEhh}{Anticipated estimate of the average number of individuals belonging to the target subpopulation per household.}
  \item{old.clus.size}{Average number of households sampled per PSU in the survey you have used to compute the input value for the \code{DEFF} argument.}
  \item{new.clus.size}{Average number of households to be sampled per PSU in the survey you are planning.}
  \item{N}{The size of your target (sub)population. Only needed if you don't want to neglect finite population corrections (which is the default behaviour if \code{N = NULL}).}
  \item{alpha}{Significance level used to build confidence intervals (with confidence level equal to \code{1 - alpha}) or probability of Type I error (False Positive) for hypothesis testing.}
  \item{verbose}{Print on screen some description or just return the results?}
  \item{n}{The sample size for which the expected precision or power or MDE must be calculated (see \sQuote{Details}).}
  \item{P1}{Anticipated estimate of the proportion for group 1.}
  \item{P2}{Anticipated estimate of the proportion for group 2 (only needed when it differs from \code{P1}).}
  \item{MDE}{Minimum Detectable Effect for the difference between proportions.}
  \item{K1}{Ratio between the sample size of group 1 and the total sample size of groups 1 and 2 (defaults to \code{1/2}, i.e. equal sized groups).}
  \item{beta}{Probability of Type II error (False Negative) for hypothesis testing (so that \code{Power = 1 - beta}).}
  \item{sides}{Do you want a \code{'two-tailed'} test (the default) or a \code{'one-tailed'} test?}
  \item{pooled.variance}{Do you want to use the pooled estimator of the variance (the default) or the unpooled one?}
}

\details{
These functions are intended as simple everyday tools for basic sampling design decisions and survey planning. They address sample size requirements and power calculations for \emph{proportions} in surveys (typically household surveys) adopting one- or two-stage sampling designs. The proportions acting as estimation targets can be defined either in terms of individuals or in terms of households. Specific arguments (such as \code{hhSize}, \code{AVEhh}, \code{old.clus.size}, and \code{new.clus.size}) can be tweaked to cover different designs (one or two stages) and target populations (proportions of individuals or households), see the \sQuote{Examples} section.

The formulas needed to guesstimate the minimum sample size required to satisfy specific precision constraints in the estimation of a proportion (function \code{\link{n.prop}}) or the expected precision of the estimator of the proportion given a specified sample size (function \code{\link{prec.prop}}) are available in many textbooks and spreadsheet templates (see, e.g., [Rosner 2006], [Lance, Hattori 2016], [ILO 2014], [UNICEF 2023]).

The formula needed to guesstimate the minimum sample size required to attain specified levels of significance and power in a statistical test that compares two proportions (function \code{\link{n.comp2prop}}) can be derived from equation 10.14 on page 381 of [Rosner 2006], by extending it to real-world complex sampling designs (see, e.g., arguments \code{RR} and \code{DEFF}). Functions \code{\link{pow.comp2prop}} and \code{\link{mde.comp2prop}} solve the resulting equation for \code{Power} (= \code{1 - beta}) and \code{MDE} given a specified \emph{total} sample size, \code{n}.

\strong{NOTE:} For both \code{\link{pow.comp2prop}} and \code{\link{mde.comp2prop}}, argument \code{n} must be the \emph{total} sample size, namely the sum of the sample sizes of the two groups (group 1 and group 2) whose estimated proportions need to be compared. The way the total sample size is allocated to the two groups is controlled by argument \code{K1}.

\strong{NOTE:} Despite most of the arguments to the above five functions come with default values, those defaults must not be considered as suggested values: in general, users will need to specify actual values that suit their specific design needs. For instance, as illustrated in the \sQuote{Examples} section, the Minimum Detectable Effect does not necessarily need to be equal to its default expression (\code{MDE = abs(P2 - P1)}).

\strong{NOTE:} Arguments \code{F}, \code{hhSize}, and \code{AVEhh} are necessarily related to each other, as it is illustrated by the expression of the default of \code{AVEhh}, thus it is not required to specify all of them. If all are explicitly passed, the function will check that they comply with the identity \code{AVEhh = F * hhSize}.

\strong{NOTE:} Most of the arguments to the above five functions are vectorized, meaning that users can pass \code{numeric} vectors to them to investigate multiple design scenarios in a single shot (e.g. via \code{MDE = c(0.06, 0.08, 0.10)}). When vectors of length greater than 1 are passed to multiple arguments, \R recycling rule will kick in.

\strong{NOTE:} Input values for sample size, \code{n}, and population size (if any), \code{N}, should be positive integers. The functions will silently make them so, by using their \code{round(abs( ))} values. When both \code{n} and \code{N} are specified, an error will be raised if \code{n > N}.

\strong{NOTE:} Since, as shown in the \sQuote{Examples} section, the functions can be tweaked to cover different designs (one or two stages) and target populations (proportions of \emph{individuals} or \emph{households}), users who do not want to neglect finite population corrections must specify the population size \code{N} in terms of the \emph{analytical units} that are relevant to their analysis.

\strong{NOTE:} Both functions \code{\link{n.prop}} and \code{\link{n.comp2prop}} return the \emph{ceiling} of the analytically calculated minimum sample sizes (which are not integers in general). This is conservative and will result in slightly better sampling performances than the nominally requested ones (e.g. higher precision or power, and smaller confidence intervals or MDEs).
}

\value{
For \code{n.prop}, \code{prec.prop}, \code{pow.comp2prop}, and \code{mde.comp2prop}, a \code{numeric} vector whose length depends on the length of the numeric inputs passed to the function (under \R recycling rule).

For \code{n.comp2prop}, a list of length 3 with names \code{n1}, \code{n2}, and \code{n}, providing the required sample sizes for group 1, group 2, and the overall sample. Each of the elements of the list is a \code{numeric} vector whose length depends on the length of the numeric inputs passed to the function (under \R recycling rule).
}

\author{
Diego Zardetto
}

\references{
UNICEF (2023) MICS7 TOOLS. URL: \url{https://mics.unicef.org/tools#survey-design}.

ILO (2014) ILO-IPEC TOOLS. URL: \url{https://www.ilo.org/resource/training-material/interactive-tools-sampling-household-based-child-labour-surveys}.

ILO (2014) ILO-IPEC TEMPLATE. URL: \url{https://www.ilo.org/sites/default/files/2024-09/Child_labour_Sampling_Tools_Tool_01_Sample_Size_and_Margin_of_Error_Template_2014.xlsx}.

Rosner, B. (2006) Fundamentals of biostatistics (7th ed.). Boston, MA: Brooks/Cole.

Lance, P., Hattori, A. (2016) Sampling and evaluation. Chapel Hill, North Carolina: MEASURE Evaluation, University of North Carolina (2016).
}

\seealso{
Functions \code{\link{n.mean}}, \code{\link{prec.mean}}, \code{\link{n.comp2mean}}, \code{\link{pow.comp2mean}}, and \code{\link{mde.comp2mean}} provide the same functionalities documented here, when the statistics of interest are \emph{means} rather than proportions.
}

\examples{
##########################
# Reproducible example 1 #
##########################
# Solve the illustrating example problem stored in sheet 'Sample Size(SS) for one domain'
# of the MICS7 Excel template (which can be freely downloaded from:
# https://mics.unicef.org/tools?round=mics7):

n.prop(prec= 0.12, prec.ind= "RME", P= 0.2, DEFF= 1.5, F= 0.09, hhSize= 5.0, RR= 0.9)

# NOTE: We get 3953. The result of the MICS7 template is a bit larger: 4115. This is 
#       simply because the MICS7 template approximates the 1 - 0.05/2 quantile of the
#       standard normal distribution with 2. Indeed, one can see that 4115 is exactly
#       equal to:
floor( 3953 * (2 / 1.96)^2 )


##########################
# Reproducible example 2 #
##########################
# Sheet 'RME and Expected Cases given SS' of the MICS7 Excel template contains one example
# of the reverse problem of finding the expected standard error for a given sample size.
# You can solve it as follows:

prec.prop(n= 4000, prec.ind= "SE", P= 0.2, DEFF= 1.5, F= 0.09, hhSize= 5.0, RR= 0.9)

# NOTE: This time, the result of the MICS7 template (0.01217) is *exactly* the same, since
#       the standard error is a measure of precision that does *not* depend on the
#       confidence level (so that MICS' approximation qnorm(1 - 0.05/2) = 2 plays no role
#       here).


##########################
# Reproducible example 3 #
##########################

# Solve the problem stored as illustrating example in the ILO-IPEC Excel template (which
# can be freely downloaded from the ILO-IPEC TEMPLATE URL in the References section
# above):

n.prop(prec= 0.03, prec.ind= "ME", P= 0.1, DEFF= 4, F= 0.15, hhSize= 5.0, RR= 0.9)

# NOTE: We get 2277. The result of the ILO-IPEC template is a bit larger: 2370. This
#       simply because the ILO-IPEC template too approximates the 1 - 0.05/2 quantile of
#       the standard normal distribution with 2. Indeed, one can see that 2370 is exactly
#       equal to:
floor( 2277 * (2 / 1.96)^2 )


# By using the additional input that the number of households per PSU was 16 for the
# survey that generated the pre-estimates reported above, the ILO-IPEC template also
# returns the estimated intraclass correlation coefficient (ICC) of the binary indicator
# variable underlying the proportion: 27.3%.

# To obtain the same result, just use the arguments 'old.clus.size' and 'new.clus.size'
# as follows:
n.prop(prec= 0.03, prec.ind= "ME", P= 0.1, DEFF= 4, F= 0.15, hhSize= 5.0, RR= 0.9,
       old.clus.size= 16, new.clus.size= 16)

# NOTE: Once more, the result is in agreement with the ILO-IPEC template.

# NOTE: By making explicit the cluster size, we obtained a slightly larger sample size
#       than before (2288 vs 2277). This is because, when the user specifies a value for
#       'new.clus.size', the function returns a sample size that is guaranteed to be an
#       integer multiple of it.

# NOTE: Following Kish's naming convention, ReGenesees uses the label 'ROH' (which stands
#       for 'rate of homogeneity') for the ICC.


####################################################################
# Impact of adjusting the survey design to accommodate a different #
# number of households per sampled PSU.                            #
####################################################################
# More generally, the argument 'new.clus.size' can be used to predict the impact on sample
# size requirements of changing the number of households per PSU. This is performed by
# leveraging the so-called *ICC portability assumption*.

# Suppose the survey you are designing will use 20 households per PSU, rather than the 16
# used by the survey from which you obtained the anticipated estimates used so far. Since
# this setting is expected to increase the DEFF, a larger sample size will be needed to
# the meet the previous precision constraint of ME = 0.03:
n.prop(prec= 0.03, prec.ind= "ME", P= 0.1, DEFF= 4, F= 0.15, hhSize= 5.0, RR= 0.9, 
       old.clus.size= 16, new.clus.size= 20)

# NOTE: As can be seen, when going from 16 to 20 households per PSU, the DEFF is expected 
#       to increase from 4 to 4.8, thus bringing the sample size required to achieve the
#       desired margin of error from 2288 to 2760 households.


##################################################
# Impact of Finite Population Corrections (fpc). #
##################################################
# The sample size above was obtained by neglecting finite population corrections. Now let's
# assume we know the population size is in the ballpark of 100 thousand households. If we
# correctly factor in the fpc terms, then the required sample size is reduced as follows:
n.prop(prec= 0.03, prec.ind= "ME", P= 0.1, DEFF= 4, F= 0.15, hhSize= 5.0, RR= 0.9, 
       old.clus.size= 16, new.clus.size= 20, N= 100000)



#######################################################################
# Power calculations for the comparison of two proportions: example 1 #
#######################################################################
# Suppose one wants to plan two household surveys to evaluate whether a proposed new
# targeting strategy will increase the proportion of beneficiaries of a social assistance
# program that are poor. The surveys will be conducted before (T0) and after (T1) the
# implementation of the new strategy. From a previous assessment (T < T0, with a response
# rate of 0.9), it is known that the proportion of poor households among beneficiaries was
# around 73% and that the DEFF of this estimated proportion was 1.4. The new strategy
# would be considered successful if, after the implementation of the new targeting
# strategy, the proportion of poor households among beneficiaries equals or exceeds 80%.
# Households will be selected from the list of beneficiaries (so F = 1) with the same
# one-stage sampling design as the previous assessment, and the baseline and endline
# surveys will have equal sample sizes (so K1 = 1/2). Moreover, it is believed that the new
# strategy cannot result in harming the targeting performance (so a one-tailed test is
# deemed appropriate). How many households must be sampled to achieve 80% power for the
# detection of the sought-after effect at significance level 5%?

n.comp2prop(P1 = 0.73, P2 = 0.80, sides = "one-tailed", DEFF = 1.4, RR = 0.9)

# NOTE: The required sample size is 1410 households, 705 at baseline (T0) and 705 at
#       endline (T1).

# NOTE: The default value hhSize = 1 (which implies AVEhh = F = 1) was used because the
#       target proportion was defined in terms of households, rather than in terms of
#       individuals.

# Suppose the available budget is not enough to survey the required 705 households per
# round, whereas a sample of 600 households per round (i.e. 1200 overall) would be feasible.
# What would be the expected power in that case?

pow.comp2prop(n = 1200, P1 = 0.73, P2 = 0.80, sides = "one-tailed", DEFF = 1.4, RR = 0.9)

# NOTE: The power would be reduced to about 74.2%.


#######################################################################
# Power calculations for the comparison of two proportions: example 2 #
#######################################################################
# An assistance program aims to reduce by 20% the prevalence of stunted children under age
# 5 in a given target area. The historical prevalence before the intervention is estimated
# at 35%, so bringing it down to 28% (a 20% reduction, resulting in a 7 percentage points
# decrease) would be considered a full success. To evaluate the impact of the intervention,
# two household surveys will be conducted before (T0) and after (T1) the implementation.
# The surveys will have the same sample sizes, and it is desired that the study manages to
# detect a possible reduction effect of 5 percentage points (thus somewhat *smaller* than
# the official target). In the target area, the proportion of children under age 5 in the
# general population is estimated to be 15% and the average household size to be 5.6. The
# latest household survey that estimated the prevalence of stunted children under age 5
# reported a DEFF of 2.70, used a two-stage stratified cluster sampling design with 20
# households per each sampled EA, and had a response rate of 92%. If the same sampling
# design were adopted, how many households would be needed to achieve 80% power for the
# detection of the sought-after minimum effect at significance level 5%?

n.comp2prop(P1 = 0.35, MDE = 0.05, sides = "one-tailed", DEFF = 2.7, RR = 0.92,
            F = 0.15, hhSize = 5.6, old.clus.size = 20, new.clus.size = 20)

# NOTE: The required sample size is 7880 households, 3940 at baseline (T0) and 3940 at 
#       endline (T1). The EAs to be visited are 3940 / 20 = 197 per round.

# What would be the implications of using a design with 15 households per EA, instead of
# 20?

n.comp2prop(P1 = 0.35, MDE = 0.05, sides = "one-tailed", DEFF = 2.7, RR = 0.92,
            F = 0.15, hhSize = 5.6, old.clus.size = 20, new.clus.size = 15)

# NOTE: When going from 20 to 15 households per EA, the DEFF is expected to decrease from
#       2.70 to 2.25. As a result, the required sample size would become 6570 households,
#       3285 at baseline (T0) and 3285 at endline (T1). The EAs to be visited would become
#       3285 / 15 = 219 per round.
#       Therefore, the proposed design would save 655 (= 3940 - 3285) households per round,
#       but would require visiting 22 (= 219 - 197) additional EAs per round.

# Suppose that the available budget only allows to survey 5000 households, 2500 per round,
# and that logistic constraints dictate to survey 15 households per EA. What would be the
# expected MDE of the study in that case, if a 80% power is desired?

mde.comp2prop(n= 5000, P1 = 0.35, sides = "one-tailed", DEFF = 2.7, RR = 0.92,
              F = 0.15, hhSize = 5.6, old.clus.size = 20, new.clus.size = 15)

# NOTE: The MDE would be about 5.7 percentage points. This means that, if the assistance
#       program actually results in the hoped-for 7 percentage points reduction, the study
#       will detect its effect as statistically significant with high probability (>= 80%).
#       However, should the program result in a smaller, yet not negligible, impact (e.g.
#       a reduction of 5 percentage points), the study will have less than 80% probability
#       to detect it, because of the sample size constraint.
}

\keyword{survey}
