\name{n.mean}
\alias{n.mean}
\alias{prec.mean}
\alias{n.comp2mean}
\alias{pow.comp2mean}
\alias{mde.comp2mean}

\title{
Sample Size Requirements and Power Calculations for Means
}

\description{
These functions estimate the minimum sample size required to (i) satisfy specific precision constraints in the estimation of means and to (ii) attain specified levels of significance and power in a statistical test that compares two means. The inverse problems of finding, given a specified sample size, (iii) the expected precision of the estimator of the mean and (iv) the expected power or (v) minimum detectable effect for the test that compares two means are also addressed.
}

\usage{
n.mean(prec, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
       DEFF = 1, RR = 1,
       F = 1, hhSize = 1, AVEhh = F * hhSize,
       old.clus.size = NULL, new.clus.size = NULL,
       N = NULL, alpha = 0.05, verbose = TRUE)

prec.mean(n, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
          DEFF = 1, RR = 1,
          F = 1, hhSize = 1, AVEhh = F * hhSize,
          old.clus.size = NULL, new.clus.size = NULL,
          N = NULL, alpha = 0.05, verbose = TRUE)

n.comp2mean(sigmaY, MDE, K1 = 1/2,
            alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
            DEFF = 1, RR = 1,
            F = 1, hhSize = 1, AVEhh = F * hhSize,
            old.clus.size = NULL, new.clus.size = NULL,
            verbose = TRUE)

pow.comp2mean(n, sigmaY, MDE, K1 = 1/2,
              alpha = 0.05, sides = c("two-tailed", "one-tailed"),
              DEFF = 1, RR = 1,
              F = 1, hhSize = 1, AVEhh = F * hhSize,
              old.clus.size = NULL, new.clus.size = NULL,
              verbose = TRUE)

mde.comp2mean(n, sigmaY, K1 = 1/2,
              alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
              DEFF = 1, RR = 1,
              F = 1, hhSize = 1, AVEhh = F * hhSize,
              old.clus.size = NULL, new.clus.size = NULL,
              verbose = TRUE)
}

\arguments{
  \item{prec}{The precision level you want to attain in estimation.}
  \item{prec.ind}{The precision indicator for which the specified precison level \code{prec} has to be attained. The following uncertainty measures can be used as precision indicators: Margin of Error (\code{'ME'}, the default), Relative Margin of Error (\code{'RME'}), Standard Error (\code{'SE'}), and Coefficient of Variation (\code{'CV'}, also known as Relative Standard Error).}
  \item{sigmaY}{Anticipated estimate of the standard deviation of interest variable Y for your target (sub)population.}
  \item{muY}{Anticipated estimate of the mean of interest variable Y for your target (sub)population.}
  \item{DEFF}{Anticipated estimate of the design effect of the estimator of the mean.}
  \item{RR}{Anticipated estimate of the response rate.}
  \item{F}{Anticipated estimate of the proportion of your target subpopulation in the general population that will be sampled.}
  \item{hhSize}{Anticipated estimate of the average household size.}
  \item{AVEhh}{Anticipated estimate of the average number of individuals belonging to the target subpopulation per houshold.}
  \item{old.clus.size}{Average number of households sampled per PSU in the survey you have used to compute the input value for the \code{DEFF} argument.}
  \item{new.clus.size}{Average number of households to be sampled per PSU in the survey you are planning.}
  \item{N}{The size of your target (sub)population. Only needed if you don't want to neglect finite population corrections (which is the default behaviour if \code{N = NULL}).}
  \item{alpha}{Significance level used to build confidence intervals (with confidence level equal to \code{1 - alpha}) or probability of Type I error (False Positive) for hypothesis testing.}
  \item{verbose}{Print on screen some description or just return the results?}
  \item{n}{The sample size for which the expected precision or power or MDE must be calculated (see \sQuote{Details}).}
  \item{MDE}{Minimum Detectable Effect for the difference between means.}
  \item{K1}{Ratio between the sample size of group 1 and the total sample size of groups 1 and 2 (defaults to \code{1/2}, i.e. equal sized groups).}
  \item{beta}{Probability of Type II error (False Negative) for hypothesis testing (so that \code{Power = 1 - beta}).}
  \item{sides}{Do you want a \code{'two-tailed'} test (the default) or a \code{'one-tailed'} test?}
}

\details{
These functions are intended as simple everyday tools for basic sampling design decisions and survey planning. They address sample size requirements and power calculations for \emph{means} in surveys (typically household surveys) adopting one- or two-stage sampling designs. The means acting as estimation targets can be defined either in terms of individuals or in terms of households. Specific arguments (such as \code{hhSize}, \code{AVEhh}, \code{old.clus.size}, and \code{new.clus.size}) can be tweaked to cover different designs (one or two stages) and target populations (individual-level or household-level means), see the \sQuote{Examples} section.

The formulas needed to guesstimate the minimum sample size required to satisfy specific precision constraints in the estimation of a mean (function \code{\link{n.mean}}) or the expected precision of the estimator of the mean given a specified sample size (function \code{\link{prec.mean}}) are available in many textbooks and spreadsheet templates (see, e.g., [Rosner 2006], [Lance, Hattori 2016], [ILO 2014a], [ILO 2014b]).

The formula needed to guesstimate the minimum sample size required to attain specified levels of significance and power in a statistical test that compares two means (function \code{\link{n.comp2mean}}) can be derived from equation 8.27 on page 302 of [Rosner 2006], by extending it to real-world complex sampling designs (see, e.g., arguments \code{RR} and \code{DEFF}). Functions \code{\link{pow.comp2mean}} and \code{\link{mde.comp2mean}} solve the resulting equation for \code{Power} (= \code{1 - beta}) and \code{MDE} given a specified \emph{total} sample size, \code{n}.

\strong{NOTE:} For both \code{\link{pow.comp2mean}} and \code{\link{mde.comp2mean}}, argument \code{n} must be the \emph{total} sample size, namely the sum of the sample sizes of the two groups (group 1 and group 2) whose estimated means need to be compared. The way the total sample size is allocated to the two groups is controlled by argument \code{K1}.

\strong{NOTE:} Despite most of the arguments to the above five functions come with default values, those defaults must not be considered as suggested values: in general, users will need to specify actual values that suit their specific design needs. For instance, the default value \code{muY = NULL} is not valid when relative measures of precision are used (\code{'RME'} and \code{'CV'}).

\strong{NOTE:} Arguments \code{F}, \code{hhSize}, and \code{AVEhh} are necessarily related to each other, as it is illustrated by the expression of the default of \code{AVEhh}, thus it is not required to specify all of them. If all are explicitly passed, the function will check that they comply with the identity \code{AVEhh = F * hhSize}.

\strong{NOTE:} Most of the arguments to the above five functions are vectorized, meaning that users can pass \code{numeric} vectors to them to investigate multiple design scenarios in a single shot (e.g. via \code{MDE = c(750, 1000, 1250)}). When vectors of length greater than 1 are passed to multiple arguments, \R recycling rule will kick in.

\strong{NOTE:} Input values for sample size, \code{n}, and population size (if any), \code{N}, should be positive integers. The functions will silently make them so, by using their \code{round(abs( ))} values. When both \code{n} and \code{N} are specified, an error will be raised if \code{n > N}.

\strong{NOTE:} Since, as shown in the \sQuote{Examples} section, the functions can be tweaked to cover different designs (one or two stages) and target populations (\emph{individual-level} or \emph{household-level} means), users who do not want to neglect finite population corrections must specify the population size \code{N} in terms of the \emph{analytical units} that are relevant to their analysis.

\strong{NOTE:} Both functions \code{\link{n.mean}} and \code{\link{n.comp2mean}} return the \emph{ceiling} of the analytically calculated minimum sample sizes (which are not integers in general). This is conservative and will result in slightly better sampling performances than the nominally requested ones (e.g. higher precision or power, and smaller confidence intervals or MDEs).
}

\value{
For \code{n.mean}, \code{prec.mean}, \code{pow.comp2mean}, and \code{mde.comp2mean}, a \code{numeric} vector whose length depends on the length of the numeric inputs passed to the function (under \R recycling rule).

For \code{n.comp2mean}, a list of length 3 with names \code{n1}, \code{n2}, and \code{n}, providing the required sample sizes for group 1, group 2, and the overall sample. Each of the elements of the list is a \code{numeric} vector whose length depends on the length of the numeric inputs passed to the function (under \R recycling rule).
}

\author{
Diego Zardetto
}

\references{
ILO (2014a) ILO-IPEC TOOLS. URL: \url{https://www.ilo.org/resource/training-material/interactive-tools-sampling-household-based-child-labour-surveys}.

ILO (2014b) ILO-IPEC TOOLS USER GUIDE. URL: \url{https://www.ilo.org/sites/default/files/2024-09/Child_labour_Sampling_Tools_Tool_01_Sample_Size_and_Margin_of_Error_Guidelines_2014.pdf}.

Rosner, B. (2006) Fundamentals of biostatistics (7th ed.). Boston, MA: Brooks/Cole.

Lance, P., Hattori, A. (2016) Sampling and evaluation. Chapel Hill, North Carolina: MEASURE Evaluation, University of North Carolina (2016).
}

\seealso{
Functions \code{\link{n.prop}}, \code{\link{prec.prop}}, \code{\link{n.comp2prop}}, \code{\link{pow.comp2prop}}, and \code{\link{mde.comp2prop}} provide the same functionalities documented here, when the statistics of interest are \emph{proportions} rather than means.
}

\examples{
########################
# Reproducible example #
########################
# Solve the problem used as illustrating example in Diagram 1b on page 12 of the ILO-IPEC
# TOOLS USER GUIDE (which can be freely downloaded from ILO-IPEC TOOLS USER GUIDE URL in
# the References section above).
# This is the only example provided by the document for the estimate of a *mean* (which
# the document, somewhat oddly, refers to as an <<indicator of amount>>):

n.mean(prec= 20, prec.ind= "ME", sigmaY= 250, muY= 400, DEFF= 4.0, AVEhh= 0.75, RR= 0.9)

# NOTE: We get 3557. The result of the ILO-IPEC template is a bit larger: 3704. This is
#       simply because the ILO-IPEC template approximates the 1 - 0.05/2 quantile of the
#       standard normal distribution with 2. Indeed, one can see that 3704 is exactly
#       equal to:
round( 3557 * (2 / 1.96)^2 )


##################################################
# Impact of Finite Population Corrections (fpc). #
##################################################
# The sample size above was obtained by neglecting finite population corrections. Now let's
# assume we know the population size is in the ballpark of 100 thousand households. If we
# correctly factor in the fpc terms, then the required sample size is reduced as follows:
n.mean(prec= 20, prec.ind= "ME", sigmaY= 250, muY= 400, DEFF= 4.0, AVEhh= 0.75, RR= 0.9,
       N= 100000)



######################################################
# Power calculations for the comparison of two means #
######################################################
# A new national assistance policy aims to increase by 7% the average monthly household
# consumption, which was estimated to be equal to 362000 (in the country currency) by an
# household survey conducted two years before. The DEFF of that estimate was reported to
# be 3.5, and the same survey estimated the population standard deviation of the monthly
# household consumption to be 461000. To evaluate the impact of the intervention, two
# household surveys with the same sample sizes will be conducted before (T0) and after
# (T1) the implementation of the policy. Assuming that (i) the average monthly household
# consumption at T0 will still be close to 362000 and (ii) the new policy will be the only
# driver of changes, success of the policy would require a value of at least 387340 at T1
# (387340 = 362000 * 1.07), namely an increase of 25340. However,the government would like
# the study to be able to detect a possibly *smaller* consumption increase of 20000.
# Moreover, the government is confident the new policy cannot result in harming household
# consumption (so a one-tailed test is deemed appropriate). Assuming 20 households per EA
# will be sampled as it was the case for the survey conducted two years before, and
# anticipating a response rate of 0.9, and how many households must be selected to achieve
# 80% power for the detection of the sought-after consumption increase at significance
# level 5%?

n.comp2mean(sigmaY = 362000, MDE = 20000, sides = "one-tailed", DEFF = 3.5, RR = 0.9,
            old.clus.size= 20, new.clus.size= 20)

# NOTE: The required sample size is 31520 households, 15760 at baseline (T0) and 15760 at 
#       endline (T1). The EAs to be visited are 15760 / 20 = 788 per round.

# What would be the implications of using designs with 15 or 12 households per EA, instead
# of 20?

n.comp2mean(sigmaY = 362000, MDE = 20000, sides = "one-tailed", DEFF = 3.5, RR = 0.9,
            old.clus.size= 20, new.clus.size= c(15, 12))

# NOTE: When going from 20 to 15 or 12 households per EA, the required sample size is
#       expected to decrease to 25590 households (12795 per round) and 22032 households
#       (11016 per round), respectively. Conversely, the EAs to be visited are expected to
#       increase to 12795 / 15 = 853 and 11016 / 12 = 918 per round,respectively.

# How much sample size could be saved if the government resorts to detecting as significant
# a larger consumption increase of 25000, which is much closer to the actually intended
# goal of the policy?

n.comp2mean(sigmaY = 362000, MDE = 25000, sides = "one-tailed", DEFF = 3.5, RR = 0.9,
            old.clus.size= 20, new.clus.size= 20)

# NOTE: The required sample size would become 20200 households, 10100 per round. Therefore,
#       31520 - 20200 = 11320 households would be saved, 5660 per round.

# The above reduced sample size would, of course, imply a lower statistical power (i.e.
# smaller than the ideally desired 80%) for the detection of an increase of 20000 that
# might actually be induced by the policy. What would actually be the expected power for
# the detection of an effect of 20000 for n = 20200?

pow.comp2mean(n = 20200, sigmaY = 362000, MDE = 20000, sides = "one-tailed", DEFF = 3.5,
              RR = 0.9, old.clus.size= 20, new.clus.size= 20)

# NOTE: The power would be reduced to about 63.5%.
}

\keyword{survey}
