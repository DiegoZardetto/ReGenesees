\name{n.mean}
\alias{n.mean}
\alias{prec.mean}
\alias{n.comp2mean}

\title{
Sample Size Requirements and Power Calculations for Means
}

\description{
These functions estimate the minimum sample size required to (i) satisfy specific precision constraints in the estimation of means and to (ii) attain specified levels of significance and power in a statistical test that compares two means. The inverse problem of finding the expected precision of the estimator of the mean given a specified sample size is also addressed.
}

\usage{
n.mean(prec, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
       DEFF = 1, RR = 1,
       F = 1, hhSize = 1, AVEhh = F * hhSize,
       old.clus.size = NULL, new.clus.size = NULL,
       alpha = 0.05, verbose = TRUE)

prec.mean(n, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
          DEFF = 1, RR = 1,
          F = 1, hhSize = 1, AVEhh = F * hhSize,
          old.clus.size = NULL, new.clus.size = NULL,
          alpha = 0.05, verbose = TRUE)

n.comp2mean(sigmaY, MDE, K1 = 1/2,
            alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
            DEFF = 1, RR = 1,
            F = 1, hhSize = 1, AVEhh = F * hhSize,
            old.clus.size = NULL, new.clus.size = NULL,
            verbose = TRUE)
}

\arguments{
  \item{prec}{The precision level you want to attain in estimation.}
  \item{prec.ind}{The precision indicator for which the specified precison level \code{prec} has to be attained. The following uncertainty measures can be used as precision indicators: Margin of Error (\code{'ME'}, the default), Relative Margin of Error (\code{'RME'}), Standard Error (\code{'SE'}), and Coefficient of Variation (\code{'CV'}).}
  \item{sigmaY}{Anticipated estimate of the standard deviation of interest variable Y for your target (sub)population.}
  \item{muY}{Anticipated estimate of the mean of interest variable Y for your target (sub)population.}
  \item{DEFF}{Anticipated estimate of the design effect of the estimator of the mean.}
  \item{RR}{Anticipated estimate of the response rate.}
  \item{F}{Anticipated estimate of the proportion of your target subpopulation in the general population that will be sampled.}
  \item{hhSize}{Anticipated estimate of the average household size.}
  \item{AVEhh}{Anticipated estimate of the average number of individuals belonging to the target subpopulation per houshold.}
  \item{old.clus.size}{Average number of households sampled per PSU in the survey you have used to compute the input value for the \code{DEFF} argument.}
  \item{new.clus.size}{Average number of households to be sampled per PSU in the survey you are planning.}
  \item{alpha}{Significance level used to build confidence intervals (with confidence level equal to \code{1 - alpha}) or probability of Type I error (False Positive) for hypothesis testing.}
  \item{verbose}{Print on screen some description or just return the results?}
  \item{n}{The sample size for which the expected precision must be calculated.}
  \item{MDE}{Minimum Detectable Effect for the difference between means.}
  \item{K1}{Ratio between the sample size of group 1 and the total sample size of groups 1 and 2.}
  \item{beta}{Probability of Type II error (False Negative) for hypothesis testing.}
  \item{sides}{Do you want a \code{'two-tailed'} test (the default) or a \code{'one-tailed'} test?}
}

\details{
Bla bla bla.
}

\value{
Bla bla bla.
}

\author{
Diego Zardetto
}

\seealso{
Bla bla bla.
}

\examples{
# Bla bla bla.
}

\keyword{survey}
