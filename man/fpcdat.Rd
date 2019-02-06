\name{fpcdat}
\alias{fpcdat}
\docType{data}

\title{
A Small But Not Trivial Artificial Sample Data Set
}

\description{
A small dataset mimicking sample data selected with a 2-stage, stratified, cluster sampling without replacement. Allows to run \R code contained in the \sQuote{Examples} section of the \pkg{ReGenesees} package help pages.
}
\usage{data(fpcdat)}
\format{
  A data frame with 28 observations on the following 12 variables.
  \describe{
    \item{\code{psu}}{Identifier of the primary sampling units, \code{numeric}}
    \item{\code{ssu}}{Identifier of the second stage sampling units, \code{numeric}}
    \item{\code{stratum}}{Stratification Variable, a \code{factor} with \code{5} levels: \code{S.1}, \code{S.2}, \code{S.3}, \code{S.4}, \code{S.5}}
    \item{\code{sr}}{Strata type, \code{integer} with values \code{0} (NSR strata) and \code{1} (SR strata)}
    \item{\code{fpc1}}{First stage finite population corrections, given as population sizes (in terms of psu clusters) inside strata, \code{numeric}}
    \item{\code{fpc2}}{Second stage finite population corrections, given as population sizes (in terms of ssu clusters) inside the corresponding sampled psu, \code{numeric}}
    \item{\code{x}}{A \code{numeric} variable}
    \item{\code{y}}{A \code{numeric} variable}
    \item{\code{dom1}}{A variable defining unplanned estimation domains, \code{factor} with \code{3} levels: \code{A}, \code{B}, \code{C}}
    \item{\code{dom2}}{A variable defining unplanned estimation domains, \code{factor} with \code{6} levels: \code{a}, \code{b}, \code{c}, \code{d}, \code{e}, \code{f}}
    \item{\code{w}}{Direct weights, \code{numeric}}
    \item{\code{z}}{A \code{numeric} variable}
    \item{\code{pl.domain}}{A variable defining planned estimation domains, \code{factor} with \code{3} levels: \code{pd.1}, \code{pd.2}, \code{pd.3}}
  }
}
\details{
Though very small, the \code{fpcdat} dataset concentrates a lot of interesting features. The sampling design is a complex one, with both self-representing (SR) and not-self-representing (NSR) strata. Sampling fractions are deliberately not negligible, in order to stress the effects of finite population corrections on variance estimation. Moreover, being the observations so few, performing computations on the \code{fpcdat} dataset allows to check and understand easily all the effects of setting/changing the global variance estimation options of the \pkg{ReGenesees} package (see e.g. \code{\link{ReGenesees.options}}).
}

\seealso{
\code{\link{ReGenesees.options}} for setting/changing variance estimation options.
}

\examples{
data(fpcdat)
head(fpcdat)
str(fpcdat)
}

\keyword{datasets}