\name{Delta.clus}
\alias{Delta.clus}
\alias{sclus1}
\alias{sclus2}
\docType{data}

\title{
Two Artificial Cluster Samples for Estimation of Change
}

\description{
Two small, partially overlapping datasets, built to mimick non-independent PSU samples selected with a one- or multi-stage, stratified, cluster sampling design (but identifiers of SSUs etc. are not reported). Allow to run \R code contained in the \sQuote{Examples} section of \pkg{ReGenesees} function \code{\link{svyDelta}}.
}

\usage{data(Delta.clus)}

\format{
Two data frames, \code{sclus1} and \code{sclus2}, with 6 PSUs each (and 20 and 22 final units, respectively), and the following 5 variables.

\itemize{
 \item{For both samples \code{sclus1} and \code{sclus2}:}
}
  \describe{
    \item{\code{id}}{Identifier of sample PSUs, \code{numeric}}
    \item{\code{strata}}{Stratification variable, a \code{factor} with \code{2} levels: \code{A}, and \code{B}}
    \item{\code{w}}{Sampling weights of final units, \code{numeric}}
    \item{\code{y}}{A \code{numeric} variable}
    \item{\code{x}}{A \code{numeric} variable, correlated with \code{y}}
  }
}

\details{
The two samples, \code{sclus1} and \code{sclus2}, have 3 PSUs in common, resulting in an overlap rate of 3 / 6 = 0.5 at PSU-level. One could think of them as, e.g., two consecutive waves of a rotating panel with a 50\% overlap at PSU-level.

Common PSUs are unambigously identified by variable \code{id}.

The stratification is \emph{static}: (1) \code{sclus1} and \code{sclus2} use the same strata (i.e. levels \code{A}, and \code{B}), and (2) no common PSUs changed stratum from \code{sclus1} to \code{sclus2}.

The \sQuote{Examples} section of \code{\link{svyDelta}} will illustrate the effect of \emph{dynamic} stratification by injecting new strata and stratum-changer units in the samples. 
}

\seealso{
\code{\link{svyDelta}} for calculating estimates and sampling errors of Measures of Change from two not necessarily independent samples, and \code{\link{Delta.el}} for 2 artificial overlapping samples of elementary units.
}

\examples{
data(Delta.clus)
# Have a look:
sclus1
sclus2

# Have a look at the overlap subsample of 3 PSUs (36 final units):
sc <- merge(sclus1, sclus2, by = "id", suffixes = c("1", "2"))
sc

# Have a look at the full rotation structure (50% PSUs overlap in each stratum):
s <- merge(sclus1, sclus2, by = "id", all = TRUE, suffixes = c("1", "2"))
s <- s[order(s$strata1, s$strata2), ]
s

# As anticipated, strata are static:
with(s, table(strata1, strata2, useNA = "ifany"))

}

\keyword{datasets}