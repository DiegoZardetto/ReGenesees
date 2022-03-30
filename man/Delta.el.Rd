\name{Delta.el}
\alias{Delta.el}
\alias{s1}
\alias{s2}
\docType{data}

\title{
Two Artificial Samples of Elementary Units for Estimation of Change
}

\description{
Two small, partially overlapping datasets, built to mimick non-independent samples selected with a one-stage, stratified, element sampling design. Allow to run \R code contained in the \sQuote{Examples} section of \pkg{ReGenesees} function \code{\link{svyDelta}}.
}

\usage{data(Delta.el)}

\format{
Two data frames, \code{s1} and \code{s2}, with 20 observations each and the following 5 variables.

\itemize{
 \item{For both samples \code{s1} and \code{s2}:}
}
  \describe{
    \item{\code{id}}{Identifier of sample units, \code{numeric}}
    \item{\code{strata}}{Stratification variable, a \code{factor} with \code{2} levels: \code{A}, and \code{B}}
    \item{\code{w}}{Sampling weights, \code{numeric}}
    \item{\code{y}}{A \code{numeric} variable}
    \item{\code{x}}{A \code{numeric} variable, correlated with \code{y}}
  }
}

\details{
The two samples, \code{s1} and \code{s2}, have 8 units in common, resulting in an overlap rate of 8 / 20 = 0.4. One could think of them as, e.g., two consecutive waves of a rotating panel with a 40\% overlap.

Common units are unambigously identified by variable \code{id}.

The stratification is \emph{static}: (1) \code{s1} and \code{s2} use the same strata (i.e. levels \code{A}, and \code{B}), and (2) no common units changed stratum from \code{s1} to \code{s2}.

The \sQuote{Examples} section of \code{\link{svyDelta}} will illustrate the effect of \emph{dynamic} stratification by injecting new strata and stratum-changer units in the samples. 
}

\seealso{
\code{\link{svyDelta}} for calculating estimates and sampling errors of Measures of Change from two not necessarily independent samples, and \code{\link{Delta.clus}} for 2 artificial overlapping samples of clusters.
}

\examples{
data(Delta.el)
# Have a look:
s1
s2

# Have a look at the overlap subsample of 8 units:
sc <- merge(s1, s2, by = "id", suffixes = c("1", "2"))
sc

# Have a look at the full rotation structure (40% overlap in each stratum):
s <- merge(s1, s2, by = "id", all = TRUE, suffixes = c("1", "2"))
s <- s[order(s$strata1, s$strata2), ]
s

# As anticipated, strata are static:
with(s, table(strata1, strata2, useNA = "ifany"))

}

\keyword{datasets}