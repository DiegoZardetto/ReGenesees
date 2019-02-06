\name{find.lon.strata}
\alias{find.lon.strata}
\title{Find Strata with Lonely PSUs}
\description{
Checks whether a stratified design object contains lonely PSUs: if this is the
case, returns the lonely strata levels.
}
\usage{
find.lon.strata(design)
}
\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling
        design metadata.}
}
\details{
Lonely PSUs (i.e. PSUs which are alone inside a not self-representing stratum)
are a concern from the viewpoint of variance estimation. The suggested
\pkg{ReGenesees} facility to handle the lonely PSUs problem is the strata
aggregation technique provided in function \code{\link{collapse.strata}}
(for further alternatives, see also \code{\link{ReGenesees.options}}).

Function \code{find.lon.strata} (originally a private function intended to be
called only by \link{collapse.strata}) is a simple diagnostic tool whose purpose
is to identify the levels of the strata containing lonely PSUs (lonely strata
for short).
}
\value{
The lonely strata levels, if \code{design} actually contains lonely PSUs;
\code{invisible(NULL)} otherwise.
}
\author{Diego Zardetto}
\seealso{
\code{\link{collapse.strata}} for the suggested way of handling lonely PSUs,
\code{\link{ReGenesees.options}} for a different way to face the same problem
(namely by setting variance estimation options), and \code{\link{fpcdat}} for
useful data examples.
}
\examples{
# Load sbs data:
data(fpcdat)

# A negative example first:

  # Build a design object:
    fpcdes<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
            fpc=~fpc1+fpc2)
    fpcdes

  # Find lonely strata:
    find.lon.strata(fpcdes)


# Recall that the difference between certainty PSUs (those sampled with
# probability 1, contained inside self-representing strata) and lonely PSUs
# rests on the fpc information passed to e.svydesign, e.g.:

  # Build a new design object with the same data, now IGNORING fpcs:
    fpcdes.nofpc<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,
                  weights=~w)
    fpcdes.nofpc

  # Find lonely strata:
    find.lon.strata(fpcdes.nofpc)


# A trivial check: collapsing strata eliminates lonely PSUs

  # Apply the collapse strata technique:
    fpcdes.nofpc.clps<-collapse.strata(fpcdes.nofpc)
    fpcdes.nofpc.clps
    clps.strata.status

  # Find lonely strata:
    find.lon.strata(fpcdes.nofpc.clps)

# ...as it must be.
}
\keyword{survey}