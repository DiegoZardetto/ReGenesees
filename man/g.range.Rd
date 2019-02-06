\name{g.range}
\alias{g.range}
\title{Range of g-Weights}
\description{
Computes the range of the ratios between calibrated weights and initial weights (\emph{g-weights}).
}
\usage{
g.range(cal.design)
}
\arguments{
  \item{cal.design}{Object of class \code{cal.analytic}.}
}
\details{
This function computes the smallest interval which contains the ratios between calibrated weights and initial weights.
}
\value{
A numeric vector of length 2.
}
\note{
If \code{cal.design} has undergone \code{k} subsequent calibration steps (with \code{k >= 2}), the function will return the range of the ratios between the output weights of calibration steps \code{k} and \code{k - 1}.
}
\author{Diego Zardetto}
\seealso{
\code{\link{weights}} to extract the weights from a design object, \code{\link{e.calibrate}} for calibrating weights and \code{\link{bounds.hint}} to obtain a hint for calibration problems where range restrictions are imposed on the \emph{g-weights}.
}
\examples{
# Creation of the object to be calibrated:
data(data.examples)
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Calibration (partitioned solution) on the marginal distribution
# of age in 5 classes (age5c) inside provinces (procod)
# (totals in pop06p) with bounds=c(0.5, 1.5):
descal06p<-e.calibrate(design=des,df.population=pop06p,
           calmodel=~age5c-1,partition=~procod,calfun="logit",
           bounds=c(0.5, 1.5),aggregate.stage=2)

# Now let's verify the actual range of the obtained g-weights:
g.range(descal06p)

# which indeed is covered by c(0.5, 1.5), as required.

# Now calibrate once again, this time on the joint distribution of sex
# and marstat (totals in pop03) with the global solution:
descal2<-e.calibrate(design=descal06p,df.population=pop03,
         calmodel=~marstat:sex-1,calfun="linear",bounds=bounds)

# Notice that the print method correctly takes the calibration chain
# into account:
descal2

# The range of the g-weights for the twice calibrated object is:
g.range(descal2)

#... which is equal to:
range(weights(descal2)/weights(descal06p))

#... and must not be confused with:
range(weights(descal2)/weights(des))
}
\keyword{survey}