\name{check.cal}
\alias{check.cal}
\title{Calibration Convergence Check}
\description{
Checks whether Calibration Constraints are fulfilled; if not, assesses constraints violation degree.
}
\usage{
check.cal(cal.design)
}
\arguments{
  \item{cal.design}{Object of class \code{cal.analytic}.}
}
\details{
The function verifies if all the imposed Calibration Constraints are actually
fulfilled by object \code{cal.design}. If it is not the case, the function 
evaluates the degree of violation of the constraints and prints a summary of
the mismatches between population totals and achieved estimates (see also Section
'Calibration process diagnostics' in the help page of \code{\link{e.calibrate}}).
}
\value{
The main purpose of the function is to print on screen; anyway a list
is invisibly returned, which summarizes the results of the check.
}
\author{Diego Zardetto}
\seealso{
\code{\link{e.calibrate}} for calibrating weights (in particular, Section 'Calibration process diagnostics').
}
\examples{
# Load sbs data:
data(sbs)

# Build a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

## Example 1
 # Build template...
   pop<-pop.template(sbsdes,~emp.num:emp.cl+ent-1,~region)
 # Fill template...
   pop<-fill.template(sbs.frame,pop)
 # Calibrate...
   sbscal<-e.calibrate(sbsdes,pop,sigma2=~emp.num)
 # Check calibration... 
   check.cal(sbscal)

## Example 2
 # Build template...
   pop<-pop.template(sbsdes,~emp.num+ent-1,~area)
 # Fill template...
   pop<-fill.template(sbs.frame,pop)
 # Calibrate with tight bounds...
   sbscal<-e.calibrate(sbsdes,pop,sigma2=~emp.num,bounds=c(0.8,1.2))
 # Check calibration... 
   check.cal(sbscal)
 
 # Now try to calibrate with suggested bounds... 
   hint <- bounds.hint(sbsdes,pop)
   sbscal<-e.calibrate(sbsdes,pop,sigma2=~emp.num,bounds=hint)
 # Check calibration... 
   check.cal(sbscal)
}
\keyword{survey}