\name{pop.desc}
\alias{pop.desc}
\alias{pop.desc.pop.totals}
\alias{pop.desc.spc.pop}

\title{Natural Language Description of Known Totals Templates}

\description{
Provides a natural language description of a known totals data frame to be used for a calibration task.
}

\usage{
pop.desc(pop.totals, ...)

\method{pop.desc}{pop.totals}(pop.totals, ...)

\method{pop.desc}{spc.pop}(pop.totals, verbose = FALSE, ...)
}

\arguments{
  \item{pop.totals}{An object of class \code{pop.totals}, be it a template or the actual known totals data frame for the calibration task.}
  \item{verbose}{Fully describe the control totals of a special purpose calibration task?}
  \item{\dots}{Parameters for future extensions (currently unused).}
}

\details{
Function \code{\link{pop.template}} generates a \emph{template} (i.e. \emph{empty}) data frame of class \code{pop.totals}, which is appropriate to store the known totals of a given calibration task. Afterwards, the template data frame must be filled with actual figures.

When the sampling frame of the survey is available and the actual population totals can be calculated from this source, function \code{\link{fill.template}} (i) automatically computes the totals of the auxiliary variables from the sampling frame, (ii) safely arranges and formats these values according to the \code{template} structure. Therefore, function \code{\link{fill.template}} avoids any need for the user to understand, comply with, or even be aware of, the structure of the template that is being filled.

On the contrary, when the population totals are available to the user as such, that is in the form of already computed aggregated values (e.g. because they come from an external source, like a Population Census), it is up to the user to correctly fill the template, that is to put the \emph{right values} in the \emph{right slots} of the prepared template.

Function \code{pop.desc} has been designed for users who cannot take advantage of function \code{fill.template}, to help them understand the structure of the known totals template, in order to safely fill it with actual figures.

Invoking \code{pop.desc} will print on screen a detailed natural language description of the structure of the input \code{pop.totals} object. Such description will clarify how known totals are organized inside the template slots.
}

\value{
The main purpose of the function is to print on screen, anyway it returns invisibly the input \code{pop.totals} object (as \code{print} would do).
}

\author{Diego Zardetto}

\seealso{
\code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for the definition of the class \code{pop.totals} and to build a template data frame for known population totals, \code{\link{fill.template}} to automatically fill the template when a sampling frame is available.
}

\examples{
## First prepare some design objects to work with:

# Load household data:
data(data.examples)
# Build a design object:
exdes<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
       weights=~weight)

# Load sbs data:
data(sbs)
# Build a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

## Now build some known totals templates that (after having been filled by
## actual figures) could be used to calibrate the design objects above,
## and explore the corresponding natural language description:

  ####################################
  ## Some simple and small examples ##
  ####################################
expop<-pop.template(exdes,calmodel=~1)
expop
pop.desc(expop)

expop<-pop.template(exdes,calmodel=~sex)
expop
pop.desc(expop)

# equivalent to the one above
expop<-pop.template(exdes,calmodel=~sex-1)
expop
pop.desc(expop)

expop<-pop.template(exdes,calmodel=~x1+x2+x3,partition=~procod)
expop
pop.desc(expop)

expop<-pop.template(exdes,calmodel=~sex:marstat-1)
expop
pop.desc(expop)

# equivalent to the one above
expop<-pop.template(exdes,calmodel=~sex*marstat-1)
expop
pop.desc(expop)

# equivalent to the one above
expop<-pop.template(exdes,calmodel=~1,partition=~sex:marstat)
expop
pop.desc(expop)

expop<-pop.template(exdes,calmodel=~x1+age5c:marstat-1,partition=~regcod:sex)
expop
pop.desc(expop)


  ############################################
  ## Some more involved and bigger examples ##
  ############################################
expop<-pop.template(exdes,calmodel=~sex:age10c:regcod + sex:age5c:procod - 1)
expop
pop.desc(expop)

# equivalent to the one above (because procod is nested into regcod and
# age10c is nested into age5c)
expop<-pop.template(exdes,calmodel=~age10c+procod-1,
       partition=~regcod:sex:age5c)
expop
pop.desc(expop)
# NOTE: Most of the entries of the template above will be structural zeros,
#       as can be seen in what follows:
expop.HT<-aux.estimates(exdes, template=expop)
expop.HT
sum(expop.HT==0)


# Switch to sbs data
sbspop<-pop.template(sbsdes,
        calmodel=~(emp.num + ent):(nace.macro + emp.cl) - 1, partition=~region)
 
# Can fill the template using the sampling frame...
sbspop<-fill.template(universe=sbs.frame,template=sbspop)
sbspop

# ...and invoke function pop.desc on the filled known totals data frame:
pop.desc(sbspop)

sbspop <- pop.template(data=sbsdes,
                       calmodel=~((emp.num + ent):(nace2 + emp.cl:nace.macro))-1,
                       partition=~region:public)
sbspop
pop.desc(sbspop)

}
\keyword{survey}