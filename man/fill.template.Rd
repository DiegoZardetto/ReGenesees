\name{fill.template}
\alias{fill.template}

\title{Fill the Known Totals Template for a Calibration Task}

\description{
Given a template prepared to store the totals of the auxiliary variables for a specific calibration task, computes the actual values of such totals from a sampling frame.
}

\usage{
fill.template(universe, template, mem.frac = 10)
}

\arguments{
  \item{universe}{Data frame containing the complete list of the units belonging to the target population, along with the corresponding values of the auxiliary variables (the sampling frame).}
  \item{template}{The template for the calibration task, an object of class \code{pop.totals}.} 
  \item{mem.frac}{A \code{numeric} and non-negative value (the default is \code{10}). It triggers a memory-efficient algorithm when universe is really huge (see \sQuote{Details} and \sQuote{Performance}).}
}

\details{
Recall that a \code{template} object returned by function \code{pop.template} has a structure that complies with the standard required by \code{e.calibrate}, but is \emph{empty}, in the sense that all the known totals it must be able to store are missing (\code{NA}). Whenever these totals are available to the user as such, that is in the form of already computed aggregated values (e.g. because they come from an external source, like a Population Census), the \pkg{ReGenesees} package cannot automatically fill the template. Stated more explicitly: the user himself has to bear the responsibility of putting the \emph{right values} in the \emph{right slots} of the prepared \code{template} data frame. To this end, function \code{\link{pop.desc}} could be very helpful.

A lucky alternative arises when a \emph{"sampling frame"} (that is a data frame containing the complete list of the units belonging to the target population, along with the corresponding values of the auxiliary variables) is available. In such cases, indeed, the \code{fill.template} function is able to: (i) automatically compute the totals of the auxiliary variables from the \code{universe} data frame, (ii) safely arrange and format these values according to the \code{template} structure.

Notice that \code{fill.template} will perform a complete coherence check between \code{universe} and \code{template}. If this check fails, the program stops and prints an error message: the meaning of the message should help the user diagnose the cause of the problem. Should empty levels be present in any factor variable belonging to \code{universe}, they would be dropped.

Argument \code{mem.frac} (whose value must be numeric and non-negative) triggers a memory-efficient algorithm when universe is \emph{really huge}. The \emph{only} sound reason to ever change the value of this argument from its default (\code{mem.frac=10}) is that an invocation of \code{fill.template} caused a memory-failure (i.e. a messages beginning \code{cannot allocate vector of size ...}) on your machine. In such a case, \emph{increasing} the value of \code{mem.frac} (e.g. \code{mem.frac=20}) will provide a better chance of succeeding (for more details, see \sQuote{Performance} section below).
}

\section{Performance}{
Real-world calibration tasks (e.g. in the field of Official Statistics) can simultaneously involve hundreds of auxiliary variables and refer to target populations of several million units. In such circumstances, the naive aggregation of the calibration \code{model.matrix} of \code{universe} may turn out to be too memory-demanding (at least in ordinary PC environments) and determine a memory-failure error.

The alternative implemented in \code{fill.template} is to: (i) split \code{universe} in chunks, (ii) compute partial sums of auxiliary variables chunk-by-chunk, (iii) update \code{template} by adding progressively such partial sums. This alternative is triggered by parameter \code{mem.frac}, which also implicitly controls the number of chunks. The function estimates the memory that would be used to store the \emph{full} \code{model.matrix} of \code{universe} and compares it to the maximum memory allocable on the machine (as returned by \code{\link{memory.limit}}): if the resulting ratio is bigger than \code{1/mem.frac}, the memory-efficient algorithm starts; the number of chunks in which \code{universe} will then be split is determined in such a way that the memory needed to store the \code{model.matrix} of \emph{each} chunk does not exceed a fraction \code{1/mem.frac} of the maximum allocable memory.

Whenever \code{fill.template} switches to the memory-efficient "chunking" algorithm, a warning message will signal it and will specify as well the number of chunks that are being processed.
}

\value{
An object of class \code{pop.totals} storing the \emph{actual} values of the population totals for the specified calibration task, ready to be safely passed to \code{\link{e.calibrate}}.}

\author{Diego Zardetto}

\seealso{
\code{\link{e.calibrate}} for calibrating weights, \code{\link{pop.template}} for the definition of the class \code{pop.totals} and to build a "template" data frame for known population totals, \code{\link{pop.desc}} to provide a natural language description of the template structure, and \code{\link{\%into\%}} for the compression operator for nested factors.
}

\examples{
# Load sbs data:
data(sbs)

# Build a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)


###########################
# A simple example first. #
###########################

# Suppose you want to calibrate on the enterprise counts inside areas
  # 1) Build the population totals template:
pop<-pop.template(sbsdes, calmodel=~area-1)

 # Note: given the dimension of the obtained template...
dim(pop)

 # ...the number of known totals to be stored is 24 (one for each area).
 
 # 2) Use the fill.template function to (i) automatically compute
 #    such 24 totals from the universe (sbs.frame) and (ii) safely fill
 #    the template:
pop<-fill.template(universe=sbs.frame,template=pop) 
pop

 # 3) Lastly calibrate, e.g. with the unbounded linear distance and
 #    heteroskedastic effects proportional to emp.num:
sbscal<-e.calibrate(sbsdes,pop,sigma2=~emp.num,bounds=c(-Inf,Inf)) 


########################################
# A more involved (two-sided) example. #
########################################

# Now suppose you have to perform a calibration process which
# exploits as auxiliary information the total number of employees (emp.num)
# and enterprises (ent) inside the domains obtained by:
#  i) crossing nace2 and region;
# ii) crossing emp.cl, region and nace.macro;

# Due to the fact that nace2 is nested into nace.macro,
# the calibration model can be efficiently factorized as follows:
## 1) Add to the design object and universe the new compressed
 #    factor variable involving nested factors, namely:
sbsdes<-des.addvars(sbsdes,nace2.in.nace.macro=nace2 \%into\% nace.macro)
sbs.frame$nace2.in.nace.macro<-sbs.frame$nace2 \%into\% sbs.frame$nace.macro

  # 2) Build the template exploiting the new variable:
pop<-pop.template(sbsdes,
     calmodel=~(emp.num+ent):(nace2.in.nace.macro + emp.cl)-1,
     partition=~nace.macro:region)

 # Note: given the dimension of the obtained template...
dim(pop)

 # ...the number of known totals to be stored is 792.
 
 # 3) Use the fill.template function to (i) automatically compute
 #    such 792 totals from the universe (sbs.frame) and (ii) safely fill
 #    the template:
pop<-fill.template(universe=sbs.frame,template=pop)

 # Note: out of the 792 known totals in pop, only non-zero entries are actually
 # relevant

 # 4) Lastly calibrate, e.g. with the unbounded linear distance and
 #    heteroskedastic effects proportional to emp.num:
sbscal<-e.calibrate(sbsdes,pop,sigma2=~emp.num,bounds=c(-Inf,Inf))

# Note: a global calibration task would have led to identical calibrated
# weights, but in a more memory-hungry and time-consuming way, as you can
# verify:
  # 1) Build template:
pop.g<-pop.template(sbsdes,
       calmodel=~(emp.num+ent):(nace2:region + emp.cl:nace.macro:region)-1)
dim(pop.g)

  # 2) Fill template:
pop.g <- fill.template(sbs.frame,pop.g)

  # 3) Calibrate globally:
\dontrun{
sbscal.g<-e.calibrate(sbsdes,pop.g,sigma2=~emp.num,bounds=c(-1E6,1E6))

  # 4) Compare calibrated weights (factorized vs. global solution):
range(weights(sbscal)/weights(sbscal.g))

  # ... they are equal.
}


###########################################################
# Just a single example of the memory-efficient algorithm #
# triggered by argument 'mem.frac'.                       #
###########################################################
\dontrun{
 # First artificially increase the size of the sampling frame (e.g.
 # up to 5 million rows):
sbs.frame.HUGE<-sbs.frame[sample(1:nrow(sbs.frame),5000000,rep=TRUE),]
dim(sbs.frame.HUGE)
 
 # Build the template:
pop<-pop.template(sbsdes,
     calmodel=~(emp.num+ent):(nace2.in.nace.macro + emp.cl)-1,
     partition=~nace.macro:region)
dim(pop)

 # Fill the template by using the HUGE universe:
pop<-fill.template(universe=sbs.frame.HUGE,template=pop)
}
}
\keyword{survey}