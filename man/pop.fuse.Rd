\name{pop.fuse}
\alias{pop.fuse}

\title{Fuse Control Totals Data Frames for Special Purpose and Ordinary Calibration Tasks}

\description{
Allows to solve jointly a special purpose calibration task and an ordinary calibration task by fusing their respective control totals.
}

\usage{
pop.fuse(spc.pop, pop, design) 
}

\arguments{
  \item{spc.pop}{A control totals data frame prepared for a \emph{special purpose calibration} task. Must be of class \code{spc.pop}.}
  \item{pop}{A known totals data frame for an \emph{ordinary calibration} task. Must be of class \code{pop.totals}.}
  \item{design}{A design object \emph{prepared} for a special purpose calibration task.}
}

\details{
\pkg{ReGenesees 2.1} introduced support for \emph{\sQuote{special purpose calibration}} tasks, i.e. facilities to calibrate survey weights so as to match complex, non-linear population parameters, instead of ordinary population totals.

Currently, \pkg{ReGenesees}' support for special purpose calibration tasks is limited to Multiple Regression Coefficients (see \code{\link{prep.calBeta}} and \code{\link{pop.calBeta}}), which includes calibration on Means as a notable special case. Further support will likely be provided in future extensions.

Function \code{pop.fuse} allows you to run a calibration task that \emph{simultaneously} involves as benchmarks:
\itemize{
  \item {Complex population parameters (e.g. multiple regression coefficients).}
  \item {Ordinary population totals.}
}

To achieve this goal, \code{pop.fuse} simply \emph{\dQuote{fuses}} the corresponding control totals data frames. The resulting fused control totals data frame is indeed enough to automatically instruct function \code{\link{e.calibrate}} to run the joint calibration task.

Argument \code{spc.pop} must be an object of class \code{spc.pop}, namely a control totals data frame prepared for a special purpose calibration task (e.g. via function \code{\link{pop.calBeta}} for the case of calibration on multiple regression coefficients).

Argument \code{pop} must be an object of class \code{pop.totals}, namely a known totals data frame for an ordinary calibration task (e.g. generated using functions \code{\link{pop.template}} and \code{\link{fill.template}}).

Argument \code{design} must be the survey design object that:\cr
(i)  was already \emph{prepared} for the special purpose calibration task at hand, and\cr
(ii) you want to calibrate simultaneously also on the ordinary population totals at hand.

Note that condition (i) requires that object \code{design} has actually been used to build object \code{spc.pop} (see the \sQuote{Examples} section). Note, moreover, that condition (ii) requires that you eventually run \code{\link{e.calibrate}} on object \code{design} (see the \sQuote{Examples} section).
}

\value{
A \emph{fused} data frame, with class \code{spc.pop}, encompassing control totals for both the special purpose calibration task and the ordinary calibration task.

Note that printing this control totals data frame might not be very telling: to better understand its structure you should instead leverage function \code{\link{pop.desc}}, for which a method dedicated to class \code{spc.pop} is available.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{e.calibrate}} to calibrate weights, functions \code{\link{prep.calBeta}} and \code{\link{pop.calBeta}} to prepare survey data and control totals for calibration on multiple regression coefficients, functions \code{\link{pop.template}} and \code{\link{fill.template}} to generate and fill population totals templates for ordinary calibration tasks, \code{\link{pop.desc}} to obtain a natural language description of control totals data frames.
}

\examples{
# Function pop.fuse allows you to run a calibration task that simultaneously
# involves as benchmarks:
# (A) complex population parameters (e.g. multiple regression coefficients)
# (B) ordinary population totals
# You just have to:
## 1) prepare the survey desing and control totals data frame for (A)
## 2) create and fill the known totals data frame for (B)
## 3) fuse the control totals data frames produced in steps 1) and 2)!

# Load sbs data:
data(sbs)
# Create a design object:
sbsdes <- e.svydesign(data = sbs, ids = ~id, strata = ~strata,
                      weights = ~weight, fpc = ~fpc)

# (A) Suppose you know with satisfactory accuracy from some external source
#     the regression coefficients of the following model:
model <- va.imp2 ~ emp.num + emp.cl

# Here, use the sbs sampling frame available in ReGenesees to simulate the
# external source and compute the values of the regression coefficients:
Beta <- coef(lm(model, data = sbs.frame))
Beta

## 1) Prepare the survey design and control totals for calibration (A):
sbsdes.A <- prep.calBeta(sbsdes, model, Beta)
pop.A <- pop.calBeta(sbsdes.A)

# (B) Suppose you know the number of enterprises and employees by economic
#     activity macro sectors
## 2) Create and fill the known totals data frame for calibration (B)
pop.B <- pop.template(sbsdes.A, calmodel = ~(ent + emp.num):nace.macro - 1)
# Note that, to create the template above, you could have used equally well 
# the original object sbsdes.
pop.B <- fill.template(pop.B, universe = sbs.frame)

## 3) Lastly, fuse the control totals data frames produced in steps 1) and 2)
pop.AB <- pop.fuse(pop.A, pop.B, sbsdes.A)
# Note that, to create the fused control totals data frame above, you *MUST USE*
# the *PREPARED* design object sbsdes.A

# Have a look:
pop.desc(pop.AB)
# ...and recall you can set verbose = TRUE to see the full structure
# pop.desc(pop.AB, verbose = TRUE)

# Now you are ready to calibrate simultaneously on (A) and (B)
sbscal.AB <- e.calibrate(sbsdes.A, pop.AB)
# Note again that, to run the calibration on the *fused* control totals, you
# *MUST USE* the *PREPARED* design object sbsdes.A

# Now, check that all the benchmarks are indeed matched:
# (A) Multiple regression cofficients:
svystatB(sbscal.AB, model)
Beta

# (B) Ordinary population totals:
svystatTM(sbscal.AB, ~ent + emp.num, ~nace.macro)
pop.B
# OK

}
\keyword{survey}