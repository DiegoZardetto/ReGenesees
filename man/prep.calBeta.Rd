\name{prep.calBeta}
\alias{prep.calBeta}
\alias{pop.calBeta}

\title{Calibration on Multiple Regression Coefficients}

\description{
Prepare survey data and control totals to run a calibration task on multiple regression coefficients.
}

\usage{
prep.calBeta(design, model, Beta,
             by = NULL, partition = FALSE, drop.z = TRUE)

pop.calBeta(design) 
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata. For function \code{pop.calBeta}, a \emph{prepared} design object as returned by function \code{prep.calBeta}.}
  \item{model}{Formula (or list of formulas) specifying the linear model(s) whose regression coefficients will be used as calibration benchmarks (see \sQuote{Details}).}
  \item{Beta}{Numeric vector (or list of numeric vectors) of regression coefficients to be used as calibration benchmarks (see \sQuote{Details}).}
  \item{by}{If regression coefficients are known at domain level (i.e. within subpopulations), a formula specifying the domains (see \sQuote{Details}). Specify \code{NULL} (the default option) if the regression coefficients are known at the overall population level.}
  \item{partition}{In case domains have been specified through argument \code{by}, should a partitioned calibration task be performed? The default is \code{FALSE}, which selects a global calibration task.}
  \item{drop.z}{Can the prepared calibration variables (see \sQuote{Details}) be dropped upon completion of the calibration task? Specify \code{TRUE} (the default) if you want to save memory space.}
}
 
\details{
\strong{\emph{Special Purpose Calibration} tasks}\cr
\pkg{ReGenesees 2.1} introduced support for \emph{\sQuote{special purpose calibration}} tasks, i.e. facilities to calibrate survey weights so as to match complex population parameters, instead of ordinary population totals.

When calibration benchmarks come in the form of population parameters that are complex, \emph{non-linear} functions of auxiliary variables (like, in the present case, multiple regression coefficients), calibration constraints have to be linearized first. This generates \emph{synthetic linearized variables} \strong{\emph{z}}, which are the \emph{actual calibration variables} the calibration algorithm will eventually work on. Typically, \emph{control totals} for these synthetic linearized variables \strong{\emph{z}} will \emph{all} be \emph{zero}. See, e.g. [Lesage, 2011].

Put briefly:
\itemize{
  \item {Function \code{prep.calBeta} generates and binds to \code{design} the \emph{synthetic linearized variables} \strong{\emph{z}} needed for the calibration task.}
  \item {Then, given the \emph{prepared} \code{design} object returned by \code{prep.calBeta}, function \code{pop.calBeta} generates the \emph{control totals} data frame for those \strong{\emph{z}} variables.}
}
Of course, once prepared, survey data and control totals returned by the above functions will be fed in input to function \code{\link{e.calibrate}}, which will run the calibration task.\cr\cr

\strong{Function \code{prep.calBeta()}}\cr
Function \code{prep.calBeta} makes it possible to:
\enumerate{
  \item {Calibrate on regression coefficients of \emph{different linear models}, each known at the \emph{overall population level}.}
  \item {Calibrate on regression coefficients of a \emph{single linear model}, known possibly for \emph{different subpopulations}.}
}
Note that, as detailed below, the key argument to switch from case 1. to case 2. is \code{by}. For extensive illustration of both the above cases 1. and 2., please see the \sQuote{Examples} section.

The mandatory argument \code{model} specifies the linear model(s) whose regression coefficients will be used as calibration benchmarks. Use a \code{formula} object to specify a \emph{single} linear model, or a \code{list} of \code{formula} objects to specify \emph{several} different linear models. For details on model specification, see e.g. \code{\link{lm}}.\cr
The \code{design} variables referenced by \code{model} formula(s) must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}).

The mandatory argument \code{Beta} specifies the vector(s) of regression coefficients that will be used as calibration benchmarks. Use a \code{numeric} vector to specify regression coefficients of a \emph{single} linear model, or a \code{list} of \code{numeric} vectors to specify regression coefficients of \emph{several} different linear models. The \code{Beta} vector(s) must not contain any missing value (\code{NA}).

If argument \code{model} is passed as a list, argument \code{Beta} must be passed as a list too, and the \code{length} of both must be the same. If this is the case, \code{Beta} vectors will be \emph{positionally} tied to linear model formulas contained in \code{model}.

Each vector of regression coefficients specified through \code{Beta} must be \emph{consistent} with the corresponding \code{model}, namely match its model matrix columns (as generated using actual \code{design} data). Function \code{prep.calBeta} will check for this consistency and raise an informative error message in case of failure.\cr
Note that, in order to ensure consistency with \code{model}, not only the \code{length}, but also the \emph{order} of \code{Beta} elements \emph{matters}. If in doubt, you can easily learn about the \emph{right} ordering of \code{Beta} coefficients, given \code{model}, by calling function \code{\link{svystatB}} as follows: \code{svystatB(design, model)}. This will, of course, return the estimated regression coefficients of \code{model} \emph{before} calibration.\cr
Note that each \code{Beta} vector can have \code{names} or not. If a \code{Beta} vector has \emph{names} that match the expected ones (given the corresponding \code{model} formula), but appear in a \emph{different order}, then function \code{prep.calBeta} will suitably re-order them and inform you with a \code{warning} message.

As anticipated, argument \code{by} can be used to enable calibration on regression coefficients known at \emph{domain level} (i.e. within subpopulations).

If passed, \code{by} must be a formula: for example, the statement \code{by=~B1:B2} defines as domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. Notice that a formula like \code{by=~B1+B2} will be automatically translated into the factor-crossing formula \code{by=~B1:B2}. The \code{design} variables referenced by \code{by} (if any) should be of type \code{factor}, otherwise they will be coerced. Note that, to prevent obvious collinearity issues, the variables referenced by argument \code{by} must \emph{not} appear in the input \code{model} formula: otherwise, the program will stop and print an error message.

If you specify domains through argument \code{by}, you will be allowed to specify \emph{just a single} linear model through argument \code{model}. Instead, through argument \code{Beta}, you will have to specify domain-level vectors of regression coefficients. Therefore, \code{Beta} will be necessarily passed as a \emph{list}. Note that, in this case, the \code{Beta} list must have as many components as the domains defined through argument \code{by}, and the two will be matched \emph{positionally} by function \code{prep.calBeta}. Therefore, the \emph{order} of elements in the \code{Beta} list \emph{matters}. Specifically, \code{Beta} vectors must appear in the same order as the domains specified by argument \code{by}. You can easily learn about the \emph{right} ordering of \code{by} domains (and hence \code{Beta} elements) by calling function \code{\link{svystatB}} as follows: \code{svystatB(design, model, by)}. This will, of course, return the estimated regression coefficients of \code{model} within \code{by} domains \emph{before} calibration.

Lastly, if you specify domains through argument \code{by}, you will be allowed to decide whether \code{design} should be prepared for a \emph{global} or a \emph{partitioned} calibration task (see \code{\link{e.calibrate}}). Recall that partitioned calibration tasks are computationally more efficient, but produce exactly the same results. Note that, if you select \code{partition = TRUE}, than the \emph{calibration domains} used to split the global calibration task will be the same domains specified via argument \code{by}. More explicitly: function \code{\link{e.calibrate}} will eventually be called on \code{prep.calBeta}'s output object silently assuming \code{partition = by}.

The synthetic linearized variables \strong{\emph{z}} prepared by \code{prep.calBeta} will eventually be used by function \code{\link{e.calibrate}} to solve the calibration task. Argument \code{drop.z} allows you to instruct \code{\link{e.calibrate}} to drop these variables from its output object - or rather keep them within it - upon completion of the calibration task. The default has been set to \code{TRUE} to reduce memory usage. In case you want to be able to consistently trim weights after calibration via function \code{\link{trimcal}}, you \emph{must} specify \code{drop.z = FALSE}.\cr\cr

\strong{Function \code{pop.calBeta()}}\cr
Given the \emph{prepared} \code{design} object returned by \code{prep.calBeta}, function \code{pop.calBeta} generates the \emph{control totals} data frame needed to run the calibration task. This dataframe suitably accomodates the control totals of the \emph{synthetic linearized variables} \strong{\emph{z}} prepared by \code{prep.calBeta}.

Note that the data frame object returned by \code{pop.calBeta} is \emph{already filled}, with control totals that are \emph{all zero}, and it is ready to be directly passed to function \code{\link{e.calibrate}} (see the \sQuote{Examples} section).

Note, lastly, that printing this control totals data frame might not be very telling: to better understand its structure you should instead leverage function \code{\link{pop.desc}}, for which a method dedicated to \emph{\sQuote{special purpose calibration}} tasks is available (see the \sQuote{Examples} section).
}

\section{Linear Algebra Remark}{
Contrary to ordinary calibration tasks, the control totals of a \sQuote{special purpose calibration} task are typically all zero. Therefore, calibration constraints of such tasks - unlike ordinary ones - define a system of linear equations that is \emph{homogeneous}. Homogeneous systems always admit the trivial zero solution, which, in calibration terms, would mean output weights that are identically zero (and thus statistically unacceptable). For this reason, the possibility that function \code{\link{e.calibrate}} ends up with a \emph{false convergence} of the special purpose calibration task cannot, in general, be ruled out. Note that functions \code{\link{e.calibrate}} and \code{\link{check.cal}} are able to detect such false convergence events and warn the user about it.

A satisfactory countermeasure to this issue is to set calibration bounds to any interval that does not include zero (see argument \code{bounds} of \code{\link{e.calibrate}}).\cr
A second, very attractive alternative would be to run the special purpose calibration task \emph{jointly} with some ordinary calibration task, as the joint calibration constraints would then define a \emph{non-homogeneous} system. This solution can be obtained straightforwardly using function \code{\link{pop.fuse}}.
}

\value{
\itemize{
  \item {For function \code{pop.calBeta}, an object of the same class as \code{design}, storing the freshly created \emph{synthetic linearized variables} \strong{\emph{z}} as columns of its \code{$variables} slot (see the \sQuote{Examples} section).}
  \item {For function \code{pop.calBeta}, a \emph{control totals} data frame for those \strong{\emph{z}} variables, with class \code{spc.pop} (see the \sQuote{Examples} section).}
}
}

\author{
Diego Zardetto
}

\references{
Lesage, E. (2011). \emph{\dQuote{The use of estimating equations to perform a calibration on complex parameters}}. Survey Methodology. 37.
}

\seealso{
\code{\link{e.calibrate}} to calibrate weights, \code{\link{svystatB}} to compute estimates and sampling errors of Multiple Regression Coefficients, \code{\link{pop.desc}} to obtain a natural language description of control totals (including those for \emph{special purpose calibration} tasks), \code{\link{pop.fuse}} to fuse population totals prepared for ordinary and special purpose calibration tasks.
}

\examples{
# Load sbs data:
data(sbs)
# Create a design object:
sbsdes <- e.svydesign(data = sbs, ids = ~id, strata = ~strata,
                      weights = ~weight, fpc = ~fpc)

#############################################################################
# Calibrate on the regression coefficients of a *single model* known at the #
# *overall population level*                                                #
#############################################################################
# Suppose you know the coefficients of the following linear model, which you
# obtained fitting the model on some external source (e.g. Census or register
# data):
model0 <- y ~ emp.num*emp.cl

# Here, use the sbs sampling frame available in ReGenesees to simulate the
# external source and compute the values of the regression coefficients:
B0 <- coef(lm(model0, data = sbs.frame))
B0

# Now, prepare the survey design for calibration:
sbsdes0 <- prep.calBeta(design = sbsdes, model = model0, Beta = B0)

# Have a look at the freshly created *synthetic* auxiliary variables: 
head(sbsdes0$variables)

# Then, prepare the control totals dataframe for the calibration task:
pop0 <- pop.calBeta(sbsdes0)

# Have a look...
class(pop0)
pop.desc(pop0)
# ...and note that all the control totals are zero!
pop0

# Lastly, calibrate:
sbscal0 <- e.calibrate(sbsdes0, pop0)

# Check that calibration estimates of regression coefficients now match the
# input B0 values derived from the external source:
svystatB(sbscal0, model0)
B0
# OK


#########################################################################
# Calibrate simultaneously on the regression coefficients of *different #
# models* known at the *overall population level*                       #
#########################################################################
# Suppose you know the coefficients of the following linear models, which you
# obtained fitting the models on some external sources
model1 <- va.imp2 ~ emp.num:emp.cl + nace.macro - 1
model2 <- y ~ dom3 - 1

# Here, use the sbs sampling frame available in ReGenesees to simulate the
# external sources and compute the values of the regression coefficients:
B1 <- coef(lm(model1, data = sbs.frame))
B2 <- coef(lm(model2, data = sbs.frame))

  ## First, just for illustration, calibrate only on B1 regression coefficients:
  sbsdes1 <- prep.calBeta(sbsdes, model = model1, Beta = B1)
  pop1 <- pop.calBeta(sbsdes1)
  sbscal1 <- e.calibrate(sbsdes1, pop1)

  # Check that calibration estimates of regression coefficients now match the
  # input B1 values derived from the external source...
  svystatB(sbscal1, model1)
  B1
  # ...but, of course, do *not* match those of B2:
  svystatB(sbscal1, model2)
  B2
  # OK

  ## Second, just for illustration, calibrate only on B2 regression coefficients:
  sbsdes2 <- prep.calBeta(sbsdes, model = model2, Beta = B2)
  pop2 <- pop.calBeta(sbsdes2)
  sbscal2 <- e.calibrate(sbsdes2, pop2)

  # Check that calibration estimates of regression coefficients now match the
  # input B2 values derived from the external source...
  svystatB(sbscal2, model2)
  B2
  # ...but, of course, do *not* match those of B1:
  svystatB(sbscal2, model1)
  B1
  # OK

## Now, calibrate *simultaneously* on *B1 and B2* regression coefficients
# Prepare the survey design for the joint calibration task:
sbsdes1_2 <- prep.calBeta(sbsdes, model = list(model1, model2), Beta = list(B1, B2))

# Prepare the control totals dataframe for the joint calibration task:
pop1_2 <- pop.calBeta(sbsdes1_2)

# Have a look to the control totals (note the presence of two models and
# Beta vectors):
pop.desc(pop1_2)
pop1_2

# Lastly, run the calibration:
sbscal1_2 <- e.calibrate(sbsdes1_2, pop1_2)

# Check that calibration estimates of regression coefficients now match *both*
# the B1 and B2 values derived from the external sources:
svystatB(sbscal1_2, model1)
B1
svystatB(sbscal1_2, model2)
B2
# OK


###############################################################################
# Calibrate simultaneously on the regression coefficients of a *single model* #
# known for *different subpopulations*                                        #
###############################################################################
# NOTE: In this case, both *global* and *partitioned* calibration tasks are
#       possible, and both will be illustrated below.

# Suppose you know the coefficients of the following linear model, which you
# obtained fitting the model on some external source (e.g. Census or register
# data) *within subpopulations* defined by some factor variable(s):
model <- va.imp2 ~ emp.num:emp.cl + nace.macro - 1

# Here, use the sbs sampling frame available in ReGenesees to simulate the
# external source and suppose the subpopulations are defined by variable 'dom3'.
# Thus, compute the values of the regression coefficients as follows:
B <- by(sbs.frame, sbs.frame$dom3, function(df) coef(lm(model, data = df)))
B

## Let's start with the *global* solution
# Prepare the survey design for the calibration task (note that the 'by'
# argument is used):
sbsdes.g <- prep.calBeta(sbsdes, model, Beta = B, by = ~dom3)

# Prepare the control totals for the calibration task:
pop.g <- pop.calBeta(sbsdes.g)

# Have a look to the control totals (note the presence of one Beta vector *for
# each domain*):
pop.desc(pop.g)
pop.g

# Run the calibration:
sbscal.g <- e.calibrate(sbsdes.g, pop.g)

# Check that calibration estimates of regression coefficients now match the
# input B values derived from the external source *for each domain*:
svystatB(sbscal.g, model, ~dom3)
B
# OK

## Let's proceed with the *partitioned* solution
# Prepare the survey design for the calibration task (note that 'by' and
# 'partition' arguments are used):
sbsdes.p <- prep.calBeta(sbsdes, model, Beta = B, by = ~dom3, partition = TRUE)

# Prepare the control totals for the calibration task:
pop.p <- pop.calBeta(sbsdes.p)

# Have a look to the control totals (note the presence of one Beta vector *for
# each domain*):
pop.desc(pop.p)
pop.p

# Run the calibration:
sbscal.p <- e.calibrate(sbsdes.p, pop.p)

# Check that calibration estimates of regression coefficients now match the
# input B values derived from the external source *for each domain*:
svystatB(sbscal.p, model, ~dom3)
B
# OK

## Lastly, check that calibration weights obtained using the *global* and 
## *partitioned* solution are the same:
g.range(sbscal.g)
g.range(sbscal.p)
all.equal(weights(sbscal.g), weights(sbscal.p))
# OK


###########################################################
# BONUS TIP: Calibration on the mean (or on domain means) #
#            of one variable or multiple variables.       #
###########################################################
# Since the domain mean of a numeric variable can be thought as a
# regression coefficient (see the 'Examples' section of ?svystatB),
# you can use the ReGenesees facilities documented above to
# *calibrate on the mean (or on domain means)* of *one variable
# or multiple variables*.
# NOTE: The examples below cover the following cases of calibration on:
#       (i)   The overall mean of a single variable.
#       (ii)  The means of a single variable within domains of just
#             one type.
#       (iii) The domain means of several variables, with multiple
#             and different domain types.

# Load artificial household survey data and define a survey design:
data(data.examples)
exdes <- e.svydesign(data = example, ids = ~towcod + famcod,
                     strata = ~SUPERSTRATUM, weights = ~weight)
exdes <- des.addvars(exdes, ones = 1)

## CASE (i): Calibrate on the overall mean of a single variable
# Suppose you know with satisfactory accuracy the *average income* of your
# target population (but you do *not* have reliable information on the
# *total income*, nor on the *total number of individuals*):
income.AVG <- 1270

# You can calibrate on *average income* as follows:
exdes.new <- prep.calBeta(exdes, income ~ 1, Beta = income.AVG)
pop.new <- pop.calBeta(exdes.new)
excal.new <- e.calibrate(exdes.new, pop.new)

# Now, check that calibration estimate of average income now match the known
# value derived from the external source *without residual uncertainty*:
svystatTM(excal.new, ~income, estimator = "Mean")
income.AVG

# ...while there is *residual uncertainty* in the estimates of the numerator and
# denominator totals:
svystatTM(excal.new, ~income + ones, estimator = "Total")
# OK

## CASE (ii): Calibrate on the means of a single variable within domains
#             of just one kind
# You can calibrate on *domain means* along the lines illustrated above (note,
# however, that argument 'by' would provide an alternative way to achieve
# the same result).
# Suppose you know with satisfactory accuracy the *average income* by the
# crossclassification of sex and marital status:
income.AVG.sex.marstat <- c(f.married   = 1310, m.married   = 1260,
                            f.unmarried = 1150, m.unmarried = 1200,
                            f.widowed   = 1380, m.widowed   = 1300)

# Run the calibration on *average income by sex:marstat* as follows:
exdes.new <- prep.calBeta(exdes, income ~ sex:marstat -1,
                          Beta = income.AVG.sex.marstat)
pop.new <- pop.calBeta(exdes.new)
excal.new <- e.calibrate(exdes.new, pop.new)

# Now, check that calibration estimates of average income by domains now match
# the known values derived from the external source *without residual
# uncertainty*:
svystatTM(excal.new, ~income, ~sex:marstat, estimator = "Mean")
income.AVG.sex.marstat

# ...while there is *residual uncertainty* in the estimates of the numerator and
# denominator totals:
svystatTM(excal.new, ~income + ones, ~sex:marstat, estimator = "Total")
# OK

## CASE (iii): Calibrate on the domain means of several variables, with multiple
#              and different domain types.
# Suppose you know with satisfactory accuracy:
# - the average income by sex:
# - the average income by marstat:
# - the average of variable z by age (variable 'age5c', 5 classes):
income.AVG.sex <- c("f" = 1245, "m" = 1250)
income.AVG.marstat <- c("married" = 1260, "unmarried" = 1230, "widowed" = 1290)
z.AVG.age5c <- c("1" = 125, "2" = 130, "3" = 135, "4" = 125, "5" = 140)

# Run the calibration as follows:
exdes.new <- prep.calBeta(exdes, model = list(income ~ sex -1,
                                              income ~ marstat -1,
                                              z ~ age5c -1),
                                 Beta = list(income.AVG.sex,
                                             income.AVG.marstat,
                                             z.AVG.age5c)
                          )
pop.new <- pop.calBeta(exdes.new)
excal.new <- e.calibrate(exdes.new, pop.new)

# Now, check that calibration estimates match the known domain means derived
# from the external source:
svystatTM(excal.new, ~income, ~sex, estimator = "Mean")
income.AVG.sex
svystatTM(excal.new, ~income, ~marstat, estimator = "Mean")
income.AVG.marstat
svystatTM(excal.new, ~z, ~age5c, estimator = "Mean")
z.AVG.age5c
# OK

}
\keyword{survey}