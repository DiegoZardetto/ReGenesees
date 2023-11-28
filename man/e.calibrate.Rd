\name{e.calibrate}
\alias{e.calibrate}
\alias{cal.analytic}
\alias{ecal.status}

\title{Calibration of Survey Weights}

\description{
Adds to an \code{analytic} object the calibrated weights column. 
}

\usage{
e.calibrate(design, df.population,
            calmodel = if (inherits(df.population, "pop.totals"))
                       attr(df.population, "calmodel"),
            partition = if (inherits(df.population, "pop.totals"))
                        attr(df.population, "partition") else FALSE,
            calfun = c("linear", "raking", "logit"),
            bounds = c(-Inf, Inf), aggregate.stage = NULL,
            sigma2 = NULL, maxit = 50, epsilon = 1e-07, force = TRUE)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{df.population}{Data frame containing the known population totals for the auxiliary variables.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model (see \sQuote{Details}); \code{FALSE} (the default) implies no calibration domains.}
  \item{calfun}{\code{character} specifying the distance function for the calibration process; the default is
                \code{'linear'}.}  
  \item{bounds}{Allowed range for the ratios between calibrated and initial weights; the default is
                \code{c(-Inf,Inf)}.}  
  \item{aggregate.stage}{An integer: if specified, causes the calibrated weights to be constant within sampling
                         units at this stage.}
  \item{sigma2}{Formula specifying a possible heteroskedasticity effect in the calibration model;
                \code{NULL} (the default) implies homoskedasticity.}
  \item{maxit}{Maximum number of iterations for the Newton-Raphson algorithm; the default is \code{50}.}
  \item{epsilon}{Tolerance for the absolute relative differences between the population totals and the corresponding
                 estimates based on the calibrated weights; the default is \code{10^-7}.}
  \item{force}{If \code{TRUE}, whenever the calibration algorithm does not converge, forces the function to return
               a value (see \sQuote{Details} and \sQuote{Calibration process diagnostics}); the default is \code{TRUE}.}
}

\details{
This function creates an object of class \code{cal.analytic}. A \code{cal.analytic} object makes it possible to compute estimates and standard errors of calibration estimators [Deville, Sarndal 92] [Deville, Sarndal, Sautory 93].

The mandatory argument \code{calmodel} symbolically defines the calibration model you intend to use, that is - in the language of the Generalized Regression Estimator - the assisting linear regression model underlying the calibration problem. More specifically, the \code{calmodel} formula identifies the auxiliary variables and the constraints for the calibration problem, with a notation inspired by [Wilkinson, Rogers 73]. For example, \code{calmodel=~(X+Z):C+(A+B):D-1} defines the calibration problem in which constraints are imposed: (i) on the totals of auxiliary (quantitative) variables \code{X} and \code{Z} within the subpopulations identified by the (qualitative) classification variable \code{C} and, at the same time, (ii) on the absolute frequency of the (qualitative) variables \code{A} and \code{B} within the subpopulations identified by the (qualitative) classification variable \code{D}.\cr
The \code{design} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}).

Problems for which one or more qualitative variables can be \emph{\dQuote{factorized}} in the formula that specifies the calibration model, are particularly interesting. These variables split the population into non-overlapping subpopulations known as \emph{\dQuote{calibration domains}} for the model. An example is provided by the statement \code{calmodel=~(A+B+X+Z):D-1} in which the variable that identifies the calibration domains is \code{D}; similarly, the formula \code{calmodel=~(A+B+X+Z):D1:D2-1} identifies as calibration domains the subpopulations determined by crossing the modalities of \code{D1} and \code{D2}. The interest in models of this kind lies in the fact that the \emph{global} calibration problem they describe can, actually, be broken down into \emph{local} subproblems, one per calibration domain, which can be solved separately [Vanderhoeft 01]. Thus, for example, the global problem defined by \code{calmodel=~(A+B+X+Z):D-1} is equivalent to the sequence of problems defined by the \emph{\dQuote{reduced model}} \code{calmodel=~A+B+X+Z-1} in each of the domains identified by the modalities of \code{D}. The opportunity to separately solve the subproblems related to different calibration domains achieves a significant reduction in computation complexity: the gain increases with increasing survey data size and (most importantly) with increasing auxiliary variables number.

The optional argument \code{partition} makes it possible to choose, in cases in which the calibration problem can be factorized, whether to solve the problem globally or in a partitioned way (that is, separately for each calibration domain). The global solution (which is the default option) can be selected invoking the \code{e.calibrate} function with \code{partition=FALSE}. To request the partitioned solution - a strongly recommended option when dealing with a lot of auxiliary variables and big data sizes - it is necessary to specify via \code{partition} the variables defining the calibration domains for the model. If a formula is passed through the \code{partition} argument (for example: \code{partition=~D1:D2}), the program checks that \code{calmodel} actually describes a "reduced model" (for example: \code{calmodel=~A+B+X+Z-1}), that is it does \emph{not} reference any of the partition variables; if this is not the case, the program stops and prints an error message. Notice that a formula like \code{partition=~D1+D2} will be automatically translated into the factor-crossing formula \code{partition=~D1:D2}.\cr
The \code{design} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).

The mandatory argument \code{df.population} is used to specify the known totals of the auxiliary variables referenced by \code{calmodel} within the subpopulations (if any) identified by \code{partition}. These known totals must be stored in a data frame whose structure (i) depends on the values of \code{calmodel} and \code{partition} and (ii) must conform to a standard. In order to facilitate understanding of and compliance with this standard, the \pkg{ReGenesees} package provides the user with four functions: \code{\link{pop.template}}, \code{\link{population.check}}, \code{\link{pop.desc}} and \code{\link{fill.template}}. The \code{pop.template} function is able to guide the user in constructing the known totals data frame for a specific calibration problem, the \code{\link{pop.desc}} function provides a natural language description of the template structure, the \code{\link{fill.template}} function can be exploited to automatically fill the template when a sampling frame is available, while the \code{population.check} function allows to check whether a known totals data frame conforms to the standard required by \code{e.calibrate}. In any case, if the \code{df.population} data frame does not comply with the standard, the \code{e.calibrate} function stops and prints an error message: the meaning of the message should help the user diagnose the cause of the problem.

The \code{calfun} argument identifies the distance function to be used in the calibration process. Three built-in functions are provided: \code{"linear"}, \code{"raking"}, and \code{"logit"} (see [Deville, Sarndal, Sautory 93]). The default is \code{"linear"}, which corresponds to the euclidean metric and yields the Generalized Regression Estimator (provided that no range restrictions are imposed on the g-weights). The \code{"raking"} distance corresponds to the \emph{\dQuote{multiplicative method}} of [Deville, Sarndal, Sautory 93]. 

The \code{bounds} argument allows to add \emph{\dQuote{range constraints}} to the calibration problem. To be precise, the interval defined by \code{bounds} will contain the values of the ratios between final (calibrated) and initial (direct) weights. The default value is \code{c(-Inf,Inf)}, i.e. no range constraints are imposed. These constraints are optional unless the \code{"logit"} function is selected: in the latter case the range defined by \code{bounds} has to be finite (see, again, [Deville, Sarndal, Sautory 93]).       

The value passed by the \code{aggregate.stage} argument must be an integer between \code{1} and the number of sampling stages of \code{design}. If specified, causes the calibrated weights to be constant within sampling units selected at the \code{aggregate.stage} stage (actually this is only allowed if the initial weights had already this property, as it is sometimes the case in multistage cluster sampling). If not specified, the calibrated weights may differ even for sampling units with identical initial weights. The same holds if some final units belonging to the same cluster selected at the stage \code{aggregate.stage} fall in distinct calibration domains (i.e. if the domains defined by \code{partition} "cut across" the \code{aggregate.stage}-stage clusters).

The argument \code{sigma2} can be used to take into account a possible heteroskedasticity effect in the (assisting linear regression model underlying the) calibration problem. In such cases, \code{sigma2} must identify some variable to which the variances of the error terms are believed to be proportional. Notice that \code{sigma2} can also be interpreted from a "purely calibration-based" point of view: it corresponds to the \eqn{1/q_k} unit-weights appearing inside the distance measures of [Deville, Sarndal 92] [Deville, Sarndal, Sautory 93]. The final effect is, on average, that \emph{calibrated weights associated to higher values of} \code{sigma2} \emph{tend to stay closer to their corresponding initial weights}.\cr
Note that it is technically possible to exploit this behaviour in order to \emph{prevent some subset of the initial weights from being altered by calibration}. The trick is simple: just build a convenience \code{sigma2} variable whose values are set to some very high value (e.g. \code{1E12}) for those units whose initial weight must be preserved, and to \code{1} otherwise (see \sQuote{Examples}). Nevertheless, this trick should be used sparingly and very carefully, as otherwise it may: (i) cause the calibration algorithm to not converge, (ii) result in introducing bias in calibration estimates. In particular, with respect to bias, one should not select the units whose weight must be preserved on the basis of the current sample.\cr
The \code{sigma2} formula can reference just a single \code{design} variable: such variable must be \code{numeric}, strictly positive and must not contain \code{NA}s. If \code{aggregate.stage} is specified, \code{sigma2} must obviously be constant inside \code{aggregate.stage}-stage clusters (otherwise the function stops and prints an error message).

The \code{maxit} argument sets the maximum number of iteration for the Newton-Raphson algorithm that is used to solve the calibration problem (the only exception being \emph{unbounded  linear} calibration, i.e. \code{calfun='linear'} and \code{bounds=c(-Inf, Inf)}, which is actually handled by directly solving a linear problem). The default value of \code{maxit} is \code{50}.

The \code{epsilon} argument determines the convergence criterion for the optimization algorithm: it fixes the maximum allowed absolute value for the relative differences between the population totals and the corresponding estimates based on the calibrated weights. The default value is \code{10^-7}.

The calibrated weights computed by \code{e.calibrate} must ensure that the calibration estimators of the auxiliary variables \emph{exactly} match the corresponding known population totals. It is, however, possible (more likely when range constraints are imposed) that, for a specific calibration problem and for given values of \code{epsilon} and \code{maxit}, the solving algorithm does not converge. In this case, if \code{force = FALSE}, \code{e.calibrate} stops and prints an error message. If - on the contrary - \code{force = TRUE}, the function is forced to return the best approximation achieved for the calibrated weights, nevertheless signaling the calibration failure by a warning (see also Section 'Calibration process diagnostics').
}

\section{Calibration Process Diagnostics}{
When, dealing with a factorizable calibration problem, the user selects the partitioned solution, the global calibration problem gets split into as many \emph{sub-problems} as the number of subpopulations defined by \code{partition}. In turn, each one of these calibration sub-problems can end without convergence on any one of the involved auxiliary variables. A calibration process with such a complex structure needs some ad hoc tool for error diagnostics. For this purpose, every call to \code{e.calibrate} creates, by side effect, a dedicated data structure named \code{ecal.status} into the \code{.GlobalEnv}.

\code{ecal.status} is a list with up to three components: the first, \code{"call"}, identifies the call to \code{e.calibrate} that generated the list, the second, \code{return.code}, is a matrix each element of which identifies the return code of a specific calibration sub-problem. The meaning of the return codes is as follows:
\preformatted{
CODE        MEANING
  -1........not yet tackled sub-problem;
   0........solved sub-problem (convergence achieved);
   1........unsolved sub-problem (no convergence): output forced.
}

Recall that the latter return code (\code{1}) may only occur if \code{force = TRUE}.\cr
If any \code{return.code} equal to \code{1} exists, the \code{ecal.status} list gains a third component named \code{"fail.diagnostics"} which is itself a list; its components correspond to sub-problems for which convergence was not achieved, and store useful information about the auxiliary variables for which calibration constraints are violated. Therefore, users can exploit \code{ecal.status} to identify sub-problems and variables from which errors stemmed, hence taking a step forward to eliminate them.

Notice, lastly, that the \code{ecal.status} list will also be persistently bound to the \code{e.calibrate} return object, stored inside a dedicated attribute. For the inspection of such diagnostics information the \code{\link{check.cal}} function is available.   
}

\value{
An object of class \code{cal.analytic}. The data frame it contains includes (in addition to the data already stored in \code{design}) the calibrated weights column. The name of this column is obtained by pasting the name of the initial weights column with the string \code{".cal"}.
}

\note{
The \code{cal.analytic} class is a specialization of the \code{\link{analytic}} class; this means that an object created by \code{e.calibrate} inherits from the \code{\link{analytic}} class and you can use on it all methods defined on the latter class, e.g. \code{print}, \code{summary}, \code{\link{weights}}. Moreover, a calibrated design can be passed again to \code{e.calibrate}, thus undergoing further calibration steps.
}

\author{Diego Zardetto}

\references{
Deville, J.C., Sarndal, C.E. (1992) \emph{\dQuote{Calibration Estimators in Survey Sampling}}, Journal of the American Statistical Association, Vol. 87, No. 418, pp. 376-382.

Deville, J.C., Sarndal, C.E., Sautory, O. (1993) \emph{\dQuote{Generalized Raking Procedures in Survey Sampling}}, Journal of the American Statistical Association, Vol. 88, No. 423, pp. 1013-1020.

Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.

Wilkinson, G.N., Rogers, C.E. (1973) \emph{\dQuote{Symbolic Description of Factorial Models for Analysis of Variance}}, Journal of the Royal Statistical Society, series C (Applied Statistics), Vol. 22, pp. 181-191.

Vanderhoeft, C. (2001) \emph{\dQuote{Generalized Calibration at Statistic Belgium}}, Statistics Belgium Working Paper n. 3.

Sarndal, C.E., Lundstrom, S. (2005) Estimation in surveys with nonresponse. John Wiley & Sons.

Scannapieco, M., Zardetto, D., Barcaroli, G. (2007) \emph{\dQuote{La Calibrazione dei Dati con R: una Sperimentazione sull'Indagine Forze di Lavoro ed un Confronto con GENESEES/SAS}}, Contributi Istat n. 4., \url{https://www.istat.it/it/files/2018/07/2007_4.pdf}.
}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata, \code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}} and \code{\link{svystatL}} for calculating estimates and standard errors, \code{\link{pop.template}} for constructing known totals data frames in compliance with the standard required by \code{e.calibrate}, \code{\link{population.check}} to check that the known totals data frame satisfies that standard, \code{\link{pop.desc}} to provide a natural language description of the template structure, \code{\link{fill.template}} to automatically fill the template when a sampling frame is available, \code{\link{bounds.hint}} to obtain a hint for range restricted calibration, \code{\link{g.range}} to asses the variation of weights after calibration and \code{\link{check.cal}} to check if calibration constraints have been fulfilled.
}

\examples{
######################################################################
# Calibration of a design object according to different calibration  #
# models (the known totals data frames pop01, \ldots, pop05p and the    #
# bounds vector are all contained in the data.examples file).        #
# For the examples relating to calibration models that can be        #
# factorized both a global and a partitioned solution are given.     #
######################################################################

# Load household data:
data(data.examples)

# Creation of the object to be calibrated:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)


# 1) Calibration on the total number of units in the population
#   (totals in pop01):
descal01<-e.calibrate(design=des,df.population=pop01,calmodel=~1,
          calfun="logit",bounds=bounds,aggregate.stage=2)

# Printing descal01 immediately recalls that it is a
# "calibrated" object:
descal01

# Use the summary() function if you need some additional details, e.g.:
summary(descal01)

# Use the 'variables' slot to extract survey data, e.g.:
head(descal01$variables)

# Use the weights() function to extract weights, e.g.:
summary(weights(descal01))

# Checking the result (first add the new 'ones' variable 
# to estimate the number of final units in the population):
descal01<-des.addvars(descal01,ones=1)
svystatTM(descal01, ~ones)


# 2) Calibration on the marginal distributions of sex and marstat
#    (totals in pop02):
descal02<-e.calibrate(design=des,df.population=pop02,
          calmodel=~sex+marstat-1,calfun="logit",bounds=bounds,
          aggregate.stage=2)

# Checking the result:
svystatTM(descal02,~sex+marstat)


# 3) Calibration (global solution) on the joint distribution of sex
#    and marstat (totals in pop03):
descal03<-e.calibrate(design=des,df.population=pop03,
          calmodel=~marstat:sex-1,calfun="logit",bounds=bounds)

# Checking the result:
svystatTM(descal03,~sex,~marstat) # or: svystatTM(descal03,~marstat,~sex)

# which, obviously, is not respected by descal02 (notice the size of SE):
svystatTM(descal02,~sex,~marstat)


# 3.1) Again a calibration on the joint distribution of sex and marstat
#      but, this time, with the partitioned solution (partition=~sex,
#      totals in pop03p):
descal03p<-e.calibrate(design=des,df.population=pop03p,
           calmodel=~marstat-1,partition=~sex,calfun="logit",
           bounds=bounds)

# Checking the result:
svystatTM(descal03p,~sex,~marstat)


# 4) Calibration (global solution) on the totals for the quantitative
#    variables x1, x2 and x3 in the subpopulations defined by the
#    regcod variable (totals in pop04):
descal04<-e.calibrate(design=des,df.population=pop04,
          calmodel=~(x1+x2+x3):regcod-1,calfun="logit",
          bounds=bounds,aggregate.stage=2)

# Checking the result:
svystatTM(descal04,~x1+x2+x3,~regcod)


# 4.1) Same problem with the partitioned solution (partition=~regcod,
#      totals in pop04p):
descal04p<-e.calibrate(design=des,df.population=pop04p,
           calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
           bounds=bounds,aggregate.stage=2)

# Checking the result:
svystatTM(descal04p,~x1+x2+x3,~regcod)


# 5) Calibration (global solution) on the total for the quantitative
#    variable x1 and on the marginal distribution of the qualitative
#    variable age5c, in the subpopulations defined by crossing sex
#    and marstat (totals in pop05):
descal05<-e.calibrate(design=des,df.population=pop05,
          calmodel=~(age5c+x1):sex:marstat-1,calfun="logit",
          bounds=bounds)

# Checking the result:
svystatTM(descal05,~age5c+x1,~sex:marstat)


# 5.1) Same problem with the partitioned solution (partition=~sex:marstat,
#      totals in pop05p):
descal05p<-e.calibrate(design=des,df.population=pop05p,
           calmodel=~age5c+x1-1,partition=~sex:marstat,
           calfun="logit",bounds=bounds)

# Checking the result:
svystatTM(descal05p,~age5c+x1,~sex:marstat)

# Notice that 3.1 and 5.1) 5.2) do not impose the aggregate.stage=2
# condition. This condition cannot, in fact, be fulfilled because
# in both cases the domains defined by partition "cut across"
# the des second stage clusters (households). To compare the results,
# the same choice was also made for 3) and 5).


# 5.2) Just a single example to inspect the ecal.status list generated
#      for diagnostics purposes.
#      Let's shrink the bounds in order to prevent perfect convergence
#      (recall that force=TRUE by default):
approx.cal<-e.calibrate(design=des,df.population=pop05p,
            calmodel=~age5c+x1-1,partition=~sex:marstat,
            calfun="logit",bounds=c(0.95,1.05))

# ...now use check.cal function to assess the amount of calibration
# constraints violation:
check.cal(approx.cal)

# ...or (equivalently) inspect directly ecal.status:
ecal.status


#################################################
# Some examples illustrating how calibration    #
# can be exploited to reduce nonresponse bias   #
# (see, e.g. [Sarndal, Lundstrom 05]).          #
#################################################

# Load sbs data:
data(sbs)

  #######################
  # Full-response case. #
  #######################

  # Create a full-response design object:
  sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

  # Now estimate the average value added and its 95% confidence interval:
  mean.VA<-svystatTM(design=sbsdes,y=~va.imp2,estimator="Mean",vartype= "cvpct",
           conf.int=TRUE,conf.lev=0.95)
  mean.VA

  # Compare the obtained estimate with the true population parameter:
  MEAN.VA<-mean(sbs.frame$va.imp2)
  MEAN.VA

  # We get a small overestimation of about 4%...
  100*(coef(mean.VA)-MEAN.VA)/MEAN.VA
  
  # which, anyway, doesn't indicate a significant bias for the
  # full-response sample, because the 95% confidence interval
  # covers the true value.

  
  ##################################################
  # Nonresponse case: assume a response propensity #
  # which increases with enterprise size.          #
  ##################################################

  # Set bigger response probabilities for bigger firms,
  # e.g. exploiting available information about the
  # number of employees (emp.cl):
  levels(sbs$emp.cl)
  p.resp <- c(.4,.6,.8,.95,.99)

  # Tie response probabilities to sample observations:
  pr<-p.resp[unclass(sbs$emp.cl)]

  # Now, randomly select a subsample of responding units from sbs:
  set.seed(12345)          # (fix the RNG seed for reproducibility)
  rand<-runif(1:nrow(sbs))
  sbs.nr<-sbs[rand<pr,]

  # This implies an overall response rate of about 73%:
  nrow(sbs.nr)/nrow(sbs) 
 
  # Treat the non-response sample as it was complete: this should
  # lead to biased estimates of value added, as the latter is
  # positively correlated with firms size...
  sbsdes.nr<-e.svydesign(data=sbs.nr,ids=~id,strata=~strata,weights=~weight)
  
  #...indeed:
  old.op <- options("RG.lonely.psu"="adjust")   # (prevent lonely-PSUs troubles)
  mean.VA.nr<-svystatTM(design=sbsdes.nr,y=~va.imp2,estimator="Mean",
              vartype= "cvpct",conf.int=TRUE,conf.lev=0.95)
  mean.VA.nr

  # and, comparing with the true population average, we see a
  # significant overestimation effect, with the 95% confidence
  # interval not even covering the parameter:
  MEAN.VA
  
  # Nonresponse bias can be effectively reduced by calibrating
  # on variables explaining the response propensity: e.g., in
  # the present example, on the population distribution of emp.cl:
    # Prepare the known totals template...
    N.emp.cl<-pop.template(data=sbs.nr,calmodel=~emp.cl-1)
    N.emp.cl

    # Fill it by using the sampling frame...
    N.emp.cl<-fill.template(sbs.frame,N.emp.cl)
    N.emp.cl

    # Lastly calibrate:
      # Get a hint on the calibration bounds:
      hint<-bounds.hint(sbsdes.nr,N.emp.cl)
      sbscal.nr<-e.calibrate(design=sbsdes.nr,df.population=N.emp.cl,
                 bounds=hint)
      sbscal.nr

  # Now estimate the average value added on the calibrated design:
  mean.VA.cal.nr<-svystatTM(design=sbscal.nr,y=~va.imp2,estimator="Mean",
              vartype= "cvpct",conf.int=TRUE,conf.lev=0.95)

  # options(old.op)   # (reset variance estimation options)

  # As expected, we see a significant bias reduction:
  MEAN.VA
  mean.VA.nr
  mean.VA.cal.nr

  # Even if the 95% confidence interval still doesn't cover the
  # true value, by calibration we passed from an initial overestimation
  # of about 33% to a 7% one: 
  100*(coef(mean.VA.nr)-MEAN.VA)/MEAN.VA
  100*(coef(mean.VA.cal.nr)-MEAN.VA)/MEAN.VA


  #################################################
  # A multi-step calibration example showing that #
  # a calibrated object can be calibrated again   #
  # (this can be sometimes useful in practice):   #
  # Step 1: calibrate to reduce nonresponse bias; #
  # Step 2: calibrate again to gain efficiency.   #
  #################################################

  # Suppose you already performed a first calibration step,
  # as shown in the example above, with the aim of softening
  # nonresponse bias:
  sbscal.nr

  # Now you may want to calibrate again in order to reduce
  # estimators variance, by using further available auxiliary
  # information, e.g. the total number of employees (emp.num)
  # and enterprises (ent) inside the domains obtained
  # by crossing nace.macro and region:

    # Build the second step population totals template:
    pop2<-pop.template(sbscal.nr,
          calmodel=~emp.num+ent-1,
          partition=~nace.macro:region)

    # Use the fill.template function to (i) automatically compute
    # the totals from the universe (sbs.frame) and (ii) safely fill
    # the template:
    pop2<-fill.template(universe=sbs.frame,template=pop2)

    # Now perform the second calibration step:
        # Get a hint on the calibration bounds:
        hint2<-bounds.hint(sbscal.nr,pop2)
        sbscal.nr2<-e.calibrate(design=sbscal.nr,df.population=pop2,
                    bounds=hint2)

  # Notice that printing sbscal.nr2 you immediately understand
  # that it is a "twice-calibrated" object: 
  sbscal.nr2

  # Notice also that, even if the second calibration step causes
  # sbscal.nr2 to be no more exactly calibrated with respect to 
  # emp.cl (look at the cvpct values)...
  old.op <- options("RG.lonely.psu"="adjust")   # (prevent lonely-PSUs troubles)
  svystatTM(design=sbscal.nr2,y=~emp.cl,vartype="cvpct")

  # ...the nonresponse bias has not been resurrected (i.e. it gets stuck
  # to its previous 7%):
  mean.VA.cal.nr2<-svystatTM(design=sbscal.nr2,y=~va.imp2,estimator="Mean",
                   vartype= "cvpct",conf.int=TRUE,conf.lev=0.95)

  options(old.op)   # (reset variance estimation options)

  mean.VA.cal.nr2
  100*(coef(mean.VA.cal.nr2)-MEAN.VA)/MEAN.VA


  ##############################################################
  # Provided the auxiliary variables are chosen in a smart way #
  # a single calibration step can simultaneously succeed in:   #
  # (i)  softening nonresponse bias;                           #
  # (ii) reducing estimators variance.                         #
  ##############################################################

  # Let's come back to the original design with nonresponse:
  sbsdes.nr

  # Now, let's try to calibrate simultaneously on (see examples above):
  # (i)  the population distribution of emp.cl;
  # (ii) the total number of employees (emp.num) and enterprises (ent)
  #      inside the domains obtained by crossing nace.macro and region:

    # Build the population totals template (notice that we are now forced
    # to a global calibration, as we are assuming to ignore emp.cl counts
    # inside domains obtained by crossing nace.macro and region):
    pop1<-pop.template(sbs.nr,
          calmodel=~emp.cl+(emp.num+ent):nace.macro:region-1)

    # Use the fill.template function to (i) automatically compute
    # the totals from the universe (sbs.frame) and (ii) safely fill
    # the template:
    pop1<-fill.template(universe=sbs.frame,template=pop1)

    # Now perform the single calibration step:
        # Get a hint on the calibration bounds:
        hint1<-bounds.hint(sbsdes.nr,pop1)
        sbscal.nr1<-e.calibrate(design=sbsdes.nr,df.population=pop1,
                    bounds=hint1)

    sbscal.nr1

  # Now:
  # (i) verify the nonresponse bias reduction effect:
    old.op <- options("RG.lonely.psu"="adjust")  #(prevent lonely-PSUs troubles)
    mean.VA.cal.nr1<-svystatTM(design=sbscal.nr1,y=~va.imp2,estimator="Mean",
                     vartype= "cvpct",conf.int=TRUE,conf.lev=0.95)
    options(old.op)

    mean.VA.cal.nr1
    100*(coef(mean.VA.cal.nr1)-MEAN.VA)/MEAN.VA

  # thus we are back to ~7%, as for the previous 2-step calibration example.

  # (ii) compare cvpct with the previous 2-step calibration example:
  mean.VA.cal.nr1
  mean.VA.cal.nr2

  # hence, both bias reduction and efficiency are almost the same in 2-step and
  # single step calibration (auxiliary information being equal): the choice
  # will often depend on practical considerations (e.g. convergence, computation
  # time).  


############################################################################
# Example with heteroskedastic assisting linear model: shows how to obtain #
# the ratio estimator of a total by calibration.                           #
############################################################################

# Load sbs data:
data(sbs)

# Create the design object to be calibrated:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

# Suppose you have to calibrate on the total amount of employees:
# Prepare the template:
pop<-pop.template(data=sbsdes,calmodel=~emp.num-1)
pop

# Fill it by using the sampling frame (sbs.frame)...
pop<-fill.template(sbs.frame,pop)
pop

# ... thus the total number of employees is 984394.
# Now calibrate assuming that error terms variances are proportional
# to emp.num:
sbscal<-e.calibrate(design=sbsdes,df.population=pop,sigma2=~emp.num)

# Now compute the calibration estimator of the total
# of value added (i.e. variable va.imp2)...
VA.tot.cal<-svystatTM(design=sbscal,y=~va.imp2)
VA.tot.cal

#... and observe that this is identical to the ratio estimator of the total...
TOT.emp.num <- pop[1, 1]
VA.ratio<-svystatL(design=sbsdes, expression(TOT.emp.num * (va.imp2/emp.num)))
VA.ratio

# ...as it must be.

# Recall that, for the calibration problem above, one must expect, by virtue of
# simple theoretical arguments, that the g-weights are constant and equal to the
# ratio between the known total of emp.num (984394) and its HT estimate.
# This property is exactly satisfied by our numerical results, see below:
pop[1, 1]/coef(svystatTM(sbsdes, ~emp.num))
g.range(sbscal)

# ...as it must be.


#########################################################################
# A second example of calibration with heteroskedastic assisting linear #
# model. Shows that calibrated weights associated to higher values of   #
# sigma2 tend to stay closer to their corresponding initial weights.    #
#########################################################################

# Perform a calibration process which exploits as auxiliary
# information the total number of employees (emp.num)
# and enterprises (ent) inside the domains obtained by:
#  i) crossing nace2 and region;
# ii) crossing emp.cl, region and nace.macro;

# Build the population totals template:
pop<-pop.template(sbsdes,
     calmodel=~(emp.num+ent):(nace2+emp.cl:nace.macro)-1,
     partition=~region)

# Use the fill.template function to (i) automatically compute
# the totals from the universe (sbs.frame) and (ii) safely fill
# the template:
pop<-fill.template(universe=sbs.frame,template=pop)

# Now calibrate:
  # 1) First, without any heteroskedasticy effect
  sbscal1<-e.calibrate(sbsdes,pop,calfun="linear",bounds = c(0.01, 3),
           sigma2=NULL)
  
  # 2) Then, with heteroskedastic effect proportional to emp.num:
  sbscal2<-e.calibrate(sbsdes,pop,calfun="linear",bounds = c(0.01, 3),
           sigma2=~emp.num)

# Compute the g-weights for both the calibrated objects:
g1<-weights(sbscal1)/weights(sbsdes)
g2<-weights(sbscal2)/weights(sbsdes)

# Now visually compare the absolute deviations from 1 of the g-weights
# as a function of emp.num:
plot(log10(sbs$emp.num),abs(g1-1), col="blue", pch=19, cex=0.5)
points(log10(sbs$emp.num),abs(g2-1), col="red", pch=19, cex=0.5)

#...as emp.num grows red points clearly tend to stay closer to
# the horizontal axis than blue ones, as expected.


#########################################################################
# A third example. Shows how to exploit the sigma2 argument to prevent  #
# some initial weights from being altered by calibration.               #
#########################################################################

# Let's refer again to object sbsdes:
sbsdes

# Let's assume that for some reason we want to prevent the *highest* initial
# weights from being altered by calibration:
dmax <- max(weights(sbsdes))
dmax

# The relevent units are the following 4:
to.keep <- which(weights(sbsdes) == dmax)
to.keep

# Now, let's prepare a convenience variable (to be later bound to the 'sigma2'
# argument of e.calibrate) whose values are set to a very high value (say 1E12)
# for those units whose initial weight must be preserved, and to 1 otherwise.
# For definiteness, let's call such variable 'fixed':
sbsdes <- des.addvars(sbsdes, fixed = ifelse(weights(sbsdes) == dmax, 1E12, 1))

# Now, let's perform a calibration process which exploits as auxiliary
# information the total number of employees (emp.num)
# and enterprises (ent) inside the domains obtained by:
#  i) crossing region and emp.cl;
# ii) crossing region and nace.macro;

# Build the population totals template:
pop<-pop.template(sbsdes,
     calmodel = ~(emp.num + ent):(emp.cl + nace.macro) - 1,
     partition = ~region)

# Use the fill.template function to (i) automatically compute
# the totals from the universe (sbs.frame) and (ii) safely fill
# the template:
pop<-fill.template(universe=sbs.frame,template=pop)

# Now calibrate:
  # 1) First, *without* any heteroskedasticy effect:
       sbscal1 <- e.calibrate(sbsdes, pop, calfun = "linear", sigma2 = NULL)
       g.range(sbscal1)

       ## As expected, calibration weights of the 4 units to.keep *differ* from
       ## the corresponding initial weights:
       weights(sbsdes)[to.keep]
       weights(sbscal1)[to.keep]

  # 2) Then, *with* heteroskedasticy effect given by our convenience variable
  #    'fixed':
       sbscal2 <- e.calibrate(sbsdes, pop, calfun = "linear", sigma2 = ~fixed)
       g.range(sbscal2)

       ## Let's verify that calibration weights of the 4 units to.keep are now
       ## *equal* to the corresponding initial weights:
       weights(sbsdes)[to.keep]
       weights(sbscal2)[to.keep]

       ## ...as it must be.

# NOTE: It should be clear that the additional request to hold some weights
#       fixed while calibrating will - all other things being equal - increase
#       the probability of non-convergence of the calibration algorithm.


##############################################################
# Calibrating simultaneously on unit-level and cluster-level #
# auxiliary information: an household survey example.        #
##############################################################

# Load household data:
data(data.examples)

# Define the survey design: 
exdes<-e.svydesign(data=example,ids=~towcod+famcod,strata=~stratum,
       weights=~weight)

# Collapse strata to eliminate lonely PSUs:
exdes<-collapse.strata(design=exdes,block.vars=~sr:procod)

# Now add new convenience variables to the design object:
  ## 'houdensity': to estimate households counts
  ## 'ones':       to estimate individuals counts
exdes<-des.addvars(exdes,
                   houdensity=ave(famcod,famcod,FUN = function(x) 1/length(x)),
                   ones=1)

# Let's see how it's possible to calibrate *simultaneously* on:
  # 1. the number of *individuals* crossclassified by sex, 5 age classes,
  #    and province;
  # 2. the number of *households* by region.

# First, for the purpose of running the example, let's generate some
# artificial population totals. We have only to get HT estimates for
# the auxiliary variables and perturb them randomly:
  # Get HT estimates of auxiliary variables:
  xx<-aux.estimates(design=exdes,calmodel=~houdensity+sex:age5c:procod-1,
                    partition=~regcod)

  # Add a random uniform perturbation to these numbers:
  set.seed(12345)   # Fix the RNG seed for reproducibility
  xx[,-1]<-round(xx[,-1]*runif(prod(dim(xx[,-1])),0.8,1.2))

# Now we have at hand our artificial population totals, and
# we can proceed with the calibration task:
excal<-e.calibrate(design=exdes,df.population=xx,calfun= "linear",
                   bounds=c(0,3),aggregate.stage=2)

# To perceive the effect of calibration, let's e.g. compare the HT and
# calibrated estimates of the average number of individuals per household
# at population level:
svystatR(exdes,~ones,~houdensity,vartype="cvpct")
svystatR(excal,~ones,~houdensity,vartype="cvpct")


########################################
# Calibrating on different patterns of #
# "incomplete" auxiliary information.  #
########################################

# Usually calibration constraints involve "complete auxiliary information",
# i.e. totals which are known either:
#     (i)  for the target population as a whole (e.g. total number of
#          employees working in italian active enterprises at a given date);
# or:
#     (ii) for each subpopulation belonging to a complete partition of
#          the target population (e.g. number of male and female people
#          residing in Italy at a given date).
# 
# Anyway, it may happen sometimes that the available auxiliary information
# is actually "incomplete", i.e. one doesn't know all the totals for all the
# subpopulations in a partition, but rather only for some of them. As an
# example, suppose marital status has categories "married", "unmarried",
# and "widowed" and that one only knows the number of "unmarried" people.
#
# In what follows I show how you can use ReGenesees to handle a calibration
# task on "incomplete" auxiliary information.

  #####################
  # A simple example. #
  #####################

  # Load household data:
  data(data.examples)

  # Define the survey design:
  des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
       weights=~weight)

  # Suppose you only know the number of "unmarried" people (let's say 398240)
  # but you ignore "married" and "widowed" totals, and you want to calibrate
  # on this incomplete information.

  # First, add to the survey design a new numeric variable with value 1
  # for unmarried people and 0 otherwise:
  des<-des.addvars(des,unmarried=as.numeric(marstat=="unmarried"))

  # Second, prepare a template to store the known "unmarried" people count:
  pop<-pop.template(des,calmodel=~unmarried-1)

  # Third, fill the template with the known total:
  pop[1,1]<-398240

  # Fourth, calibrate:
  descal<-e.calibrate(des,pop)

  # Now test that only "unmarried" estimated total has 0 percent CV:
  Zapsmall(svystatTM(descal,~marstat,vartype="cvpct"))

  # ...as it must be.


  ###############################
  # A more complicated example. #
  ###############################

  # Load sbs data:
  data(sbs)

  # Define the survey design:
  sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

  # Suppose you want to calibrate on the following "incomplete" known totals:
    # 1. enterprises counts by nace.macro
    # 2. enterprises counts by dom3 ONLY inside nace.macro 'Industry'
    # 3. total of y by emp.cl ONLY inside nace.macro 'Commerce'

  # First, add to the survey design new variables identifying the domains
  # where "incomplete" totals 2. and 3. are known:
  ## 2. -> nace.macro = 'Industry'
  sbsdes<-des.addvars(sbsdes,Industry=as.numeric(nace.macro=="Industry"))
  ## 3. -> nace.macro = 'Commerce'
  sbsdes<-des.addvars(sbsdes,Commerce=as.numeric(nace.macro=="Commerce"))

  # Do the same for the sampling frame:
  ## 2. -> nace.macro = 'Industry'
  sbs.frame$Industry=as.numeric(sbs.frame$nace.macro=="Industry")
  ## 3. -> nace.macro = 'Commerce'
  sbs.frame$Commerce=as.numeric(sbs.frame$nace.macro=="Commerce")

  # Second, prepare a template to store the totals listed in 1., 2. and 3.;
  # to this purpose one can e.g. compute HT estimates of the involved auxiliary
  # variables:
  Xht<-aux.estimates(design=sbsdes,
                     calmodel=~nace.macro+Industry:dom3+Commerce:y:emp.cl-1)
  Xht

  # Third, use the structure above to compute actual population totals
  # from the sampling frame:
  pop <- fill.template(universe=sbs.frame,template=Xht)
  pop

  # Fourth, calibrate:
  sbscal <- e.calibrate(design=sbsdes,df.population=pop)

  # Test1: nace.macro counts have 0 CVs:
  test1<-svystatTM(design=sbscal,y=~nace.macro,vartype="cvpct")
  test1

  # Test2: only 'Industry' macrosector has 0 CVs for dom3 counts:
  test2<-svystatTM(design=sbscal,y=~dom3,by=~nace.macro,vartype="cvpct")
  Zapsmall(test2)

  # Test3: only 'Commerce' macrosector has 0 CVs for y total by emp.cl:
  test3<-svystatTM(design=sbscal,y=~y,by=~emp.cl:nace.macro,vartype="cvpct")
  Zapsmall(test3)
}
\keyword{survey}