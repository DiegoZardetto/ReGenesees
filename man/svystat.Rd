\name{svystat}
\alias{svystat}
\alias{gvf.input.gr}
\alias{plot.gvf.input.gr}
\alias{SE.svystat.gr}
\alias{VAR.svystat.gr}
\alias{cv.svystat.gr}
\alias{deff.svystat.gr}
\alias{confint.svystat.gr}
\alias{coef.svystat.gr}


\title{Compute Many Estimates and Errors in Just a Single Shot}

\description{
Computes many estimates and errors (e.g. for disparate estimation domains) in just a single shot, primarily to use them in fitting GVF models. Can handle estimators of all kinds.
}

\usage{
svystat(design, kind = c("TM", "R", "S", "SR", "B", "Q", "L", "Sigma", "Sigma2"),
        by = NULL, group = NULL, forGVF = TRUE,
        combo = -1, ...)

\method{plot}{gvf.input.gr}(x, ...)

\method{coef}{svystat.gr}(object, ...)
\method{SE}{svystat.gr}(object, ...)
\method{VAR}{svystat.gr}(object, ...)
\method{cv}{svystat.gr}(object, ...)
\method{deff}{svystat.gr}(object, ...)
\method{confint}{svystat.gr}(object, ...)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{kind}{\code{character} specifying the summary statistics function to call: it may be \code{'TM'} (i.e. \code{svystatTM}, the default),
              \code{'R'} (i.e. \code{svystatR}), \code{'S'} (i.e. \code{svystatS}), \code{'SR'} (i.e. \code{svystatSR}),
              \code{'B'} (i.e. \code{svystatB}), \code{'Q'} (i.e. \code{svystatQ}), \code{'L'} (i.e. \code{svystatL}),
              \code{'Sigma'} (i.e. \code{svySigma}), and \code{'Sigma2'} (i.e. \code{svySigma2}).}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL}
           (the default option) estimates refer to the whole population.}
  \item{group}{Formula specifying a partition of the population into "groups": the output will be returned separately for each group.
               If \code{NULL} (the default option) the output is returned as a whole.}
  \item{forGVF}{Select \code{TRUE} (the default) if you want to use the output to fit a GVF model. Otherwise, the output will be
                simply a set of summary statistics objects.}
  \item{combo}{An \code{integer} which is only meaningful if \code{by} is passed. Requests to compute outputs for all the domains
               determined by crossing the \code{by} variables \emph{up to} a given order (see \sQuote{Details}).}
  \item{\dots}{For function \code{svystat}, additional arguments to the summary statistic function implied by \code{kind}.
               Otherwise, further arguments passed to or from other methods.}
  \item{x}{The object of class \code{gvf.input.gr} to plot.}
  \item{object}{An object of class \code{svystat.gr} containing survey statistics.}
}

\details{
This function can compute \emph{all} the summary statistics provided by \pkg{ReGenesees}, and is principally meant to return a lot of them in just a single shot.

If \code{forGVF = TRUE} the output will be ready to feed \pkg{ReGenesees} GVF fitting infrastructure, otherwise it will consist simply of a set of summary statistic objects.

Use argument \code{kind} to specify the summary statistic you need. The default value \code{'TM'} selects function \code{\link{svystatTM}}, which yields Totals and Means. All the arguments needed by the summary statistic function implied by \code{kind} (e.g. argument \code{y} for \code{svystatTM} when \code{kind = 'TM'}) will be passed on through argument \sQuote{\code{\ldots}}.

As usual in summary statistics, argument \code{by} can be used to request domain estimates.

The \code{group} formula (if any) specifies a way of partitioning the population into groups: the output will be reported separately for each group. In the GVF context, a \emph{\dQuote{grouped}} output will permit to fit \emph{separate} GVF models inside different groups (and hence to compute separate variance predictions for different groups).

Note that \code{group} and \code{by} share identical syntax and semantics as model formulae, despite they have different purposes in function \code{svystat} (as explained above).

Parameter \code{combo} is \emph{only} meaningful if \code{by} is passed. Its purpose is to allow computing estimates and errors simultaneously for many estimation domains.

If the \code{by} formula involves \code{n} variables, specifying \code{combo = m} requests to compute outputs for all the domains determined by all the interactions of \code{by} variables \emph{up to} order \code{m} (with \code{-1 <= m <= n}), as follows:
\preformatted{
  COMBO         MEANING
  m = -1.......'no combo', i.e. treat 'by' formula as usual (the default);
  m =  0.......'order zero' combination, i.e. just a single domain:
                the whole population;
  m =  1.......'order zero' plus 'order one' combinations, the latter being
                all the marginal domains defined by 'by' variables;
  m =  n........combinations of any order, the maximum being the one with
                all 'by' variables interacting simultaneously.
}

The \code{plot} method can be used \emph{only} when \code{forGVF = TRUE} and produces a matrix (or many matrices, if \code{group} is passed) of scatterplots with polynomial smoothers.

Methods \code{\link{coef}}, \code{\link{SE}}, \code{\link{VAR}}, \code{\link{cv}}, \code{\link{deff}}, and \code{\link{confint}} can be used \emph{only} when \code{forGVF = FALSE}, to extract estimates and variability statistics. 
}

\value{
An object storing estimates and errors, whose detailed structure depends on input parameters' values.

If \code{forGVF = FALSE}, a set of summary statistics possibly stored into a list (with class \code{svystat.gr} in the most general case).

If \code{forGVF = TRUE} and argument \code{group} is \emph{not} passed, an object of class \code{\link{gvf.input}}.

If \code{forGVF = TRUE} and argument \code{group} is passed, an object of class \code{gvf.input.gr}. This is a list of objects of class \code{\link{gvf.input}}, each one pertaining to a different population group.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{estimator.kind}} to assess what kind of estimates are stored inside a survey statistic object, \code{\link{gvf.input}} as an alternative to prepare the input for GVF model fitting, \code{\link{GVF.db}} to manage \pkg{ReGenesees} archive of registered GVF models, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.

}

\examples{
# Load sbs data:
data(sbs)

# Create a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)


##########################################################################
# svystat as an alternative way to compute 'ordinary' summary statistics #
##########################################################################
## Total number of employees
svystat(sbsdes,y=~emp.num,forGVF=FALSE)
# equivalent to:
svystatTM(sbsdes,y=~emp.num)

## Average number of employees per enterprise
svystat(sbsdes,y=~emp.num,estimator="Mean",forGVF=FALSE)
# equivalent to:
svystatTM(sbsdes,y=~emp.num,estimator="Mean")

## Average value added per employee by economic activity macro-sector
## (nace.macro):
svystat(sbsdes,kind="R",num=~va.imp2,den=~emp.num,by=~nace.macro,forGVF=FALSE)
# equivalent to:
svystatR(sbsdes,num=~va.imp2,den=~emp.num,by=~nace.macro)

## Counts of employees by classes of number of employees (emp.cl) crossed
## with economic activity macro-sector (nace.macro):
svystat(sbsdes,y=~emp.num,by=~emp.cl:nace.macro,forGVF=FALSE)
# equivalent to:
svystatTM(sbsdes,y=~emp.num,by=~emp.cl:nace.macro)

## Provided forGVF = FALSE, you can use estimator.kind on svystat output:
stat<-svystat(sbsdes,kind="R",num=~va.imp2,den=~emp.num,by=~emp.cl:nace.macro,
      group=~region,forGVF=FALSE)
stat
estimator.kind(stat,sbsdes)


##########################################################
# Understanding syntax and semantics of argument 'combo' #
##########################################################
# Load household data:
data(data.examples)

# Create a design object:
houdes<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
        weights=~weight)

# Add convenience variable 'ones' to estimate counts:
houdes<-des.addvars(houdes,ones=1)

## To facilitate understanding, let's for the moment keep forGVF = FALSE.
## Let's use estimates and errors of counts of individuals by sex and
## five age classes (age5c):
svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE)

## Now let's play with argument 'combo':
   # combo = -1
   # -> 'no combo', i.e. treat 'by' formula as usual
   svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE,combo=-1)

   # combo = 0
   # -> 'order zero' combination, i.e. just a single domain: the whole population
   svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE,combo=0)

   # combo = 1
   # -> 'order zero' plus 'order one' combinations, the latter being all the
   #    marginal domains defined by 'by' variables
   svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE,combo=1)

   # combo = 2
   # -> since 'by' has 2 variables, this means combinations of any order up to
   #    the maximum
   svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE,combo=2)

   # combo = 3
   # -> yields an error, as 'combo' cannot exceed the number of 'by' variables
   #    (2 in this example)
\dontrun{
   svystat(houdes,y=~ones,by=~age5c:sex,forGVF=FALSE,combo=3)
}


######################################################################
# svystat as an alternative way to prepare input data for GVF models #
######################################################################
## The same estimates and errors of the last example above, now with
## forGVF = TRUE: note the different output data format
svystat(houdes,y=~ones,by=~age5c:sex,combo=2)

## Note that the agile command above is indeed equivalent to the following
## lengthier, cumbersome statement:
gvf.input(houdes,
          svystatTM(houdes,y=~ones),
          svystatTM(houdes,y=~ones,by=~age5c),
          svystatTM(houdes,y=~ones,by=~sex),
          svystatTM(houdes,y=~ones,by=~age5c:sex)
          )


  ################################################
  # Using argument 'group' to prepare input data #
  # for separate GVF models                      #
  ################################################
  ## The same estimates and errors of the last example above, now prepared
  ## separately for different regions (regcod):
  svystat(houdes,y=~ones,by=~age5c:sex,combo=2,group=~regcod)

  ## Again the same estimates and errors, prepared separately for groups
  ## defined crossing marital status (marstat) and region:
  svystat(houdes,y=~ones,by=~age5c:sex,combo=2,group=~marstat:regcod)
 
  ## NOTE: Output has class "gvf.input.gr". This will tell ReGenesees' GVF
  ##       fitting facilities to handle estimates and errors pertaining to
  ##       different groups independently of each other.


## NOTE: Parameter combo allows svystat to gather a huge amount of estimates and
##       errors in just a single slot, as the number of estimation domains grows
##       exponentially with the number of by variables.
##       See, for instance, the following example:
out <- svystat(houdes,y=~ones,by=~age5c:marstat:sex:regcod,combo=4)
dim(out)
head(out)
plot(out)


##################################################
# Minor details: accessor functions and plotting #
##################################################
  ## Accessor functions work only when forGVF = FALSE
# Average value added per employee by nace.macro:
out <- svystat(sbsdes,kind="R",num=~va.imp2,den=~emp.num,by=~nace.macro,forGVF=FALSE)
out
# Access CV values and confidence intervals:
cv(out)
confint(out)

# The same as above, separately for regions:
out <- svystat(sbsdes,kind="R",num=~va.imp2,den=~emp.num,by=~nace.macro,group=~region,forGVF=FALSE)
out
# Access CV values and confidence intervals:
cv(out)
confint(out)


  ## Plot function works only when forGVF = TRUE
# Counts of individuals by sex, marstat and age5c, and all their interactions:
out <- svystat(houdes,y=~ones,by=~age5c:marstat:sex,combo=3)
# Plot GVF input:
plot(out)

# The same as above, grouped by region:
out <- svystat(houdes,y=~ones,by=~age5c:marstat:sex,combo=3,group=~regcod)
# Plot GVF inputs, separately by groups (regions):
plot(out)
}

\keyword{survey}