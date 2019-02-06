\name{contrasts.RG}
\alias{contrasts.RG}
\alias{contrasts.off}
\alias{contrasts.reset}
\alias{contr.off}
\alias{calmodel}

\title{Set, Reset or Switch Off Contrasts for Calibration Models}

\description{
These functions control the way \pkg{ReGenesees} translates a symbolic calibration model (as specified by the \code{calmodel} formula in \code{\link{e.calibrate}}, \code{\link{pop.template}}, \code{\link{fill.template}}, \code{\link{aux.estimates}}, \ldots) to its numeric encoding (i.e. the model-matrix used by the internal algorithms to perform actual computations).
}
\usage{
contrasts.RG()
contrasts.off()
contrasts.reset()
contr.off(n, base = 1, contrasts = TRUE, sparse = FALSE)
}
\arguments{
  \item{n}{Formally as in function \code{contr.treatment} (see \sQuote{Details}).}
  \item{base}{Formally as in function \code{contr.treatment} (see \sQuote{Details}).}
  \item{contrasts}{Fictitious, but formally as in function \code{contr.treatment}. (see \sQuote{Details})}
  \item{sparse}{Formally as in function \code{contr.treatment}. (see \sQuote{Details})}
}
\details{
All the calibration facilities in package \pkg{ReGenesees} transform symbolic calibration models (as specified by the user via \code{calmodel}) into numeric model-matrices. Factor variables occurring in \code{calmodel} play a special role in such transformations, as the encoding of a factor can (and, by default, do) \emph{depend} on the \emph{structure} of the formula in which it occurs. The \pkg{ReGenesees} functions documented below control the way factor levels are translated into auxiliary variables and mapped to columns of population totals data frames. The underlying technical tools are \code{\link{contrasts}} handling functions (see Section 'Technical Remarks and Warnings' for further details).

Under the calibration perspective, ordered and unordered factors appearing in \code{calmodel} must be treated the same way. This obvious constraint defines the \pkg{ReGenesees} default for contrasts handling. Such a default is silently set when loading the package. Moreover, you can set it also by calling \code{contrasts.RG()}. As can be understood by reading Section 'Technical Remarks and Warnings' below, the default setup can be seen as \strong{\dQuote{efficient-but-slightly-risky}}.

A call to \code{contrasts.off()} simply disables all contrasts and imposes a complete dummy coding of factors. Under this setup, all levels of factors occurring in \code{calmodel} generate a distinct model-matrix column, \emph{even if some of these columns can be linearly dependent}. To be very concise, the \code{contrasts.off()} setup can be seen as \strong{\dQuote{safe-but-less-efficient}} as compared to the default one (read Section 'Technical Remarks and Warnings' for more details).

Function \code{contr.off} is not meant to be called directly by users: it serves only the purpose of enabling the \code{contrasts.off()} setup.

A call to \code{contrasts.reset()} restores \R factory-fresh defaults for contrasts (which do distinguish ordered and unordered factors). Users may want to use this function after having completed a \pkg{ReGenesees} session, e.g. before switching to other \R functions relying on contrasts (such as \code{\link{lm}}, \code{\link{glm}}, \ldots).
}
\section{Technical Remarks and Warnings}{

\dQuote{\emph{[...] the corner cases of model.matrix and friends is some of the more impenetrable code in the R sources.}}\cr
\href{http://tolstoy.newcastle.edu.au/R/e17/help/12/01/2497.html}{Peter Dalgaard}\cr

Contrasts handling functions tell \R how to encode the model-matrix associated to a given model-formula on specific data (see, e.g., \code{\link{contr.treatment}}, \code{\link{contrasts}}, \code{\link{model.matrix}}, \code{\link{formula}}, and references therein). More specifically, contrasts control the way factor-terms and interaction-terms occurring in formulae get actually represented in the model matrix. For instance, \R (by default) avoids the complete dummy coding of a factor whenever it is able to understand, on the basis of the structure of the model-formula, that some of the factor levels would generate linearly dependent (i.e. redundant) columns in the model-matrix (see Section \sQuote{Examples}).

The usage of contrasts to build smaller, full-rank calibration model-matrices would be a good opportunity for \pkg{ReGenesees}, provided it comes \emph{without any information loss}. Indeed, smaller model-matrices mean less population totals to be provided by users, and higher efficiency in computations.

Unfortunately, few controversial cases have been signalled in which \R ability to "simplify" a model-matrix on the basis of the structure of the related model-formula seems to lead to strange, unexpected results (see, e.g., \href{http://tolstoy.newcastle.edu.au/R/e17/help/12/01/2375.html}{this R-help thread}). No matter whether such \R behaviour is or not an actual bug with respect to its impact on \R linear model fitting or ANOVA facilities, it surely represents a concern for \pkg{ReGenesees} with respect to calibration (see Section \sQuote{Examples}). The risk is the following: there could be rare cases in which exploiting \R contrasts handling functions inside \pkg{ReGenesees} ends up with a \emph{wrong} (i.e. incomplete) population totals template, and (eventually) with \emph{wrong} calibration results.

Though one could adopt several ad-hoc countermeasures to sterilize the risk described above while still taking advantage of contrasts (see Section \sQuote{Examples}), the choice of completely disabling contrasts via \code{contrasts.off()} would result in a \strong{100\% safety guarantee}. If computational efficiency is not a serious concern for you, switching off contrasts may determine the best \pkg{ReGenesees} setup for your analyses.
}
\author{
Diego Zardetto
}
\references{
Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical models}. Chapter 2 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

\href{http://tolstoy.newcastle.edu.au/R/e17/help/12/01/2375.html}{"Why does the order of terms in a formula translate into different models/model matrices?"}, R-help thread.
}
\seealso{
\code{\link{e.calibrate}}, \code{\link{pop.template}}, \code{\link{fill.template}}, and \code{\link{aux.estimates}} for the meaning and the usage of \code{calmodel} in \pkg{ReGenesees}. \code{\link{formula}}, \code{\link{model.matrix}}, \code{\link{contrasts}}, and \code{\link{contr.treatment}} to understand the role of contrasts in \R.
}
\examples{
######################
# Easy things first: #
######################

  # 1) When ReGenesees is loaded, its standard way of handling contrasts
  #    (i.e. no ordered-unordered factor distinction) is silently set:
options("contrasts")

  # 2) To switch off contrasts (i.e. apply always dummy coding to factors),
  #    simply type:
contrasts.off()

  # 3) To restore R factory-fresh defaults for contrasts, simply type:
contrasts.reset()

  # 4) To switch on again standard ReGenesees contrasts, simply type:
contrasts.RG()


#############################################################
# A simple calibration example to understand the effects of #
# switching off contrasts.                                  #
#############################################################

# Load sbs data:
data(sbs)

# Create a design object:
sbsdes<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,fpc=~fpc)

# Suppose you want to calibrate on the marginals of 'region' (a factor
# with 3 levels: "North", "Center", and "South") and 'dom3' (a factor
# with 4 levels: "A", "B", "C", and "D").
# Let's see how things go under the 'contrast on' (default) and 'contrasts off'
# setups:

  ########################################
  # 1) ReGenesees default: contrasts ON. #
  ########################################
  # As you see contrasts are ON:
  options("contrasts")

  # Build and fill the population totals template:
  temp1<-pop.template(data=sbsdes,calmodel=~region+dom3-1)
  pop1<-fill.template(universe=sbs.frame,template=temp1)

  # Now inspect the obtained known totals data.frame:
  pop1
  
  # As you see: (i) it has only 6 columns, and (ii) the "A" level of
  # factor 'dom3' is missing. This is because contrasts are ON, so that
  # R is able to understand that only 6 out of the 3 + 4 marginal counts
  # are actually independent. Indeed, the "A" counts...
  sum(sbs.frame$dom3=="A")
  
  # ...are actually redundant, since they can be deduced by pop1:
  sum(pop1[,1:3])-sum(pop1[,4:6]) 
  
  # Now calibrate:
  cal1<-e.calibrate(sbsdes,pop1)

  
  ##################################################
  # 2) Switch OFF contrasts: dummy coding for all! #
  ##################################################
  # To switch off contrasts simply call:
  contrasts.off()

  # Build and fill the population totals template:
  temp2<-pop.template(data=sbsdes,calmodel=~region+dom3-1)
  pop2<-fill.template(universe=sbs.frame,template=temp2)
  
  # Now inspect the obtained known totals data.frame:
  pop2
  
  # As you see: (1) it has now 7 columns, and (2) the "A" level of factor
  # 'dom3' has been resurrected. This is because contrasts are OFF,
  # so that each level of factors in calmodel are coded to dummies.
  
  # Now calibrate. Since only 6 out of 7 dummy auxiliary variables are
  # actually independent, the model.matrix computed by e.calibrate will not be
  # full-rank. As a consequence, e.calibrate would use the Moore-Penrose
  # generalized inverse (in practice, this could depend on the machine R
  # is running on):
  cal2<-e.calibrate(sbsdes,pop2)
  
  # Compare the calibration weights generated under setups 1) and 2):
  all.equal(weights(cal2),weights(cal1))
  
  # Lastly set back contrasts to ReGenesees default:
  contrasts.RG()


#############################################
# Weird results, risks and countermeasures. #
# ("When the going gets tough...")          #
#############################################

# Suppose you want to calibrate on: (A) the joint distribution of 'region' (a
# factor with 3 levels: "North", "Center", and "South") and 'nace.macro' (a
# factor with 4 levels: "Agriculture", "Industry", "Commerce", and "Services")
# and, at the same time, on (B) the total number of employees ('emp.num', a
# numeric variable) by 'nace.macro'.
#
# You rightly expect that 3*4 + 4 = 16 population totals are needed for such a
# calibration task. Indeed, knowing the enterprise counts for the 3*4 cells of
# the joint distribution (A) doesn't tell anything on the number of employees
# working in the 4 nace macrosectors (B), and vice-versa.
#
# Moreover, you might expect that calibration models:
  # (i)  calmodel = ~region:nace.macro  + emp.num:nace.macro - 1
  # (ii) calmodel = ~emp.num:nace.macro + region:nace.macro  - 1
#  
# should produce the same results.
# Unfortunately, WHEN CONTRASTS ARE ON, this is not the case: only model (i)
# leads to the expected, right results. Let's see.

  ###########################################
  # A strange result when contrasts are ON: #
  # the order of terms in calmodel matters! #
  ###########################################
  # As you see contrasts are ON:
  options("contrasts")

  # Start with (i) calmodel = ~region:nace.macro + emp.num:nace.macro - 1
  # Build and fill the population totals template:
  temp1<-pop.template(data=sbsdes,~region:nace.macro+emp.num:nace.macro-1)
  pop1<-fill.template(universe=sbs.frame,template=temp1)

  # Now inspect the obtained known totals data.frame:
  pop1

  # and verify it stores the right, expected number of totals (i.e. 16):
  dim(pop1)

  # Now calibrate:
  cal1<-e.calibrate(sbsdes,pop1)


  # Now compare with (ii) calmodel = ~emp.num:nace.macro + region:nace.macro - 1
  # Build and fill the population totals template:
  temp2<-pop.template(data=sbsdes,~emp.num:nace.macro+region:nace.macro-1)
  pop2<-fill.template(universe=sbs.frame,template=temp2)

  # First check if it stores the right, expected number of totals (i.e. 16):
  dim(pop2)

  # Apparently 4 totals are missing; let's inspect the known totals data.frame
  # to understand which ones:
  pop2

  # Thus we are missing the 4 'nace.macro' totals for 'region' level "North".
  # Everything goes as if R contrasts functions mistakenly treated the term
  # emp.num:nace.macro as a factor-factor interaction (i.e. a 2 way joint
  # distribution), which would have justified to eliminate the 4 missing totals
  # as redundant.

  # Notice that calibrating on pop2 would generate wrong results...
  cal2<-e.calibrate(sbsdes,pop2)

  # ...indeed the 4 estimates of 'nace.macro' for 'region' level "North" are not
  # actually calibrated (look at the magnitude of SE estimates):
  svystatTM(cal2,~region,~nace.macro)


  ################################################################
  # A possible countermeasure (still working with contrasts ON). #
  ################################################################
  # Empirical evidence tells that the weird case above is extremely rare
  # and that it manifests whenever a numeric (say X) and a factor (say F) both
  # interact with the same factor (say D), i.e. calmodel=~(X+F):D-1.
  #
  # The risky order-dependent nature of such models can be sterilized (while
  # still taking advantage of contrasts-driven simplifications for large,
  # complex calibrations) by using a numeric variable with values 1 for
  # all sample units.
  #
  # For instance, one could use variable 'ent' in the sbs data.frame, to
  # handle the (A) part of the calibration constraints. Indeed you may easily
  # verify that both the calmodel formulae below:
  # (i)  calmodel = ~ent:region:nace.macro  + emp.num:nace.macro - 1
  # (ii) calmodel = ~emp.num:nace.macro + ent:region:nace.macro  - 1
  #
  # produce exactly the same, right results.


  ##################################################################
  # THE ULTIMATE, 100% SAFE, COUNTERMEASURE: switch contrasts OFF! #
  ##################################################################
  # No contrasts means no model-matrix simplifications at all, hence
  # also no unwanted, wrong simplifications. Let's see:

  # To switch off contrasts simply call:
  contrasts.off()

  # Compare again, with contrasts OFF, the calibration models:
  # (i)  calmodel = ~region:nace.macro  + emp.num:nace.macro - 1
  # (ii) calmodel = ~emp.num:nace.macro + region:nace.macro  - 1
  
  # Build and fill the population totals templates:
  temp1<-pop.template(data=sbsdes,~region:nace.macro+emp.num:nace.macro-1)
  pop1<-fill.template(universe=sbs.frame,template=temp1)
  pop1

  temp2<-pop.template(data=sbsdes,~emp.num:nace.macro+region:nace.macro-1)
  pop2<-fill.template(universe=sbs.frame,template=temp2)
  pop2

  # Verify they store the same, right number of totals (i.e. 16):
  dim(pop1)
  dim(pop2)

  # Verify they lead to right calibrated objects...
  cal1<-e.calibrate(sbsdes,pop1)
  cal2<-e.calibrate(sbsdes,pop2)

  # ...with the same calibrated weights:
  all.equal(weights(cal2),weights(cal1))
  
  # Lastly set back contrasts to ReGenesees default:
  contrasts.RG()
}