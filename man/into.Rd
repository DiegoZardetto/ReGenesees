\name{\%into\%}
\alias{\%into\%}

\title{Compress Nested Factors}

\description{
The special binary operator \code{\%into\%} transforms nested factors in such a way as to reduce the dimension and/or the sparsity of the model matrix of a calibration problem.
}

\usage{
inner \%into\% outer
"\%into\%"(inner, outer)
}

\arguments{
  \item{inner}{Factor with levels nested into \code{outer} (see \sQuote{Details}).}
  \item{outer}{Factor whose levels are an aggregation of those in \code{inner} (see \sQuote{Details}).}
}

\details{
Arguments \code{inner} and \code{outer} must be both \code{factor}s and must have the same \code{length}. Moreover, \code{inner} has to be \emph{strictly nested} into \code{outer}. Nesting is defined by treating elements in \code{inner} and \code{outer} as if they were positionally tied (i.e. as if they belonged to columns of a given data frame). The definition is as follows:\cr
\cr
\code{inner} and \code{outer} are strictly nested if, and only if, 1) every set of equal elements in \code{inner} correspond to a set of equal elements in \code{outer}, and 2) \code{inner} has \emph{more} non-empty levels than \code{outer}.\cr
\cr
If \code{inner} and \code{outer} do not fulfill the conditions above, evaluating \code{inner \%into\% outer} gives an error.

Suppose \code{inner} is actually nested into \code{outer} and define \code{inner.in.outer <- inner \%into\% outer}. The output factor \code{inner.in.outer} is built by recoding \code{inner} levels in such a way that each of them is mapped into the integer which represents its order inside the corresponding level of \code{outer} (see \sQuote{Examples}). As a consequence, the levels of \code{inner.in.outer} will be \code{1:n.max}, being \code{n.max} the \emph{maximum} number of levels of \code{inner} tied to a level of \code{outer}. Since this number is generally considerably smaller than the number of levels of \code{inner}, \code{inner.in.outer} can be seen as a \emph{compressed} representation of \code{inner}. Obviously, compression comes at a price: indeed \code{inner.in.outer} can now be used to identify a level of \code{inner} only \emph{inside} a given level of \code{outer} (see \sQuote{Examples}).

The usefulness of the \code{\%into\%} operator emerges in the calibration context. As we already documented in \code{\link{e.calibrate}}, factorizing a calibration problem (i.e. exploiting the \code{partition} argument of \code{e.calibrate}) determines a significant reduction in computation complexity, especially for big surveys. Now, it is sometimes the case that a calibration model is actually factorizable, even if this property is not self-apparent, due to factor nesting. In such cases, anyway, trying naively to factorize the \code{outer} variable(s) typically leads to very big and sparse model matrices (as well as population totals data frames), with the net result of washing-out the expected efficiency gain. A better alternative is to exploit the \code{\%into\%} operator in order to \emph{compress} the \code{inner} variable in such a way that the \code{outer} variable can be actually factorized \emph{without} giving rise to huge and sparse matrices. Section \sQuote{Examples} reports some practical illustration of the above line of reasoning.
}

\value{
A factor with levels \code{1:n.max}, being \code{n.max} the \emph{maximum} number of levels of \code{inner} tied to a level of \code{outer}.
}

\author{Diego Zardetto}

\seealso{
Further examples can be found in the \code{\link{fill.template}} help page.
}

\examples{
#################################################
## General properties of the \%into\% operator.   #
#################################################
 # First build a small data frame with 2 nested factors representing
 # regions and provinces:
dd <- data.frame(
                 reg  = factor( rep(LETTERS[1:3], c(6, 3, 1)) ),
                 prov = factor( rep(letters[1:6], c(3, 2, 1, 2, 1, 1)) )
                )
dd

 # Since prov is strictly nested into reg we can compute:
prov.in.reg <- dd$prov \%into\% dd$reg
prov.in.reg

 # Note that prov.in.reg has 3 levels because, as can be seen from dd,
 # the maximum number of provinces inside regions is 3. Thus prov.in.reg
 # is actually a compressed version of dd$prov (whose levels were 6)
 # but, obviously, it can now be used to identify a province only inside
 # a given region. This means that the the two factors below are identical (up
 # to levels' labels):
dd$prov
interaction(prov.in.reg,dd$reg,drop=TRUE)
 
 # Note that all the statements below generate errors:
\dontrun{
dd$reg  \%into\% dd$prov
dd$reg  \%into\% dd$reg
dd$prov \%into\% dd$prov
}

######################################################################
## A more useful (and complex) example from the calibration context. #
######################################################################
 # First define a design object:
data(data.examples)
exdes <- e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
weights=~weight)

 # Now suppose you have to perform a calibration process which
 # exploits the following known population totals:
 # 1) Joint distribution of sex and age10c (age in 10 classes)
 #    at the region level;
 # 2) Joint distribution of sex and age5c (age in 5 classes)
 #    at the province level;
 #
 # The auxiliary variables corresponding to the population totals above
 # can be symbolically represented by a calibration model like the following:
 # ~(procod:age5c + regcod:age10c - 1):sex
 #
 # At first sight it seems that only the sex variable can be factorized
 # in the model above. However if one observe that regions are an aggregation
 # of provinces, one realizes that also the regcod variable can be factorized.
 # Similarly, since categories of age5c are an aggregation of categories of
 # age10c, age5c can be factorized too. In both cases, using the \%into\%
 # operator will save computation time and memory usage. 
 # Let us see it in practice:
 #
 ## 1) Global calibration (i.e. calmodel=~(procod:age5c + regcod:age10c - 1):sex,
  # no partition variable, known totals stored in pop07):
t<-system.time(
               cal07<-e.calibrate(design=exdes,df.population=pop07,
                      calmodel=~(procod:age5c + regcod:age10c - 1):sex,
                      calfun="logit",bounds=c(0.2,1.8))
               )

 ## 2) Partitioned calibration on the self evident variable sex only
  # (i.e. calmodel=~procod:age5c + regcod:age10c - 1, partition=~sex,
  # known totals stored in pop07p):
tp<-system.time(
                cal07p<-e.calibrate(design=exdes,df.population=pop07p,
                        calmodel=~procod:age5c + regcod:age10c - 1,partition=~sex,
                        calfun="logit",bounds=c(0.2,1.8))
                )

 ## 3) Full partitioned calibration on variables sex, regcod and age5c
  # by exploiting \%into\%.
  # First add to the design object the new compressed factor variables
  # involving nested factors, namely provinces inside regions...
exdes<-des.addvars(exdes,procod.in.regcod=procod \%into\% regcod)
  
  # ...and age10c inside age5c:
exdes<-des.addvars(exdes,age10c.in.age5c=age10c \%into\% age5c)

  # Now calibrate exploiting the new variables
  # (i.e. calmodel=~procod.in.regcod + age10c.in.age5c - 1,
  # partition=~sex:regcod:age5c, known totals stored inside cal07pp)
tpp<-system.time(
                 cal07pp<-e.calibrate(design=exdes,df.population=pop07pp,
                          calmodel=~procod.in.regcod + age10c.in.age5c - 1,
                          partition=~sex:regcod:age5c,
                          calfun="logit",bounds=c(0.2,1.8))
                )

 # Now compare execution times:
t
tp
tpp

 # thus, tpp < tp < t, as expected.
 # Notice also that we obtained identical calibrated weights:
all.equal(weights(cal07),weights(cal07p))
all.equal(weights(cal07),weights(cal07pp))

 # as it must be.
}
\keyword{survey}