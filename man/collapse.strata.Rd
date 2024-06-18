\name{collapse.strata}
\alias{collapse.strata}

\title{Collapse Strata Technique for Eliminating Lonely PSUs}

\description{
Modifies a stratified design containing lonely PSUs by collapsing its design strata into superstrata.
}

\usage{
collapse.strata(design, block.vars = NULL, sim.score = NULL)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling
                design metadata.}
  \item{block.vars}{Formula specifying blocking variables: only strata belonging to the same block will be
                    aggregated (see \sQuote{Details}). If \code{NULL} (the default option) no constraints
                    will be imposed.}
  \item{sim.score}{Formula specifying a similarity score for strata: lonely strata will be paired
                   with the most similar stratum in each block (see \sQuote{Details}). If \code{NULL} (the default option)
                   random pairs will be formed.}
}

\details{
Lonely PSUs (i.e. PSUs which are alone inside a not self-representing stratum) are a concern from the viewpoint of variance estimation. As a general solution, the \pkg{ReGenesees} package can handle the lonely PSUs problem by setting proper variance estimation options (see \code{\link{ReGenesees.options}}). The \code{collapse.strata} function implements a widely used alternative: the so called collapsed strata technique. The basic idea is to build artificial \emph{\dQuote{superstrata}} by aggregating strata containing lonely PSUs to other strata, and then to use such superstrata for variance estimation (see e.g. [Wolter 85] and [Rust, Kalton 87]).

The optional argument \code{block.vars} identifies \emph{\dQuote{blocking variables}} that can be used to constrain the way lonely strata are collapsed to form superstrata. More specifically: first, blocking variables are used to partition sample data in \emph{\dQuote{blocks}} via factor crossing, then, only lonely strata belonging to the same block are aggregated. If \code{block.vars=NULL} (the default option), no constraint will act on collapsing. The \code{design} variables referenced by \code{block.vars} (if any) should be of type \code{factor}. Errors will be raised if (i) blocks cut across strata, or (ii) \code{block.vars} generate any non-aggregable strata (i.e. lonely strata which are a singleton inside a block).

The optional argument \code{sim.score} can be used to specify a similarity score for strata aggregation. This means that each lonely stratum will be collapsed with the stratum that has the most similar value of variable \code{sim.score} inside the block. Thus the similarity of two strata is actually measured by the (absolute value of the) difference among the corresponding \code{sim.score} values. Only one \code{design} variable can be referenced by the \code{sim.score} formula: (i) it must be of type \code{numeric}, (ii) it must be constant inside each stratum, and (iii) it should be positive (otherwise its \code{abs()} will be silently used). Note that if no similarity score is specified (i.e. \code{sim.score=NULL}), the achieved strata aggregation will depend on the ordering of input sample data in \code{design}.

The collapsing algorithm will, whenever possible, build superstrata by pairing a lonely stratum to another not-yet-aggregated stratum. Therefore, in general, superstrata will contain only two design strata. Rare exceptions can arise, e.g. due to constraints, with at most three design strata inside a superstratum. The choice to collapse strata in pairs has been taken because it is known to be appropriate for large-scale surveys with many strata (at least for national level estimates, see e.g. [Rust, Kalton 87]).

The \code{collapse.strata} function handles correctly finite population corrections. If \code{design} has been built by passing strata sampling fractions via the \code{fpc} argument, the function re-computes sampling fractions inside superstrata by exploiting the achieved mapping of strata to superstrata and the \code{fpc} slot of \code{design}. 
}

\section{Strata Collapse Process Diagnostics}{
As already observed in the \sQuote{Details} Section, there are three non trivial reasons why function \code{collapse.strata} can run into errors: (1) the blocks cut across strata, (2) some blocks contain a stratum needing to be aggregated while this stratum happens to be the only one inside the block, (3) the similarity score for strata aggregation varies inside strata. In order to help the user to identify such data anomalies, hence taking a step forward to eliminate them, every call to \code{collapse.strata} generates, by side effect, a diagnostics data structure named \code{clps.strata.status} into the \code{.GlobalEnv} (see \sQuote{Examples}).

The \code{clps.strata.status} list has three components: the first reports the error message, the second stores a vector identifying the data subsets that have been hit by the anomaly, the third reports the call to \code{collapse.strata} that generated the list. For instance, when error condition (1) holds, the second element of \code{clps.strata.status} identifies the strata that are cut by blocks; if, instead, error condition (2) holds, the second element of the list identifies the blocks containing non-aggregable strata.

It must be stressed that \emph{every call} to \code{collapse.strata} generates the \code{clps.strata.status} list, \emph{even} when the strata collapsing process ends \emph{successfully}. In such cases, the first element of the list reports the number of lonely strata that have undergone aggregation, whereas the second is a useful data frame (named \code{clps.table}) mapping collapsed strata to superstrata. To be more specific: each row of \code{clps.table} identifies a stratum that has been mapped to a superstratum, while the columns of \code{clps.table} give: (i) the block to which the stratum belongs, (ii) the stratum name, (iii) a flag indicating if the stratum was lonely or not, (iv) the name of the superstratum to which it has been mapped.
}

\section{Methodological Warning}{
A warning must be emphasized: strata similarity score \code{sim.score} should be based on prior knowledge and/or on expectations on \emph{true} values of stratum means for the variable(s) to be estimated, not on current sample data. Indeed, building \code{sim.score} by estimating stratum means with the current sample can lead to severe \emph{underestimation} of sampling variance, i.e. to too tight confidence intervals.
}

\value{
An object of the same class as \code{design}, without strata containing lonely PSUs.
}

\references{
Wolter, K.M. (2007) \emph{\dQuote{Introduction to Variance Estimation}}, Second Edition, Springer-Verlag, New York.

Rust, K., Kalton, G. (1987) \emph{\dQuote{Strategies for Collapsing Strata for Variance Estimation}}, Journal of Official Statistics, Vol. 3, No. 1, pp. 69-81.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{ReGenesees.options}} for a different way to handle the lonely PSUs problem (namely by setting variance estimation options).
}

\examples{
##############################################
# Explore alternative collapsing strategies. #
##############################################

# Build a survey design with lonely PSU strata:
data(data.examples)
exdes <- e.svydesign(data= example, ids= ~ towcod+famcod,
         strata= ~ stratum, weights= ~ weight)
exdes

# Explore 3 possible collapsing strategies:
  # 1) Aggregate lonely strata by forming random pairs
exdes.clps1 <- collapse.strata(exdes)
exdes.clps1

  # 2) Aggregate lonely strata in pairs under constraints:
  #    i.  aggregated strata must be both not self-representing
  #    ii. aggregated strata must belong to the same province (which
  #        is appropriate if e.g. provinces are planned estimation domains)
exdes.clps2 <- collapse.strata(exdes,~sr:procod)
exdes.clps2

  # 3) A WRONG strategy: compute strata similarity score by using
  #    sample estimates of the interest variable (y1) inside strata:
old.op <- options("RG.lonely.psu"="remove")
stat.score <- svystatTM(design= exdes, ~y1, by= ~ stratum)
options(old.op)
exdes2<-des.addvars(exdes,
                    sim.score=stat.score[match(stratum,stat.score$stratum),2])
exdes.clps3 <- collapse.strata(exdes2,~sr:procod,~sim.score)
exdes.clps3

# Compute total estimates of y1 at the province level
# for all 3 designs with collapsed strata:
stat.clps1 <- svystatTM(design= exdes.clps1, y= ~ y1, by= ~ procod,
              estimator= "Total", vartype= "cvpct")
stat.clps2 <- svystatTM(design= exdes.clps2, y= ~ y1, by= ~ procod,
              estimator= "Total", vartype= "cvpct")
stat.clps3 <- svystatTM(design= exdes.clps3, y= ~ y1, by= ~ procod,
              estimator= "Total", vartype= "cvpct")

# Compute the same estimates by using two alternatives
# to handle lonely PSUs:
  # "adjust" option
old.op <- options("RG.lonely.psu"="adjust")
stat.adj <- svystatTM(design= exdes, y= ~ y1, by= ~ procod,
            estimator= "Total", vartype= "cvpct")
options(old.op)
 # "average" option
old.op <- options("RG.lonely.psu"="average")
stat.ave <- svystatTM(design= exdes, y= ~ y1, by= ~ procod,
            estimator= "Total", vartype= "cvpct")
options(old.op)

# Lastly, compare achieved estimates for CV percentages:
stat.clps1
stat.clps2
stat.clps3
stat.adj
stat.ave

# Thus the qualitative features are as expected: the "adjust" option
# tends to give conservative sampling variance estimates, the WRONG collapsing
# strategy 3) tends to underestimate sampling variance, while other methods
# give results in-between those extrema.


###########################################################
# A simple way for defining the strata similarity scores. #
###########################################################
# Suppose that strata have been clustered in groups of similar
# strata. You can, then, use the integer codes of the factor
# variable identifying the clusters as a similarity score.
# You can do as follows:

  # Load some data:
    data(fpcdat)

  # Build a design object:
    fpcdes<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w)
    fpcdes

  # As we deliberately omitted to specify fpcs, this design
  # has 2 lonely strata out of 5:
    find.lon.strata(fpcdes)

  # Now, suppose that factor variable pl.domain identifies clusters of
  # similar strata... 
    table(fpcdat$stratum,fpcdat$pl.domain)

  # ...hence, the similarity score can be obtained simply...
    fpcdes<-des.addvars(fpcdes,score=unclass(pl.domain))

  # ...and readily be used to drive the strata collapsing:
    fpcdes.clps<-collapse.strata(fpcdes,sim.score=~score)
    fpcdes.clps
    clps.strata.status

  # As we expected from the groups defined by pl.domain, lonely stratum S.2
  # has been paired to S.3, and lonely stratum S.5 to S.4.

  # Should we have omitted to specify a similarity score, we would have
  # obtained different superstrata:
    fpcdes.clps2<-collapse.strata(fpcdes)
    fpcdes.clps2
    clps.strata.status


#################################################################
# Few examples to inspect the clps.strata.status list generated # 
# for diagnostics purposes.                                     #
#################################################################

  # 1) Ill defined blocks: cutting across strata:
\dontrun{
clps.err1 <- collapse.strata(exdes,~sex)
}
clps.strata.status

  # 2) Ill defined blocks: generating non-aggregable strata
\dontrun{
clps.err2 <- collapse.strata(exdes,~regcod:stratum)
}
clps.strata.status

  # 3) Successful collapsing: explore strata to superstrata mapping
exdes.ok <- collapse.strata(exdes, ~sr:regcod:procod)
clps.strata.status
}
\keyword{survey}