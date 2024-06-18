\name{ReGenesees.options}
\alias{ReGenesees.options}
\alias{RG.lonely.psu}
\alias{RG.ultimate.cluster}
\alias{RG.adjust.domain.lonely}
\title{Variance Estimation Options for the ReGenesees Package}
\description{
  This help page documents the options that control the behaviour of
  the \pkg{ReGenesees} package with respect to standard error estimation.
  }
  
  \details{
    The \pkg{ReGenesees} package provides four options for variance
    estimations which can be freely set and modified by the user:\cr
    
    - \code{RG.ultimate.cluster}\cr
    - \code{RG.lonely.psu}\cr
    - \code{RG.adjust.domain.lonely}\cr
    - \code{RG.warn.domain.lonely}\cr
    
    When \code{options("RG.ultimate.cluster")} is \code{TRUE},
    the \pkg{ReGenesees} package adopts the so called \emph{\dQuote{Ultimate
    Cluster Approximation}} [Kalton 79]. Under this approximation, the overall
    sampling variance for a multistage sampling design is estimated by
    taking into account only the contribution arising from the estimated PSU
    totals (thus simply ignoring any available information about subsequent
    sampling stages). For without replacement sampling designs, this approach
    is known to underestimate the true multistage variance, while - at the
    same time - overestimating its true first-stage component. Anyway, the
    underestimation error becomes negligible if the PSUs' sampling fractions
    across strata are very small. When sampling with replacement, the Ultimate
    Cluster approach is no longer an approximation, but rather an exact result.
    Hence, be \code{options("RG.ultimate.cluster")} \code{TRUE} or \code{FALSE},
    if one does not specify first-stage finite population corrections, \pkg{ReGenesees}
    will produce exactly the same variance estimates.
    
    When \code{options("RG.ultimate.cluster")} is \code{FALSE},
    each sampling stage contributes and variances get estimated by means
    of a recursive algorithm [Bellhouse, 85] inherited and adapted from
    package \pkg{survey} [Lumley 06]. Notice that the results obtained
    by choosing this option can differ from the one that would be obtained
    under the "Ultimate Cluster Approximation" \emph{only if} first-stage
    finite population corrections are specified.
    
    Lonely PSUs (i.e. PSUs which are alone inside a not self-representing
    stratum) are a concern from the viewpoint of variance estimation. The
    suggested \pkg{ReGenesees} facility to handle the lonely PSUs problem is
    the strata aggregation technique (see e.g. [Wolter 07] and [Rust, Kalton 87])
    provided in function \code{\link{collapse.strata}}.
    As a possible alternative, you can get rid of lonely PSUs also by setting
    proper variance estimation options via \code{options("RG.lonely.psu")}.
    The default setting is \code{"fail"}, which raises an error if a lonely PSU
    is met. Option \code{"remove"} simply causes the software to ignore lonely PSUs
    for variance computation purposes. Option \code{"adjust"} means that
    deviations from the \emph{population mean} will be used in variance
    estimation formulae, instead of deviations from the stratum mean
    (a conservative choice). Finally, option \code{"average"} causes the
    software to replace the variance contribution of the stratum by the average
    variance contribution across strata (this can be appropriate e.g. when one
    believes that lonely PSU strata occur at random due to uniform nonresponse
    among strata).
    
    The variance formulae for domain estimation give well-defined,
    positive results when a stratum contains only a single PSU with
    observations falling in the domain, but are not unbiased.\cr
    If \code{options("RG.adjust.domain.lonely")} is \code{TRUE}
    and \code{options("RG.lonely.psu")} is \code{"average"} or
    \code{"adjust"} the same adjustment for lonely PSUs will be used
    within a domain. Note that this adjustment is not available for
    calibrated designs.
    
    If \code{options("RG.warn.domain.lonely")} is set to \code{TRUE}, a
    warning message is raised whenever an estimation domain happens to
    contain just a single PSU belonging to a stratum. The default is \code{FALSE}.
    }

   \references{
    Kalton, G. (1979). \emph{\dQuote{Ultimate cluster sampling}}, Journal of the Royal Statistical Society, Series A, 142, pp. 210-222.

    Bellhouse, D. R. (1985). \emph{\dQuote{Computing Methods for Variance Estimation in Complex Surveys}}. Journal of Official Statistics, Vol. 1, No. 3, pp. 323-329.

    Lumley, T. (2006) \emph{\dQuote{survey: analysis of complex survey samples}}, \url{https://CRAN.R-project.org/package=survey}.

    Wolter, K.M. (2007) \emph{\dQuote{Introduction to Variance Estimation}}, Second Edition, Springer-Verlag, New York.

    Rust, K., Kalton, G. (1987) \emph{\dQuote{Strategies for Collapsing Strata for Variance Estimation}}, Journal of Official Statistics, Vol. 3, No. 1, pp. 69-81.
    }

  \seealso{
    \code{\link{e.svydesign}} and its \code{self.rep.str} argument for a 
    "compromise solution" that can be adopted when the sampling design
    involves self-representing (SR) strata, \code{\link{collapse.strata}}
    for the suggested way of handling lonely PSUs, and \code{\link{fpcdat}}
    for useful data examples.
    }

  \examples{
# Define a two-stage stratified cluster sampling without
# replacement:
data(fpcdat)
des<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
     fpc=~fpc1+fpc2)

# Now compare SE (or CV%) sizes under different settings:

## 1) Default setting, i.e. Ultimate Cluster Approximation is off
svystatTM(des,~x+y+z,vartype=c("se","cvpct"))

## 2) Turn on the Ultimate Cluster Approximation, thus missing
##    the variance contribution from the second stage
##    (hence SR strata give no contribution at all):
old.op <- options("RG.ultimate.cluster"=TRUE)
svystatTM(des,~x+y+z,vartype=c("se","cvpct"))
options(old.op)

## 3) The "compromise solution" (see ?e.svydesign) i.e. retaining
##    only the leading contribution to the sampling variance (namely
##    the one arising from SSUs in SR strata and PSUs in not-SR strata):
des2<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,weights=~w,
      fpc=~fpc1+fpc2, self.rep.str=~sr)
svystatTM(des2,~x+y+z,vartype=c("se","cvpct"))

# Therefore, sampling variances come out in the expected
# hierarchy: 1) > 3) > 2).


# Under default settings lonely PSUs produce errors in standard
# errors estimation (notice we didn't pass the fpcs):
data(fpcdat)
des.lpsu<-e.svydesign(data=fpcdat,ids=~psu+ssu,strata=~stratum,
          weights=~w)
\dontrun{
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
}

# This can be circumvented in different ways, namely:
old.op <- options("RG.lonely.psu"="adjust")
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
options(old.op)

# or:
old.op <- options("RG.lonely.psu"="average")
svystatTM(des.lpsu,~x+y+z,vartype=c("se","cvpct"))
options(old.op)

# or otherwise by collapsing strata inside planned
# estimation domains:
des.clps<-collapse.strata(design=des.lpsu,block.vars=~pl.domain)
svystatTM(des.clps,~x+y+z,vartype=c("se","cvpct"))
}
  \keyword{survey}