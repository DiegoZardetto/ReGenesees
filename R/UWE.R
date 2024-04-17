# References
# - [Kish 92]
# Kish, L. (1992). Weighting for unequal Pi. Journal of Official Statistics, 8, 183-200.

`UWE` <- function(design, by = NULL) {
############################################################################
# Compute the UWE for the current (w) and initial (w0) weights of a design # 
# object, plus the corresponding VAR.infl factor induced by w0 -> w.       #
#                                                                          #
# NOTE: Adopts Kish original definition of  UWE, which used the 1/n        #
#       version of the sample variance.                                    #
#                                                                          #
# NOTE: Initial weights w0 are:                                            #
#       - equal to w if design is an initial design object                 #
#       - equal to the *starting weights* (i.e. the reciprocals of the     #
#         $allprob slot of design) if design was obtained by the           #
#         application of an arbitrary *chain* of ReGenesees functions that #
#         modify the weights (e.g. smooth.strat.jump, e.calibrate,         #
#         ext.calibrated, trimcal, ...).                                   #
#                                                                          #
# NOTE: Levels of by variables follow the same ordering that svyby would   #
#       generate.                                                          #
############################################################################
# Check on design
if (!inherits(design, "analytic")) {
     stop("Object 'design' must inherit from class analytic")
    }

  w <- as.numeric(weights(design))
  w0 <- 1 / design$allprob[[1]]

  if (is.null(by)) {
     out <- UWE.all(w, w0)
    } else {
     if (!inherits(by, "formula")) {
             stop("If specified, 'by' must be supplied as a formula")
        }

     # Check on NAs in by variables
     na.Fail(design$variables, all.vars(by))

     ## All combinations that actually occur in this design
     byfactors <- model.frame(by, model.frame(design), na.action = na.pass)

     ## Order combinations as svyby would do (this requires the rev() calls below)
     byfactors.s <- byfactors[do.call(order, rev(byfactors)), , drop = FALSE]
     w <- w[do.call(order, rev(byfactors))]
     w0 <- w0[do.call(order, rev(byfactors))]

     byfactor <- do.call("interaction", byfactors.s)
     uniquelevels <- unique(byfactor)
     u.byfactors.s <- unique(byfactors.s)

#     uniques <- match(uniquelevels, byfactor)

     out <- sapply(uniquelevels,
                     function(lev){
                         UWE.all( w[byfactor == lev], w0[byfactor == lev])
                        }
                    )

     out <- cbind(u.byfactors.s, t(out))
     out <- data.frame(lapply(out, unlist)) # Debug 17/04/2024: same output, just better
                                            #                   subsetting behavior
     rownames(out) <- uniquelevels
    }

  return(out)
}

`UWE.all` <- function(w, w0) {
#################################################
# Compute UWE and for w and w0 and the VAR.infl #
# factor induced by w0 -> w.                    #
#                                               #
# NOTE: Kish original definition, which used    #
#       the 1/n version of the sample variance. #
#################################################

  if (!identical(length(w), length(w0))){
         stop("Current and initial weights length differ!")
    }

  n <- length(w)

  w.avg <- mean(w)
  w.var <- ((n - 1) / n) * var(w)
  w.UWE <- 1 + w.var / (w.avg^2)

  w0.avg <- mean(w0)
  w0.var <- ((n - 1) / n) * var(w0)
  w0.UWE <- 1 + w0.var / (w0.avg^2)

  out <- data.frame("UWE.curr" = w.UWE, "UWE.ini" = w0.UWE, "VAR.infl" = w.UWE / w0.UWE)
  rownames(out) <- "ALL"
  return(out)
}


`G.wadj` <- function(design, eps = 1e-12) {
###################################################################
# Get the current (w) and initial (w0) weights of a design object #
# that was generated in a weights changing pipeline, and compute  #
# the overall pipeline-level weight-adjustments G = w / w0.       #
#                                                                 #
# NOTE: Initial weights w0 are:                                   #
#       - equal to w if design is an initial design object        #
#       - equal the *starting weights* (i.e. the reciprocals of   #
#         the $allprob slot of design) if design was obtained by  #
#         the application of an arbitrary *chain* of ReGenesees   #
#         functions that modify the weights (e.g.                 #
#         smooth.strat.jump, e.calibrate, ext.calibrated,         #
#         trimcal, ...).                                          #
#                                                                 #
# NOTE: Value eps is a small cutoff to avoid NaN results from     #
#       ratios w / w0 that would be equal to 0 / 0.               #
###################################################################
  w <- as.numeric(weights(design))
  is.w.eps <- (abs(w) < eps)
  ww.eps <- w
  ww.eps[is.w.eps] <- eps * ifelse(sign(ww.eps[is.w.eps]) == 0, 1, sign(ww.eps[is.w.eps]))

  w0 <- 1 / design$allprob[[1]]
  is.w0.eps <- (abs(w0) < eps)
  ww0.eps <- w0
  ww0.eps[is.w0.eps] <- eps * ifelse(sign(ww0.eps[is.w0.eps]) == 0, 1, sign(ww0.eps[is.w0.eps]))

  G <- w / w0
  is.G.eps <- (is.w.eps & is.w0.eps)
  G[is.G.eps] <- (ww.eps / ww0.eps)[is.G.eps]
  G
 }
 
 

`is.selfweighted` <- function(design, tol = 1E-6, ...) {
     w <- weights(design)
     sw.dist <- ( abs(diff(range(w))) / abs(mean(w)) )
     if ( sw.dist > tol) {
         ans <- FALSE
        } else {
         ans <- TRUE
        }
     attr(ans, "sw.dist") <- sw.dist
     return(ans)
    }
