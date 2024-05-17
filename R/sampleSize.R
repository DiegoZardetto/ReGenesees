# Author: Diego Zardetto, 2023.
# I'm considering including some simple, everyday tools for basic sampling
# design decisions and survey planning into ReGenesees. Below are a couple of
# initial code snippets.
#
# NOTE: First release ReGenesees 2.3. Treatment of fpc added in ReGenesees 2.4. 
#
# Possible reference URLs for documentation/testing/validation:
# - MICS7 TOOLS, UNICEF 2023. URL: https://mics.unicef.org/tools?round=mics7
# - ILO-IPEC TOOLS, ILO 2014. URL: https://www.ilo.org/ipec/ChildlabourstatisticsSIMPOC/Manuals/WCMS_304559/lang--en/index.htm

# For hypothesis testing comparing two proportions (n.comp2prop)
# - Rosner, B. (2006) Fundamentals of biostatistics (7th ed.). Boston, MA: Brooks/Cole.
# - Equation 10.14 on page 381

# For hypothesis testing comparing two means (n.comp2mean)
# - Rosner, B. (2006) Fundamentals of biostatistics (7th ed.). Boston, MA: Brooks/Cole.
# - Equation 8.27 on page 302


n.prop <- function(prec, prec.ind = c("ME", "RME", "SE", "CV"), P = 0.5,
                   DEFF = 1, RR = 1,
                   F = 1, hhSize = 1, AVEhh = F * hhSize,
                   old.clus.size = NULL, new.clus.size = NULL,
                   N = NULL, alpha = 0.05, verbose = TRUE){
################################################################
# Sample size requirements for specified precision constraints #
# in the estimation of a PROPORTION.                           #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
Z <- qnorm(1 - alpha/2)
# Handle the precision indicator used to set the precision constraint to be
# satisfied
prec.ind <- match.arg(prec.ind)
## ME  = Z * SE
## Use ME as a pivot
# RME = ME / P                          -> ME = RME * P
# SE = ME / Z                           -> ME = SE * Z 
# CV  = SE / P = ME / (Z * P) = RME / Z -> ME = CV * (Z * P)
ME <- switch(prec.ind, ME = prec, RME = prec * P, SE = prec * Z, CV = prec * (Z * P))

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Calculate the required sample size
n.inf <- ( Z^2 * DEFF * P * (1 - P) ) / (ME^2 * RR * AVEhh)

# Address fpc
if (is.null(N)) {
     n <- ceiling(n.inf)
    } else {
     # Make sure input is a sensical number
     N <- round(abs(N))
     n <- ceiling( n.inf / ( 1 + (n.inf - 1) / N ) )
    }

if (!is.null(new.clus.size)) {
     # Make the required sample size a multiple of the requested new cluster
     # size
     n <- new.clus.size * (n %/% new.clus.size + sign(n %% new.clus.size))
    }

if (isTRUE(verbose)) {
     cat("# Precision constraint:")
     cat(paste("\n  ", prec.ind, sep = ""), "=", prec)
     if (prec.ind %in% c("ME", "RME")) {
          cat("\n  alpha =", alpha)
         }
     cat("\n# Anticipated estimates:")
     cat("\n  P =", P)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(N)) {
          cat("\n  N =", N)
         }
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Required sample size:")
     cat("\n  n =", n, "\n\n")

     return(invisible(n))
    }
return(n)
}

prec.prop <- function(n, prec.ind = c("ME", "RME", "SE", "CV"), P = 0.5,
                      DEFF = 1, RR = 1,
                      F = 1, hhSize = 1, AVEhh = F * hhSize,
                      old.clus.size = NULL, new.clus.size = NULL,
                      N = NULL, alpha = 0.05, verbose = TRUE){
################################################################
# Expected precision in the estimation of a PROPORTION for a   #
# specified sample size.                                       #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# precision at fixed sample size.                              #
################################################################
Z <- qnorm(1 - alpha/2)
# Handle the precision indicator that the user wants to predict
prec.ind <- match.arg(prec.ind)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the expected ME
ME.inf <- sqrt( ( Z^2 * DEFF * P * (1 - P) ) / (n * RR * AVEhh) )

# Address fpc
if (is.null(N)) {
     ME <- ME.inf
    } else {
     # Make sure input is a sensical number
     N <- round(abs(N))
     if ( any(n > N) ) stop("Specified sample size (n) is larger than specified population size (N), please double-check!")
     ME <- ME.inf * sqrt( (N - n) / (N - 1) )
    }

## Use ME as a pivot
# ME  = Z * SE
# RME = ME / P
# SE = ME / Z
# CV  = SE / P = ME / (Z * P) = RME / Z
prec <- switch(prec.ind, ME = ME, RME = ME / P, SE = ME / Z, CV = ME / (Z * P))

if (isTRUE(verbose)) {
     cat("# Sample size:")
     cat("\n  n =", n)
     cat("\n# Anticipated estimates:")
     cat("\n  P =", P)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(N)) {
          cat("\n  N =", N)
         }
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     if (prec.ind %in% c("ME", "RME")) {
          cat("\n# Significance:")
          cat("\n  alpha =", alpha)
         }
     cat("\n# -> Expected precision:")
     cat(paste("\n  ", prec.ind, sep = ""), "=", prec)
     cat("\n\n")

     return(invisible(prec))
    }
return(prec)
}

n.mean <- function(prec, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
                   DEFF = 1, RR = 1,
                   F = 1, hhSize = 1, AVEhh = F * hhSize,
                   old.clus.size = NULL, new.clus.size = NULL,
                   N = NULL, alpha = 0.05, verbose = TRUE){
################################################################
# Sample size requirements for specified precision constraints #
# in the estimation of the MEAN of a numeric variable Y.       #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
Z <- qnorm(1 - alpha/2)
# Handle the precision indicator used to set the precision constraint to be
# satisfied
prec.ind <- match.arg(prec.ind)
## ME  = Z * SE
## Use ME as a pivot
# RME = ME / muY                            -> ME = RME * muY
# SE = ME / Z                               -> ME = SE * Z 
# CV  = SE / muY = ME / (Z * muY) = RME / Z -> ME = CV * (Z * muY)
ME <- switch(prec.ind, ME = prec, RME = prec * muY, SE = prec * Z, CV = prec * (Z * muY))

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Calculate the required sample size
sigma2Y <- sigmaY^2
n.inf <- ( Z^2 * DEFF * sigma2Y ) / (ME^2 * RR * AVEhh)

# Address fpc
if (is.null(N)) {
     n <- ceiling(n.inf)
    } else {
     # Make sure input is a sensical number
     N <- round(abs(N))
     n <- ceiling( n.inf / ( 1 + (n.inf - 1) / N ) )
    }

if (!is.null(new.clus.size)) {
     # Make the required sample size a multiple of the requested new cluster
     # size
     n <- new.clus.size * (n %/% new.clus.size + sign(n %% new.clus.size))
    }

if (isTRUE(verbose)) {
     cat("# Precision constraint:")
     cat(paste("\n  ", prec.ind, sep = ""), "=", prec)
     if (prec.ind %in% c("ME", "RME")) {
          cat("\n  alpha =", alpha)
         }
     cat("\n# Anticipated estimates:")
     cat("\n  muY =", muY)
     cat("\n  sigmaY =", sigmaY)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(N)) {
          cat("\n  N =", N)
         }
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Required sample size:")
     cat("\n  n =", n, "\n\n")

     return(invisible(n))
    }
return(n)
}

prec.mean <- function(n, prec.ind = c("ME", "RME", "SE", "CV"), sigmaY, muY = NULL,
                      DEFF = 1, RR = 1,
                      F = 1, hhSize = 1, AVEhh = F * hhSize,
                      old.clus.size = NULL, new.clus.size = NULL,
                      N = NULL, alpha = 0.05, verbose = TRUE){
################################################################
# Expected precision in the estimation of the MEAN of a        #
# numeric variable Y for a specified sample size.              #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# precision at fixed sample size.                              #
################################################################
Z <- qnorm(1 - alpha/2)
# Handle the precision indicator used to set the precision constraint to be
# satisfied
prec.ind <- match.arg(prec.ind)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the the expected ME
sigma2Y <- sigmaY^2
ME.inf <- sqrt( ( Z^2 * DEFF * sigma2Y ) / (n * RR * AVEhh) )

# Address fpc
if (is.null(N)) {
     ME <- ME.inf
    } else {
     # Make sure input is a sensical number
     N <- round(abs(N))
     if ( any(n > N) ) stop("Specified sample size (n) is larger than specified population size (N), please double-check!")
     ME <- ME.inf * sqrt( (N - n) / (N - 1) )
    }

## Use ME as a pivot
# ME  = Z * SE
# RME = ME / muY
# SE = ME / Z
# CV  = SE / muY = ME / (Z * muY) = RME / Z
prec <- switch(prec.ind, ME = ME, RME = ME / muY, SE = ME / Z, CV = ME / (Z * muY))

if (isTRUE(verbose)) {
     cat("# Sample size:")
     cat("\n  n =", n)
     cat("\n# Anticipated estimates:")
     cat("\n  muY =", muY)
     cat("\n  sigmaY =", sigmaY)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     if (prec.ind %in% c("ME", "RME")) {
          cat("\n# Significance:")
          cat("\n  alpha =", alpha)
         }
     cat("\n# -> Expected precision:")
     cat(paste("\n  ", prec.ind, sep = ""), "=", prec)
     cat("\n\n")

     return(invisible(prec))
    }
return(prec)
}


n.comp2prop <- function(P1, P2 = P1, MDE = abs(P2 - P1), K1 = 1/2,
                        alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
                        pooled.variance = TRUE,
                        DEFF = 1, RR = 1,
                        F = 1, hhSize = 1, AVEhh = F * hhSize,
                        old.clus.size = NULL, new.clus.size = NULL,
                        verbose = TRUE){
################################################################
# Sample size requirements for attaining specified levels of   #
# significance and power in a statistical test that COMPARES   #
# TWO PROPORTIONS with the aim of detecting as statistically   #
# significant a specified level of difference that might exist #
# between the proportions (MDE, Minimum Detectable Effect).    #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))
Z.beta <- qnorm(1 - beta)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Calculate the required sample size
if (isTRUE(pooled.variance)) {
     Pbar <- P1 * K1 + P2 * (1 - K1)
     n1 <- ( ( ( Z.alpha * sqrt( Pbar * (1 - Pbar) / ( 1 - K1 ) ) +
                 Z.beta  * sqrt( P1 * ( 1 - P1 ) + P2 * ( 1 - P2 ) * ( K1 / ( 1 - K1 ) ) )
                )^2 
              ) * ( DEFF / (MDE^2 * RR * AVEhh) )
            )
     n2 <- ( ( 1 - K1 ) / K1 ) * n1

     n1 <- ceiling(n1)
     n2 <- ceiling(n2)
     n <- n1 + n2
    } else {
     n <- ( ( Z.alpha + Z.beta )^2 ) * ( P1 * ( 1 - P1 ) / K1 + P2 * ( 1 - P2 ) / ( 1 - K1 ) ) *
            ( DEFF / (MDE^2 * RR * AVEhh) )
     n1 <- ceiling( K1 * n )
     n2 <- ceiling( ( 1 - K1 ) * n )
     n <- n1 + n2
    }

if (!is.null(new.clus.size)) {
     # Make the required sample size a multiple of the requested new cluster
     # size
     n1 <- new.clus.size * (n1 %/% new.clus.size + sign(n1 %% new.clus.size))
     n2 <- new.clus.size * (n2 %/% new.clus.size + sign(n2 %% new.clus.size))
     n <- n1 + n2 
    }

out <- list(n1 = n1, n2 = n2, n = n)

if (isTRUE(verbose)) {
     cat("# Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Power:")
     cat("\n  1 - beta", "=", 1 - beta)
     cat("\n# Anticipated estimates:")
     cat("\n  P1 =", P1)
     if (!missing(P2)) {
          cat("\n  P2 =", P2)
         }
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     cat("\n# Design parameters:")
     cat("\n  K1 =", K1)
     if (!is.null(new.clus.size)) {
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Required sample size:")
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  n =", n, "\n\n")

     return(invisible(out))
    }
return(out)
}

pow.comp2prop <- function(n, P1, P2 = P1, MDE = abs(P2 - P1), K1 = 1/2,
                          alpha = 0.05, sides = c("two-tailed", "one-tailed"),
                          pooled.variance = TRUE,
                          DEFF = 1, RR = 1,
                          F = 1, hhSize = 1, AVEhh = F * hhSize,
                          old.clus.size = NULL, new.clus.size = NULL,
                          verbose = TRUE){
################################################################
# Expected power for a statistical test that COMPARES TWO      #
# PROPORTIONS given specified sample size, significance, and   #
# MDE (Minimum Detectable Effect).                             #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the expected power
if (isTRUE(pooled.variance)) {
     Pbar <- P1 * K1 + P2 * (1 - K1)
     A <- ( DEFF / (MDE^2 * RR * AVEhh) )
     R1 <- sqrt( Pbar * (1 - Pbar) / ( 1 - K1 ) )
     R2 <- sqrt( P1 * ( 1 - P1 ) + P2 * ( 1 - P2 ) * ( K1 / ( 1 - K1 ) ) )
     Power <- pnorm( ( sqrt( ( K1 * n ) / A ) - Z.alpha * R1 ) / R2 )
     n1 <- round(n * K1)
     n2 <- n - n1
    } else {
     A <- ( DEFF / (MDE^2 * RR * AVEhh) )
     B <- ( P1 * ( 1 - P1 ) / K1 + P2 * ( 1 - P2 ) / ( 1 - K1 ) )
     Power <- pnorm( sqrt( n / (B*A) ) - Z.alpha )
     n1 <- round(n * K1)
     n2 <- n - n1
    }

if (isTRUE(verbose)) {
     cat("# Sample size(s):")
     cat("\n  n =", n)
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  K1 =", K1)
     cat("\n# Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Anticipated estimates:")
     cat("\n  P1 =", P1)
     if (!missing(P2)) {
          cat("\n  P2 =", P2)
         }
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Expected Power:")
     cat("\n  1 - beta", "=", round(Power, 3), "\n\n")

     return(invisible(Power))
    }
return(Power)
}

mde.comp2prop <- function(n, P1, P2 = P1, K1 = 1/2,
                          alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
                          pooled.variance = TRUE,
                          DEFF = 1, RR = 1,
                          F = 1, hhSize = 1, AVEhh = F * hhSize,
                          old.clus.size = NULL, new.clus.size = NULL,
                          verbose = TRUE){
################################################################
# Expected MDE for a statistical test that COMPARES TWO        #
# PROPORTIONS given specified sample size, significance, and   #
# power.                                                       #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))
Z.beta <- qnorm(1 - beta)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the expected MDE
if (isTRUE(pooled.variance)) {
     Pbar <- P1 * K1 + P2 * (1 - K1)
     MDE <- sqrt( ( ( Z.alpha * sqrt( Pbar * (1 - Pbar) / ( 1 - K1 ) ) +
                      Z.beta  * sqrt( P1 * ( 1 - P1 ) + P2 * ( 1 - P2 ) * ( K1 / ( 1 - K1 ) ) )
                    )^2 
                  ) * ( DEFF / (n * K1 * RR * AVEhh) )
                )
    } else {
     MDE <- sqrt( ( ( Z.alpha + Z.beta )^2 ) * ( P1 * ( 1 - P1 ) / K1 + P2 * ( 1 - P2 ) / ( 1 - K1 ) ) *
                  ( DEFF / (n * RR * AVEhh) ) )
    }
n1 <- round(n * K1)
n2 <- n - n1

if (isTRUE(verbose)) {
     cat("# Sample size(s):")
     cat("\n  n =", n)
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  K1 =", K1)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Power:")
     cat("\n  1 - beta", "=", 1 - beta)
     cat("\n# Anticipated estimates:")
     cat("\n  P1 =", P1)
     if (!missing(P2)) {
          cat("\n  P2 =", P2)
         }
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Expected Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE, "\n\n")

     return(invisible(MDE))
    }
return(MDE)
}


n.comp2mean <- function(sigmaY, MDE, K1 = 1/2,
                        alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
                        DEFF = 1, RR = 1,
                        F = 1, hhSize = 1, AVEhh = F * hhSize,
                        old.clus.size = NULL, new.clus.size = NULL,
                        verbose = TRUE){
################################################################
# Sample size requirements for attaining specified levels of   #
# significance and power in a statistical test that COMPARES   #
# TWO MEANS with the aim of detecting as statistically         #
# significant a specified level of difference that might exist #
# between the means (MDE, Minimum Detectable Effect).          #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))
Z.beta <- qnorm(1 - beta)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Calculate the required sample size
sigma2Y <- sigmaY^2
n1 <- ( ( Z.alpha + Z.beta )^2 ) * ( sigma2Y / ( 1 - K1 ) ) * ( DEFF / (MDE^2 * RR * AVEhh) )
n2 <- ( ( 1 - K1 ) / K1 ) * n1

n1 <- ceiling(n1)
n2 <- ceiling(n2)
n <- n1 + n2

if (!is.null(new.clus.size)) {
     # Make the required sample size a multiple of the requested new cluster
     # size
     n1 <- new.clus.size * (n1 %/% new.clus.size + sign(n1 %% new.clus.size))
     n2 <- new.clus.size * (n2 %/% new.clus.size + sign(n2 %% new.clus.size))
     n <- n1 + n2
    }

out <- list(n1 = n1, n2 = n2, n = n)

if (isTRUE(verbose)) {
     cat("# Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Power:")
     cat("\n  1 - beta", "=", 1 - beta)
     cat("\n# Anticipated estimates:")
     cat("\n  sigmaY =", sigmaY)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     cat("\n# Design parameters:")
     cat("\n  K1 =", K1)
     if (!is.null(new.clus.size)) {
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Required sample size:")
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  n =", n, "\n\n")

     return(invisible(out))
    }
return(out)
}

pow.comp2mean <- function(n, sigmaY, MDE, K1 = 1/2,
                          alpha = 0.05, sides = c("two-tailed", "one-tailed"),
                          DEFF = 1, RR = 1,
                          F = 1, hhSize = 1, AVEhh = F * hhSize,
                          old.clus.size = NULL, new.clus.size = NULL,
                          verbose = TRUE){
################################################################
# Expected power for a statistical test that COMPARES TWO      #
# MEANS given specified sample size, significance, and MDE     #
# (Minimum Detectable Effect).                                 #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the expected power
sigma2Y <- sigmaY^2
A <- ( DEFF / (MDE^2 * RR * AVEhh) )
B <- ( sigma2Y / ( 1 - K1 ) )
Power <- pnorm( sqrt( ( K1 * n ) / (B*A) ) - Z.alpha )
n1 <- round(n * K1)
n2 <- n - n1

if (isTRUE(verbose)) {
     cat("# Sample size(s):")
     cat("\n  n =", n)
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  K1 =", K1)
     cat("\n# Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Anticipated estimates:")
     cat("\n  sigmaY =", sigmaY)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Expected Power:")
     cat("\n  1 - beta", "=", round(Power, 3), "\n\n")

     return(invisible(Power))
    }
return(Power)
}

mde.comp2mean <- function(n, sigmaY, K1 = 1/2,
                          alpha = 0.05, beta = 0.2, sides = c("two-tailed", "one-tailed"),
                          DEFF = 1, RR = 1,
                          F = 1, hhSize = 1, AVEhh = F * hhSize,
                          old.clus.size = NULL, new.clus.size = NULL,
                          verbose = TRUE){
################################################################
# Expected MDE for a statistical test that COMPARES TWO        #
# MEANS given specified sample size, significance, and power.  #
#                                                              #
# Also exploits the ICC *portability* assumption to guess the  #
# impact that a change of cluster size would have on the       #
# required sample size.                                        #
################################################################
sides <- match.arg(sides)
Z.alpha <- switch(sides, "two-tailed" = qnorm(1 - alpha/2), "one-tailed" = qnorm(1 - alpha))
Z.beta <- qnorm(1 - beta)

## Arguments F, hhSizehave, and AVEhh are related: which ones have been passed?
## If possible, deduce values that were not explicitly passed.
if (!missing(AVEhh)) {
     if (!missing(hhSize) & !missing(F)) {
         if (!isTRUE(all.equal(AVEhh, F * hhSize))) stop("Specified arguments do not comply with AVEhh = F * hhSize, please double-check!")
        }
     if (missing(F) & !missing(hhSize)) {
         F <- AVEhh / hhSize
        }
     if (missing(hhSize) & !missing(F)) {
         hhSize <- AVEhh / F
        }
     if (missing(hhSize) & missing(F)) {
         hhSize <- NULL
         F <- NULL
        }
    }

# Guesstimate the impact of changing cluster size (i.e. number of sampled hh per
# sampled PSU) by leveraging DEFF portability
DEFF.old <- DEFF
if (!is.null(new.clus.size)) {
     # Make sure input is a sensical number
     new.clus.size <- round(abs(new.clus.size))
     if (is.null(old.clus.size)) stop("Specify the cluster size that resulted in the input DEFF value!")
     # Make sure input is a sensical number
     old.clus.size <- round(abs(old.clus.size))
     # NOTE: The ROH estimated here is an individual-level one, because the
     #       the average cluster size is calculated in terms of individuals, not
     #       households! 
     ROH.old <- (DEFF.old - 1) / (old.clus.size * AVEhh - 1)
     DEFF <- 1 + ROH.old * (new.clus.size * AVEhh - 1)
    }

# Make sure input is a sensical number
n <- round(abs(n))

# Calculate the expected MDE
sigma2Y <- sigmaY^2
MDE <- sqrt( ( ( Z.alpha + Z.beta )^2 ) * ( sigma2Y / ( 1 - K1 ) ) * ( DEFF / (n * K1 * RR * AVEhh) ) )
n1 <- round(n * K1)
n2 <- n - n1

if (isTRUE(verbose)) {
     cat("# Sample size(s):")
     cat("\n  n =", n)
     cat("\n  n1 =", n1)
     cat("\n  n2 =", n2)
     cat("\n  K1 =", K1)
     cat("\n# Significance:")
     cat("\n  alpha", "=", alpha)
     cat("\n# Power:")
     cat("\n  1 - beta", "=", 1 - beta)
     cat("\n# Anticipated estimates:")
     cat("\n  sigmaY =", sigmaY)
     cat("\n  DEFF =", DEFF.old)
     cat("\n  RR =", RR)
     cat("\n  F =", F)
     cat("\n  hhSize =", hhSize)
     cat("\n  AVEhh =", AVEhh)
     if (!is.null(new.clus.size)) {
          cat("\n# Design parameters:")
          cat("\n  old.clus.size =", old.clus.size)
          cat("\n  ROH =", ROH.old)
          cat("\n  new.clus.size =", new.clus.size)
          cat("\n  new.DEFF =", DEFF)
         }
     cat("\n# -> Expected Minimum Detectable Effect:")
     cat("\n  MDE", "=", MDE, "\n\n")

     return(invisible(MDE))
    }
return(MDE)
}
