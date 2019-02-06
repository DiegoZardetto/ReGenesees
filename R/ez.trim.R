`ez.trim` <-
function (design, formula, population, 
          aggregate.stage = NULL, sigma2 = NULL, weights.back = NULL, 
          w.range = c(-Inf, Inf), trimfun = "linear",
          maxit = 50, epsilon = 1e-07, force = FALSE)
############################################
# ---------------------------------------- #
#    VARIANCE ESTIMATION WORKS THIS WAY    #
# ---------------------------------------- #
#                                          #
#      CAL          TRIM                   #
#   w -----> w.cal ------> w.cal.cal <--   #
#   |  g1           g2                 |   #
#   |  s2           s2                 |   #
#    ----------------------------------    #
#                CAL(+)TRIM                #
#               g.tld = g1*g2              #
#      s2.tld = g.tld*s2 = g1*g2*s2        #
#                                          #
############################################
{
    # NOTE: ONLY linear trimming metric enabled, TO DATE.
    #       This preference is essentially dictated by numerical stability
    #       considerations.
    if (is.character(trimfun)) 
        trimfun <- match.arg(trimfun)
    if (is.character(trimfun)) 
        trimfun <- switch(trimfun, linear = trim.linear, raking = trim.raking, 
            logit = trim.logit)
    else if (!inherits(trimfun, "trimfun"))
        stop("'trimfun' must be a string or of class 'trimfun'.")

    mm <- model.matrix(formula, model.frame(formula, data = design$variables))
    ww <- design$dir.weights
    # Heteroskedsticity
    if (is.null(sigma2)) {
        sigma2 <- rep(1, NROW(mm))
        }
    else {
        sigma2 <- design$variables[, all.vars(sigma2)]
        }

      # No whalf needed for calibration trimming!
      # whalf <- sqrt(ww)
    sigma <- sqrt(sigma2)
      # No QR decomposition needed for calibration trimming!
      # tqr <- qr(mm * (whalf/sigma))
    sample.total <- colSums(mm * ww)
    if (length(sample.total) != length(population)) 
        stop("Population and sample totals are not the same length.")
    if (any(sample.total == 0)) {
#       zz <- (population == 0) & (apply(mm, 2, function(x) all(x == 0))) # OLD: It's much slower than NEW
        zz <- (population == 0) & (colSums(abs(mm)) == 0)
        mm <- mm[, !zz, drop = FALSE]
        population <- population[!zz]
        sample.total <- sample.total[!zz]
    }
    # Cluster level model:
      # OLD: It's much slower than NEW
#    if (!is.null(aggregate.stage)) {
#        aggindex <- design$cluster[[aggregate.stage]]
#        mm <- apply(mm, 2, function(col) tapply(col, aggindex, sum))
#        ww <- as.numeric(tapply(ww, aggindex, mean))
#        sigma2 <- as.numeric(tapply(sigma2, aggindex, mean))
#    }
      # NEW: rowsum is definitely faster
    if (!is.null(aggregate.stage)) {
        aggindex <- design$cluster[[aggregate.stage]]
        n.aggindex <- rowsum(rep(1, NROW(mm)), aggindex)
        mm <- rowsum(mm, aggindex)
        ww <- as.numeric(rowsum(ww, aggindex) / n.aggindex)
        sigma2 <- as.numeric(rowsum(sigma2, aggindex) / n.aggindex)
    }

    ## The key point is here: compute unit-specific ranges for the g-weights,
    ## taking into account the bounds imposed on calibration weights
       #------------------------------------------------------------------------
       # The OLD WAY: g-weights range NO LONGER considered
       # bounds <- cbind(pmax(w.range[1]/ww, bounds[1]), pmin(w.range[2]/ww, bounds[2]))
       # DEBUG 10/01/2017: code line below did not correctly cope with negative ww
       # bounds <- cbind(w.range[1]/ww, w.range[2]/ww)
       #------------------------------------------------------------------------
    ## NOTE: The effective unit-specific bounds for the g-weights must ensure
    ##       that the product g1*g2 is positive, i.e. g1 and g2 must have the
    ##       same sign. If both limits specified by w.range are positive, the
    ##       rules below do the right job (even if w.cal has negative entries).
    boundsL <- ifelse(ww > 0, w.range[1]/ww, ifelse(ww <0, w.range[2]/ww, -Inf))
    boundsU <- ifelse(ww > 0, w.range[2]/ww, ifelse(ww <0, w.range[1]/ww,  Inf))
    bounds <- cbind(boundsL, boundsU)
    ## DONE

    g <- ez.grake(sample.total, mm, ww, trimfun, bounds = bounds, population = population,
                  epsilon = epsilon, maxit = maxit, sigma2 = sigma2)
    if (!force && !is.null(attr(g, "failed"))) 
        stop("Calibration failed")
    ret.code <- if (is.null(attr(g, "failed"))) 0 else 1
    fail.diagnostics <- attr(g, "fail.diagnostics")

    # For aggregated calibration must re-expand g-weights
    if (!is.null(aggregate.stage)) {
        g <- g[aggindex]
        sigma2 <- sigma2[aggindex]  # DEBUG 19/02/2017: Forgot to re-expand sigma2, fixed.
        ww <- ww[aggindex]          # DEBUG 19/02/2017: Same as above, fixed.
        mm <- mm[aggindex, ]        # DEBUG 19/02/2017: Same as above, fixed.
    }

    # Handle zero calibrated weights (if any): must provide a tiny cutoff
    # to prevent problems in svyrecvar (fictitious 0/0=NaN which causes an error
    # in qr). Results DO NOT DEPEND on the choosen cutoff value.
    # NOTE (11/07/2012): The same problem would manifest itself if some of the
    #                    direct weights is 0. This should not happen in
    #                    production (as a 0 direct weight doesn't make sense),
    #                    but can happen in research (e.g. when experimenting on
    #                    calibrated Jackknife).
    g.0 <- (abs(g) < 1E-12)
    if (any(g.0)) g[g.0] <- (1E-12) * ifelse(sign(g[g.0])==0, sign(ww[g.0]), sign(g[g.0]))   # DEBUG 14/02/2017: Those g which equal 0 STRICTLY
                                                                                             #                   must be cutoff PRESERVING the
                                                                                             #                   original sign of the corresponding ww!
                                                                                             # NOTE: This is only needed for trimming, due to how variance
                                                                                             #       estimation goes (see below)

    cal.weights <- g * design$dir.weights
    attr(cal.weights,"ret.code") <- ret.code
    attr(attr(cal.weights,"ret.code"), "fail.diagnostics") <- fail.diagnostics
      # No QR decomposition needed for calibration trimming!
      # attr(cal.weights,"qr") <- tqr
      # No gwhalf needed (enters variance estimation) for calibration trimming!
      # attr(cal.weights,"gwhalf") <- g * whalf * sigma

# ---------------------------------------- #
#    VARIANCE ESTIMATION WORKS THIS WAY    #
# ---------------------------------------- #
      weights.back <- design$variables[, all.vars(weights.back)]
      # g1 <- ww/weights.back; sigma2.tld <- g1*g*sigma2; hence:
      sigma2.tld <- (cal.weights/weights.back)*sigma2
      sigma.tld <- sqrt(sigma2.tld)
      whalf.tld <- sqrt(cal.weights)
      tqr.tld <- qr(mm * (whalf.tld/sigma.tld))
      attr(cal.weights,"qr") <- tqr.tld
      attr(cal.weights,"gwhalf") <- whalf.tld * sigma.tld 
# ---------------------------------------- #
#    VARIANCE ESTIMATION WORKS THIS WAY    #
# ---------------------------------------- #

    cal.weights
}
