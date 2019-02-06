`ez.regcalib` <-
function (design, formula, population, aggregate.stage = NULL, sigma2 = NULL,
          epsilon, force)
##########################################################################
#  Versione modificata della funzione regcalibrate del package survey.   #
#  NOTA: La funzione ritorna DIRETTAMENTE i pesi calibrati e non un      #
#        oggetto design (questo fa risparmiare memoria).                 #
#  NOTA: Per motivi di efficienza sono state eliminate alcune            #
#        funzionalita' originali NON NECESSARIE per e.calibrate          #
#  NOTE: epsilon and force (with the same meaning as in Newton-Raphson)  #
#        introduced to deal with possible approximate solutions obtained #
#        when switching to ginv due to collinearity.                     #
##########################################################################
{
    mm <- model.matrix(formula, model.frame(formula, data = design$variables))
    ww <- design$dir.weights
    # Heteroskedsticity
    if (is.null(sigma2)) {
        sigma2 <- rep(1, NROW(mm))
        }
    else {
        sigma2 <- design$variables[, all.vars(sigma2)]
        }

    # in.whalf and in.sigma have the purpose to deal with aggregate calibration
    whalf <- in.whalf <- sqrt(ww)
    sigma <- in.sigma <- sqrt(sigma2)
    tqr <- qr(mm * (whalf/sigma))
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
#        whalf <- sqrt(ww)
#        sigma2 <- as.numeric(tapply(sigma2, aggindex, mean))
#        sigma <- sqrt(sigma2)
#    }
      # NEW: rowsum is definitely faster
    if (!is.null(aggregate.stage)) {
        aggindex <- design$cluster[[aggregate.stage]]
        n.aggindex <- rowsum(rep(1, NROW(mm)), aggindex)
        mm <- rowsum(mm, aggindex)
        ww <- as.numeric(rowsum(ww, aggindex) / n.aggindex)
        whalf <- sqrt(ww)
        sigma2 <- as.numeric(rowsum(sigma2, aggindex) / n.aggindex)
        sigma <- sqrt(sigma2)
    }
    g <- rep(1, NROW(mm))
    Tmat <- crossprod(mm * (whalf/sigma))
    # Let's try a workaround for collinearity problems
    tT <- try(solve(Tmat, population - sample.total), silent = TRUE)
    if (collin <- inherits(tT, "try-error")) {
        warning("Calibration system is singular: switching to Moore-Penrose generalized inverse.")
        tT <- ginv(Tmat) %*% (population - sample.total)
    }
    # done.

    g <- drop(1 + mm %*% tT/sigma2)

    # Calibration Diagnostics:
      # Build a return code and a diagnosys dataframe
        # C1: deviations BELOW epsilon (be collin or !collin) -> ret.code = 0
        # C2: deviations ABOVE epsilon (be collin or !collin):
          # C2.1: if !force -> FAIL (no ret.code, default -1 applies)
          # C2.2: if  force -> ret.code = 1
    # NOTE: Even if !collin (that is: solve goes ok and ginv is avoided)
    #       there could be cases of missed calibration constraints, due
    #       to numerically instable solutions...
    # NOTE: 'collin' flag is currently unused: it might be stored as a
    #       further attribute of ret.code...   
    cal.sample.total <- colSums(mm * ww * g)
    misfit <- cal.sample.total - population
    achieved <- (abs(misfit)/(1 + abs(population)))
    # This is C1:
    if (all(achieved < epsilon)) {
        ret.code <- 0
    }
    # This is C2: assess deviations and prepare a diagnostic structure
    # to be sent to ecal.status, named 'fail.diagnostics'
    else {
        worst.ind <- which.max(achieved)
        worst.achieved <- achieved[worst.ind]
        warning("Calibration failed: worst achieved epsilon= ", worst.achieved,
                " (variable ",names(worst.achieved),"), see ecal.status.")
        # This is C2 when collin: unstable results may be due to ginv
        # and a viable alternative may be to resort to Newton-Raphson using 
        # very loose bounds, e.g. bounds=c(-1E12, 1E12)
        if (collin) {
            warning("Using Moore-Penrose inverse for unbounded linear calibration can sometimes generate 'false' results: you may try, instead, to calibrate with very loose bounds, e.g. bounds=c(-1E12, 1E12).")
        }
        # This is C2.1
        if (!force) {
            stop("Calibration failed")
        }
        # This is C2.2
        ret.code <- 1
        # Build diagnosys dataframe
        ko.ind <- which(achieved >= epsilon)
        reldiff <- misfit/(1 + population)
        ko.reldiff <- reldiff[ko.ind]
        ko.names <- names(ko.reldiff)
        ko.pop <- population[ko.ind]
        ko.est <- cal.sample.total[ko.ind]
        ko.diff <- ko.est - ko.pop
        diagnosys <- data.frame(`Variable` = ko.names,
                                `Population.Total` = ko.pop,
                                `Achieved.Estimate` = ko.est,
                                `Difference` = ko.diff,
                                `Relative.Difference` = ko.reldiff,
                                row.names = NULL)
        # Set missed constraints indices as useful rownames, i.e. the
        # positions of non convergence cells in df.population for each
        # partition (if any):
        rownames(diagnosys) <- ko.ind
        diagnosys <- diagnosys[order(abs(ko.reldiff), decreasing = TRUE), , drop = FALSE]
    }


    # For aggregated calibration must re-expand g-weights
    if (!is.null(aggregate.stage)) {
        g <- g[aggindex]
    }

    # Handle zero calibrated weights (if any): must provide a tiny cutoff
    # to prevent problems in svyrecvar (fictitious 0/0=NaN which causes an error
    # in qr). Results DO NOT DEPEND on the chosen cutoff value.
    g.0 <- (abs(g) < 1E-12)
    if (any(g.0)) g[g.0] <- (1E-12) * ifelse(sign(g[g.0])==0, 1, sign(g[g.0]))


    cal.weights <- g*design$dir.weights
    attr(cal.weights,"qr") <- tqr
    attr(cal.weights,"gwhalf") <- g * in.whalf * in.sigma
    attr(cal.weights,"ret.code") <- ret.code
    fail.diagnostics <- if (ret.code==1) diagnosys
    attr(attr(cal.weights,"ret.code"), "fail.diagnostics") <- fail.diagnostics
    cal.weights
}
