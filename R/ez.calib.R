`ez.calib` <-
function (design, formula, population, aggregate.stage = NULL, sigma2 = NULL,
          bounds = c(-Inf, Inf), calfun = c("linear", "raking", "logit"),
          maxit = 50, epsilon = 1e-07, force = FALSE)
########################################################################
#  Versione modificata della funzione calibrate del package survey.    #
#  NOTA: La funzione ritorna DIRETTAMENTE i pesi calibrati e non un    #
#        oggetto design (questo fa risparmiare memoria).               #
#  NOTA: Se aggregate.stage non e' NULL, per diminuire l'uso di        #
#        memoria e la complessita' computazionale, la funzione:        #
#        1) media la model matrix (mm) e somma i pesi diretti (ww)     #
#           nei cluster di stadio aggregate.stage, RIDUCENDO LA        #
#           DIMENSIONE DELL'INPUT di ez.grake;                         #
#        2) calibra i DATI RIDOTTI;                                    #
#        3) espande gli g-weights calcolati in 2) sulle unita' finali. #
#  NOTA: Per consentire la diagnosi dei processi di calibrazione       #
#        complessi e strutturati in molti sub-tasks (tipici di         #
#        e.calibrate) e' stato introdotto il NUOVO attributo           #
#        'ret.code' per l'oggetto (i pesi calibrati) ritornato dalla   #
#        funzione:                                                     #
#        - ret.code=0 se la convergenza E' ottenuta                    #
#        - ret.code=1 se la convergenza NON E' ottenuta MA force=TRUE  #
#  NOTA: ez.calib chiama ez.grake.                                     #
#  NOTA: sigma2 e' una formula che referenzia un vettore numerico a    #
#        valori finiti e strettamente positivi.                        #
########################################################################
{
    if (is.character(calfun)) 
        calfun <- match.arg(calfun)
    if (is.character(calfun) && calfun == "linear" && (bounds == 
        c(-Inf, Inf))) {
        cal.weights <- ez.regcalib(design, formula, population,
                                   aggregate.stage = aggregate.stage,
                                   sigma2 = sigma2, epsilon = epsilon,
                                   force = force)
        return(cal.weights)
    }
    if (is.character(calfun)) 
        calfun <- switch(calfun, linear = cal.linear, raking = cal.raking, 
            logit = cal.logit)
    else if (!inherits(calfun, "calfun"))
        stop("'calfun' must be a string or of class 'calfun'.")
    mm <- model.matrix(formula, model.frame(formula, data = design$variables))
    ww <- design$dir.weights
    # Heteroskedsticity
    if (is.null(sigma2)) {
        sigma2 <- rep(1, NROW(mm))
        }
    else {
        sigma2 <- design$variables[, all.vars(sigma2)]
        }

    whalf <- sqrt(ww)
    sigma <- sqrt(sigma2)
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
    g <- ez.grake(sample.total, mm, ww, calfun, bounds = bounds, population = population,
                  epsilon = epsilon, maxit = maxit, sigma2 = sigma2)
    if (!force && !is.null(attr(g, "failed"))) 
        stop("Calibration failed")
    ret.code <- if (is.null(attr(g, "failed"))) 0 else 1
    fail.diagnostics <- attr(g, "fail.diagnostics")

    # For aggregated calibration must re-expand g-weights
    if (!is.null(aggregate.stage)) {
        g <- g[aggindex]
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
    if (any(g.0)) g[g.0] <- (1E-12) * ifelse(sign(g[g.0])==0, 1, sign(g[g.0]))

    cal.weights <- g * design$dir.weights
    attr(cal.weights,"ret.code") <- ret.code
    attr(attr(cal.weights,"ret.code"), "fail.diagnostics") <- fail.diagnostics
    attr(cal.weights,"qr") <- tqr
    attr(cal.weights,"gwhalf") <- g * whalf * sigma 
    cal.weights
}
