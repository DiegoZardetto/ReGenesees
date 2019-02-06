`ez.grake` <-
function (sample.total, mm, ww, calfun, eta = rep(0, NCOL(mm)), bounds, population, 
    epsilon, maxit, sigma2)
###################################################################
#  Versione modificata della funzione grake del package survey.   #
#  NOTA: Per motivi di efficienza sono state eliminate alcune     #
#        funzionalita' originali NON NECESSARIE per e.calibrate   #
#        (require di MASS e attr(g,"eta")).                       #
###################################################################
{
    if (!inherits(calfun, "calfun")) 
        stop("'calfun' must be of class 'calfun'")
    Fm1 <- calfun$Fm1
    dF <- calfun$dF
    xeta <- drop(mm %*% eta)
    g <- 1 + Fm1(xeta, bounds, sigma2)
    iter <- 1
    repeat ({
        Tmat <- crossprod(mm * ww * dF(xeta, bounds, sigma2), mm)
        misfit <- (population - sample.total - colSums(mm * ww * 
            Fm1(xeta, bounds, sigma2)))
        deta <- ginv(Tmat, tol = 256 * .Machine$double.eps) %*% 
            misfit
        eta <- eta + deta
        xeta <- drop(mm %*% eta)
        g <- 1 + Fm1(xeta, bounds, sigma2)
        misfit <- (population - sample.total - colSums(mm * ww * 
            Fm1(xeta, bounds, sigma2)))
        if (all(abs(misfit)/(1 + abs(population)) < epsilon)) 
            break
        iter <- iter + 1
        if (iter > maxit) {
            achieved <- (abs(misfit)/(1 + abs(population)))
            worst.ind <- which.max(achieved)
            worst.achieved <- achieved[worst.ind]
            warning("Failed to converge: worst achieved epsilon= ", worst.achieved, " in ", 
                iter, " iterations (variable ",names(worst.achieved),"), see ecal.status.")
            # Build a diagnostic structure to be sent to ecal.status
            ko.ind <- which(achieved >= epsilon)
            reldiff <- -(misfit/(1 + population))
            ko.reldiff <- reldiff[ko.ind]
            ko.names <- names(ko.reldiff)
            ko.pop <- population[ko.ind]
            ko.est <- (population - misfit)[ko.ind]
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
            # Attach diagnostics to g as an attribute
            attr(g, "failed") <- achieved
            attr(g, "fail.diagnostics") <- diagnosys
            break
        }
    })
    g
}
