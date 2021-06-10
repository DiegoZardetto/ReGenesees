`check.cal` <-
function (cal.design)
#####################################################################
#  Dato un oggetto di classe cal.analytic, verifica la convergenza  #
#  dell'algoritmo di calibrazione.                                  #
#####################################################################
{
    if (!inherits(cal.design, "cal.analytic")) 
        stop("Object 'cal.design' must be of class cal.analytic")

    N.cal.constr <- attr(cal.design, "N.cal.constr")
    eps <- attr(cal.design, "epsilon")
    s <- attr(cal.design, "ecal.status")
    rc <- s$return.code
    failed <- any(rc!=0)
    if (failed) {
        diagn <- s$fail.diagnostics
        nfail <- sum(sapply(diagn, nrow))
        cat("Calibration Constraints missed (at tolerance level epsilon = ",
            eps, "): ", nfail," out of ", N.cal.constr, "\n", sep="")
        cat("- Summary of mismatches: \n\n")
        print(s[c("return.code", "fail.diagnostics")])
        cat("\n")
       }
    else {
        cat("All Calibration Constraints (", N.cal.constr,") fulfilled (at tolerance level epsilon = ",
            eps, ").\n",sep="")
        cat("\n")
       }
    s$eps <- eps

## SPC ZONE - START
# If cal.design has just undergone a *special purpose calibration task*,
# check for possible *false convergence* (i.e. all w.cal ~ 0)
    if ( isTRUE(attr(cal.design, "spc.justdone")) ) {
        # Here, the threshold for condition (all w.cal ~ 0) is set in such a way
        # that the calibrated estimate of the population count would drop below
        # one: N.cal < 1
         if ( sum(abs(weights(cal.design))) < 1 ) {
             if (!failed) { 
                 warning("All calibration weights collapsed to (nearly) zero! This indicates a FALSE CONVERGENCE of the *special purpose calibration* task.")
                } else {
                 warning("All calibration weights collapsed to (nearly) zero!")
                }
            }
        }
## SPC ZONE - END

    invisible(s)
}
