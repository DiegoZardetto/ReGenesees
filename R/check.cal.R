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
    invisible(s)
}
