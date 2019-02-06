find.lon.strata <- function(design){
#################################################
# Given a stratified design object, find strata #
# containing lonely PSUs (if any).              #
# NOTE: Returns the levels of the strata with   #
#       lonely PSUs.                            #
# NOTE: If no lonely strata are found, returns  #
#       NULL invisibly.                         #
#################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum. (affects ONLY clps.strata.status)
directly <- !( length(sys.calls()) > 1 )
this.call <- sys.call()

if (!inherits(design, "analytic")) 
    stop("Design object must be of class analytic")

# Check if design is actually stratified
if (!isTRUE(design$has.strata)) {
     stop("Design object is not stratified!")
     return(invisible(NULL))
    }

  analyze.strata <- function(design){
  ########################################################
  # Processes the fpc slot of a design (i.e. design$fpc) #
  # in such a way that a two components list storing PSU #
  # counts inside strata for both population and sample  #
  # (popsize and sampsize slots) is returned.            #
  ########################################################
  des.fpc <- design$fpc
  des.strat <- design$strata[, 1]
  strata.PSUcounts <- function(PSU, strata){
  tapply(PSU, strata, unique)
  }
  lapply(des.fpc, function(el) { if (is.null(el))
                                    # i.e. for popsize component,
                                    # when no fpc have been specified
                                    return(el)
                                 else strata.PSUcounts(el[, 1], des.strat) })
  }

analysis <- analyze.strata(design)

####################################################
# Given the output of analyze.strata, finds strata #
# with lonely PSUs.                                #
# If no lonely PSUs are found gives an error.      #
####################################################
  if (is.null(pop <- analysis$popsize)) {
     # i.e. when no fpc have been specified
     samp <- analysis$sampsize
     if  (any(lonely <- (samp==1))) {
         str.lev <- names(samp)
         return(str.lev[lonely])
        }
     else {
         cat("# No lonely PSUs found!\n\n")
         return(invisible(NULL))
        }
    }
  else {
         samp <- analysis$sampsize
         str.lev <- names(samp)
         # exlude certainty PSUs
         lonely <- ((samp==1) & (pop>1))
         if  (any(lonely)) {
             return(str.lev[lonely])
            }
         else{
              cat("# No lonely PSUs found!\n\n")
              return(invisible(NULL))
            }
    }
stop("This should not happen! Please report this message to the author (zardetto@istat.it)")
}
