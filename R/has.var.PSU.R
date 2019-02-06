has.var.PSU <- function(design) {
###################################################
# Check if the survey design has "variance PSUs", #
# i.e. has been built by specifying self.rep.str. #
###################################################
# test
has <- !identical(attr(design,"self.rep.str"), FALSE)
# check
if ( ("var.PSU" %in% names(design$variables)) && !has )
   stop("This should not happen! Please report this message to the author (zardetto@istat.it)")
has
}
