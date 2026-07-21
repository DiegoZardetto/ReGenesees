`ext.calibrated2` <- function(data, ids, strata = NULL,
                              fpc = NULL, self.rep.str = NULL, check.data = TRUE,
                              weights.cal, calmodel, partition = FALSE, sigma2 = NULL){
##############################################################################
# This is an APPROXIMATE version of function ext.calibrated.                 #
#                                                                            #
# NOTE: The difference is that this version does NOT ask the user to specify #
#       the *base weights* that were externally calibrated. This piece of    #
#       information is usually UNAVAILABLE to users working with publicly    #
#       disseminated survey microdata.                                       #
#                                                                            #
# NOTE: The PRICE to pay is that sampling variance estimates are NO longer   #
#       EXACTLY equal to those that would be computed by the data producer;  #
#       rather, they will approximate those correct estimates with an ERROR  #
#       that becomes negligible in the large sample limit (n -> Inf).        #
##############################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
# directly <- !( length(sys.calls()) > 1 )  # Not used YET

# Prevent havoc caused by tibbles:
if (inherits(data, c("tbl_df", "tbl")))
    data <- as.data.frame(data)

# Currently cannot cope with NEGATIVE externally calibrated weights, thus
# handle a dedicated error
if (!inherits(weights.cal, "formula"))
     stop("Externally calibrated weights must be passed as a formula")
weights.cal.char <- all.vars(weights.cal)
if (length(weights.cal.char) < 1) 
     stop("Externally calibrated weights formula must reference a survey data variable")
if (length(weights.cal.char) > 1) 
     stop("Externally calibrated weights formula must reference only one variable")
na.Fail(data, weights.cal.char)
if (!is.numeric(data[, weights.cal.char])) 
     stop("Externally calibrated weights variable ", weights.cal.char, " is not numeric")
if (any(data[, weights.cal.char] <= 0)) 
     stop("Currently cannot handle negative externally calibrated weights, sorry!")

### BLOCK BELOW IS COMMENTED BECAUSE IT DOES NOT WORK YET ###
# NOTE: To skip e.svydesign's check on NEGATIVE weights, attach a 'pass'
#       attribute to data
# attr(data, "negw.pass") <- TRUE
### BLOCK ABOVE IS COMMENTED BECAUSE IT DOES NOT WORK YET ###

# Define a new design object by treating weights.cal as initial weights
design.new <- e.svydesign(data = data, ids = ids, strata = strata,
                          weights = weights.cal,
                          fpc = fpc, self.rep.str = self.rep.str, check.data = check.data)

# Desume known population totals from design.new:
pop <- aux.estimates(design.new, calmodel = calmodel, partition = partition)

# APPROXIMATION: set g-weights to one (this is the large sample limit n -> Inf)
g <- 1

# Standard checks on sigma2 here, as e.calibrate - that would perform the same
# checks - is invoked later on...
if (!is.null(sigma2)){
     if (!inherits(sigma2, "formula"))
         stop("Heteroskedasticity variable must be passed as a formula")
     sigma2.char <- all.vars(sigma2)
     if (length(sigma2.char) > 1) 
         stop("Heteroskedasticity formula must reference only one variable")
     na.Fail(design.new$variables, sigma2.char)
     if (!is.numeric(variance <- design.new$variables[, sigma2.char])) 
         stop("Heteroskedasticity variable must be numeric")
     if ( any(is.infinite(variance)) || any(variance <= 0) )
         stop("Heteroskedasticity variable must have finite and strictly positive values")
     ### In case external calibration involved heteroskedasticity, must prepare
     ### new variable g*sigma2
     g.sigma2 <- g*variance
    }
else {
     ### In case external calibration DID NOT involve heteroskedasticity, must
     ### prepare new variable g*1
     g.sigma2 <- g
    }

# Temporarily add to design.new convenience variable g.sigma2
design.new$variables[["g.sigma2"]] <- g.sigma2

# Perform the 'cosmetic' calibration step
# NOTE: Here bounds serve the only purpose of avoiding warnings in case the
#       calibration constraints were linearly dependent: finite bounds imply
#       using Newton-Raphson which is not sensitive to collinearity.
design.new <- e.calibrate(design = design.new, df.population = pop,
                          calfun = "linear", bounds = c(0.9, 1.1),
                          sigma2 = as.formula("~g.sigma2", env = .GlobalEnv))

# Drop variable g.sigma2
design.new$variables[["g.sigma2"]] <- NULL

# Add a token to testify external calibration
# NOTE: THIS TOKEN COULD (AND MUST) BE REMOVED BY ANY SUBSEQUENT CALL OF
#       e.calibrate
attr(design.new, "ext.cal") <- TRUE
## END postprocessing

# Catch the actual call and use it to update call slot of design.new 
design.new$call <- sys.call()

# Return external calibrated object
design.new
}





