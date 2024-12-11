`ext.calibrated` <- function(data, ids, strata = NULL, weights,
                             fpc = NULL, self.rep.str = NULL, check.data = TRUE,
                             weights.cal, calmodel, partition = FALSE, sigma2 = NULL){
############################################################################
# Enable ReGenesses to digest calibration weights that have been computed  #
# externally (e.g. by other software), so as to provide correct variance   #
# estimates (i.e. appropriate to calibration estimators).                  #
############################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
# directly <- !( length(sys.calls()) > 1 )  # Not used YET

# Prevent havoc caused by tibbles:
if (inherits(data, c("tbl_df", "tbl")))
    data <- as.data.frame(data)

# Standard checks on weights here, as e.svydesign - that would perform the same
# checks - will be invoked later using weights.cal NOT weights...
if (!inherits(weights, "formula"))
     stop("Weights must be passed as a formula")
weights.char <- all.vars(weights)
if (length(weights.char) < 1) 
     stop("Weights formula must reference a survey data variable")
if (length(weights.char) > 1) 
     stop("Weights formula must reference only one variable")
na.Fail(data, weights.char)
if (!is.numeric(data[, weights.char])) 
     stop("Weights variable ", weights.char, " is not numeric")
if (any(data[, weights.char] <= 0)) 
     stop("Direct weights must be positive!")

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
                          weights = weights.cal,                      # here
                          fpc = fpc, self.rep.str = self.rep.str, check.data = check.data)

# Desume known population totals from design.new:
pop <- aux.estimates(design.new, calmodel = calmodel, partition = partition)

# Desume g-weights:
g <- weights(design.new)/data[, weights.char]

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
#       using Newton-Raphson which is not sensible to collinearity.
design.new <- e.calibrate(design = design.new, df.population = pop,
                          calfun = "linear", bounds = c(0.9, 1.1),
                          sigma2 = as.formula("~g.sigma2", env = .GlobalEnv))

# Drop variable g.sigma2
design.new$variables[["g.sigma2"]] <- NULL

## START postprocessing design.new metadata so as to match the standard
## behaviour one would have obtained by calibrating directly within ReGenesees
   # Get the name of the "cosmetic" calibrated weights, namely the string
   # obtained by pasting weights.cal name and ".cal" string
cosm.cal.char <- paste(all.vars(weights.cal), ".cal", sep = "")
   # Change it into the name one would have obtained by calibrating directly
   # within ReGenesees, namely the string obtained by pasting weights name and
   # ".cal" string
cosm.cal.char.new <- paste(weights.char, ".cal", sep = "")
   # 1) Change cosmetic calibration column name
      # If a column named in the same way already exist, change its name by
      # appending a .old subfix
      names(design.new$variables)[names(design.new$variables) == cosm.cal.char.new] <- paste(cosm.cal.char.new, ".old", sep = "")
names(design.new$variables)[names(design.new$variables) == cosm.cal.char] <- cosm.cal.char.new
   # 2) Change metadata and their names accordingly
   #    NOTE: Recall attribute 'allprob' stores the reciprocal of the 
   #          *initial weights* in ALL ReGenesees weights-changing functions!
design.new$allprob[[1]] <- 1 / data[, weights.char]
names(design.new$allprob) <- weights.char
attr(design.new, "weights") <- as.formula(paste("~", cosm.cal.char.new, sep = ""), env = .GlobalEnv)
   # 3) Change sigma2 metadata
attr(design.new, "sigma2") <- sigma2
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


is.ext.cal <- function(design){
  isTRUE(attr(design, "ext.cal"))
}

