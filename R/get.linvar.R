`get.linvar` <- function(design, expr, by = NULL, stack = FALSE, na.rm = FALSE) {
##########################################################################################
# This function generates the matrix of linearized variables that svystatL would use to  #
# estimate the sampling variance of a complex estimator (i.e. any smooth function of HT  #
# or CALIBRATION estimators) by domains.                                                 #
#                                                                                        #
# NOTE: Owing to the theory of Taylor linearization, each domain will correspond to a    #
#       different linearized variable, i.e. a different column of the output matrix.     #
#                                                                                        #
# NOTE: Only for functions of HT estimators the columns of the output matrix can be      #
#       stacked into one single column without information loss (as their value is zero  #
#       outside the active domain). The same result cannot be obtained for functions of  #
#       CALIBRATION estimators, because the g-weighted residuals are generally non-zero  #
#       outside the active estimation domain. Therefore, argument 'stack' does nothing   #
#       if 'design' is calibrated.                                                       #
#                                                                                        #
##########################################################################################
  # Check on design
  if (!inherits(design, "analytic")) {
       stop("Object 'design' must inherit from class analytic")
      }

  z <- get.linvarHT(design = design, expr = expr, by = by, na.rm = na.rm, stack = stack) 

  if (!inherits(design, "cal.analytic")) {
         ## - Functions of HT estimators
         return(z)
        } else {
         ## - Functions of CALIBRATION estimators
         # Avoid (unlikely) names collisions in design
         z.names <- paste("ThIs", colnames(z), sep = ".")
         colnames(z) <- z.names
         design$variables <- cbind(design$variables, z)
         for (z.var in z.names) {
                 z.var.formula <- as.formula(paste("~", z.var, sep = ""))
                 ## Missing values treatment (to be kind to get.residuals)
                 nas <- is.na(z[, z.var])
                 if (na.rm && sum(nas) > 0) {
                     design_new <- design[!nas, ] # Thus weights of obs with NAs in z are now 0
                     design_new$variables[[z.var]][nas] <- 0 # Thus NAs in z are now 0
                    } else {
                     design_new <- design
                    }
                 ge_z.var <- get.residuals(cal.design = design_new, y = z.var.formula, scale = "g")
                 z[, z.var] <- ge_z.var
                 # if (na.rm && sum(nas) > 0) {
                     # Q) Set the g-weighted residuals back to NA (as it would happen for HT)?
                     # A) No, as this is not the way svystatL (and svyrecvar) would behave.

                     # z[nas, z.var] <- NA

                     #    REMARK: Those functions would rather give to the NAs input values
                     #            g-weighted residuals obtained by (i) first, putting NAs
                     #            to 0, and (ii) then, giving 0 weight to the corresponding 
                     #            sample units.
                     #    NOTE: There are no easy alternatives to the above, as (A) internal
                     #          functions like qr.resid call Fortran routines that cannot
                     #          handle NAs, and (B) with calibration estimators, you cannot
                     #          simply subset the data by dropping the rows containing NAs.
                 #    }
                }
         colnames(z) <- gsub("ThIs.z", "gez", colnames(z))
        }

  return(z)
}


`get.linvarHT` <- function(design, expr, by = NULL, na.rm = FALSE, stack = FALSE) {
##########################################################################################
# This function generates the matrix of linearized variables that svystatL would use to  #
# estimate the sampling variance of a complex estimator (i.e. any smooth function of HT  #
# estimators) by domains.                                                                #
#                                                                                        #
# NOTE: Owing to the theory of Taylor linearization, each domain will correspond to a    #
#       different linearized variable, i.e. a different column of the output matrix.     #
#                                                                                        #
# NOTE: For complex functions of calibration estimators, the further linearization step  #
#       entailing the g-weighted residuals of z, i.e. g * e_z, is left to get.linvar.    #
#       Therefore, if design is calibrated, its calibration metadata will be stripped.   #
#       The full linearized variable of any complex function of calibration estimators   #
#       is then obtained by composing functions get.residuals and get.linvarHT, as       #
#       illustrated below:                                                               #
#                                                                                        #
# -------------------------------------------------------------------------------------- #
# sbscal <- des.addvars(sbscal, Z = get.linvarHT(sbscal, expression(va.imp2/emp.num)))   #
# sbsdes <- des.addvars(sbsdes, geZ = get.residuals(sbscal, ~Z, scale = "g"))            #
# svystatL(sbscal, expression(va.imp2/emp.num))                                          #
# svystatTM(sbsdes, ~geZ)                                                                #
# identical( as.numeric(SE(svystatL(sbscal, expression(va.imp2/emp.num)))),              #
#            as.numeric(SE(svystatTM(sbsdes, ~geZ))) )                                   #
# -------------------------------------------------------------------------------------- #
#                                                                                        #
##########################################################################################
  # Check on design
  if (!inherits(design, "analytic")) {
       stop("Object 'design' must inherit from class analytic")
      }

  # If design is calibrated, treat it as an ordinary one to make domain subsetting more
  # straightforward
  if (inherits(design, "cal.analytic")) {
       class(design) <- class(design)[class(design) != "cal.analytic"]
       design$postStrata <- NULL
       stack = FALSE
      }

  if (is.null(by)) {
     z <- getLin(design = design, expr = expr, na.rm = na.rm)
     out <- as.matrix(z)
     colnames(out) <- "z"
    } else {
     if (!inherits(by, "formula")) {
             stop("If specified, 'by' must be supplied as a formula")
        }

     # Process by domains as svyby would do
     byfactors<-model.frame(by, model.frame(design), na.action=na.pass)
     ## all combinations that actually occur in this design
     byfactor<-do.call("interaction", byfactors)
     dropped<- weights(design,"sampling")==0
     uniquelevels<-unique(byfactor[!dropped])
     uniques <- match(uniquelevels, byfactor)

     if (isFALSE(stack)) {
         out <- sapply(uniques, function(i) {
                 domain.index <- byfactor %in% byfactor[i]
                 domain.z <- rep(0, nrow(design))
                 z <- getLin(design[byfactor %in% byfactor[i], ], expr = expr, na.rm = na.rm)
                 domain.z[domain.index] <- z
                 as.matrix(domain.z)
                }
            )
         colnames(out) <- paste("z", uniquelevels, sep = "_")
         # Permute columns to reconstruct svyby's output order
         out <- out[, order(byfactor[uniques])]
        } else {
         out <- rep(0, nrow(design))
         for (i in uniques) {
                 domain.index <- byfactor %in% byfactor[i]
                 z <- getLin(design[byfactor %in% byfactor[i], ], expr = expr, na.rm = na.rm)
                 out[domain.index] <- z
                }
         out <- as.matrix(out)
         colnames(out) <- paste("z", paste(all.vars(by), collapse = "_"), sep = "_")
        }
    }

  return(out)
}


`getLin` <- function(design, expr, na.rm = FALSE, ...){
##########################################################################################
# This function "extracts" the linearized variable z that svystatL would use to estimate #
# the sampling variance of a complex estimator.                                          #
#                                                                                        #
# NOTE: The complex estimator is any smooth function of HT estimators regardless whether #
#       design is uncalibrated or calibrated.                                            #
##########################################################################################
# NOTE: Recall design can enter here only if (i) it was uncalibrated or if (ii) it had its
#       calibration metadata stripped by get.linvarHT
z <- linearize(expr = expr, design = design, na.rm = na.rm, ...)

## Missing values treatment: any NAs of input variables would be translated to NAs of the
## output linearized variable z
if (na.rm && !is.null(nas <- attr(z,"nas"))){
     z_na <- rep(NA, length(nas))
     z_na[nas == 0] <- z
     z <- z_na
    }

return(z)
}
