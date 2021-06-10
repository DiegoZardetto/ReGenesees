##
## This inaugurates ReGenesees support to *Special Purpose Calibration* tasks
##
## See e.g.:
## - Lesage, Eric. (2011). The use of estimating equations to perform a
##   calibration on complex parameters. Survey Methodology. 37.
##
## NOTE: Here, auxiliary information to be used for calibration comes in the
##       form of population parameters that are complex, non-linear functions of
##       auxiliary variables (think, e.g., of multiple regression coefficients).
##       For this reason, calibration constraints have to be linearized first.
##       This, in turn, generates *synthetic linearized variables* z that are
##       the *actual calibration variables* the calibration algorithm will
##       eventually work on. Tipically, *control totals* for these synthetic
##       linearized variables z will all be *zero*.
##

##
## Functions to prepare a survey design object to be calibrated on *Multiple
## Regression Coefficients*
##


`prep.calBeta` <- function(design, model, Beta,
                           by = NULL, partition = FALSE, drop.z = TRUE) {
################################################################################
# Prepare a survey design object to be calibrated on *Multiple Regression      #
# Coefficients*.                                                               #
#                                                                              #
# NOTE: Multiple Regression Coefficients are COMPLEX NON STANDARD BENCHMARKS   #
#       for a calibration task. As these benchmarks are population parameters  #
#       that are non-linear functions of auxiliary variables, calibration      #
#       constraints have to be linearized first. To do that, unlike for        #
#       ordinary calibration on population totals, the *actual calibration     #
#       variables* z must be computed as functions of the variables defining   #
#       the linear model y ~ x ('model') whose regression coefficients ('Beta')#
#       are known.                                                             #
#       Note also that, once more at odds with ordinary calibration, the       #
#       computation of the *actual calibration variables* z (one for each Beta #
#       component) *requires* the specification of the Beta values since the   #
#       beginning.                                                             #
#                                                                              #
# NOTE: One can calibrate                                                      #
#       (i)  on the regression coefficients of a *single linear model*,        #
#            known possibly for *different subpopulations*                     #
#       (ii) on the regression coefficients of *different linear models*,      #
#            each known at the *overall population level*                      #
#                                                                              #
# NOTE: The computed z variables will be added to design, along with useful    #
#       metadata.                                                              #
################################################################################

# NOTE: No Nas are allowed in model variables!
# NOTE: No negative weights are allowed!
# NOTE: As collinearity could in principle manifest at the sample level, despite
#       being absent in the population, one must check for it!

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Only most essential check on input type (more in inner functions...)
if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")

if (!inherits(model, "formula") && !is.list(model))
    stop("Linear regression model(s) must be specified either as one single formula or as a list of formulas")

if (!is.numeric(Beta) && !is.list(Beta))
    stop("Linear regression coefficients must be specified either as one single numeric vector or as a list of numeric vectors")

if (!is.null(by)) {
    if (!inherits(by, "formula")) 
        stop("If specified, 'by' must be supplied as a formula")
    }

if (is.list(Beta)) {
     n.Beta <- length(Beta)
     if (n.Beta == 1L)
         stop("To specify a single Beta vector, please pass it as a numeric vector")

     is.vec <- sapply(Beta, function(b) is.numeric(b))
     if (any(!is.vec))
         stop("The list of linear regression Betas must contain only numeric vectors")

     if (!is.list(model) && is.null(by)) {
         stop("If 'Beta' is a list and 'by' is not passed, 'model' must be a list too. Its length must be equal to the length of 'Beta' (", n.Beta,")")
        }
    }

if (is.list(model)) {
     n.mod <- length(model)
     if (n.mod == 1L)
         stop("To specify a single linear model, please pass it as a formula")

     is.model <- sapply(model, function(mod) inherits(mod, "formula"))
     if (any(!is.model))
         stop("The list of linear regression models must contain only formula objects")

     if (!is.list(Beta)) {
         stop("If 'model' is a list, 'Beta' must be a list too. Its length must be equal to the length of 'model' (", n.mod,")")
        }

     if (!is.null(by))
         stop("Can specify just a single model if 'by' is passed")
    }

if (is.list(model) && is.list(Beta)) {
     if (n.mod != n.Beta) stop("Lists 'model' and 'Beta' must be the same length")
    }

if (!is.null(by)) {
     # Handle by argument:
     by.vars <- all.vars(by)
     na.Fail(design$variables, by.vars)
     typetest <- sapply(by.vars, function(v) is.factor(design$variables[, v]))
     if (!all(typetest)) 
         stop("All variables in 'by' must be factors")

     # Check 'model' consistency (model can only be a formula here, otherwise an
     # error would have been raised before)
     model.vars <- all.vars(model)
     if (any(by.vars %in% model.vars)) 
         stop("Variables referenced by argument 'by' must not appear in the input 'model' formula")

     #  'interact': factor whose levels identify 'by' sub-populations
     interact <- interaction(design$variables[, rev(by.vars), drop = FALSE], drop = TRUE)
     #  'groups': list storing the row indices of units belonging to different 'by' sub-populations 
     groups <- split(1:nrow(design$variables), interact)
     n.groups <- length(groups)
     group.names <- levels(interact)

     # Check Beta consistency
     if (!is.list(Beta)) {
         stop("If 'by' is specified, 'Beta' must be a list. Its length must be equal to the number of 'by' subpopulations (", n.groups,")")
        } else {
             if (n.Beta != n.groups) stop("'Beta' length (", n.Beta, ") does not match the number of 'by' subpopulations (", n.groups,")")
            }
    }

# If partitioned calibration requested, then design must be prepared with a 'by' argument
if (!is.logical(partition)) stop("Argument 'partition' must be logical")

if (isTRUE(partition)) {
     if (is.null(by))
         stop("Only if argument 'by' is specified you can prepare 'design' for a *partitioned* calibration task")
    }

# Ensure partition is TRUE only if by is pecified
if (is.null(by)) partition <- FALSE

# For coding convenience, reduce to list inputs only
if (input.list <- is.list(model)) {
     modelS <- model
     BetaS <- Beta
    } else {
     if (!is.null(by)) {
         modelS <- rep(list(model), n.Beta)
         names(modelS) <- group.names
         n.mod <- length(modelS)
         BetaS <- Beta
         names(BetaS) <- group.names
        } else {
          modelS <- list(model)
          n.mod <- 1
          BetaS <- list(Beta)
          n.Beta <- 1
        }
    }

# Should z variables eventually be dropped after calibration
if (!is.logical(drop.z)) stop("Argument 'drop.z' must be logical")

# Initialize "spec.purp.calib" attribute of the output prepared design object
spec.purp.calib <- vector(mode = "list", length = 8)
names(spec.purp.calib) <- c("type", "benchmark", "model", "z.lin", "by", "calmodel", "partition", "drop.z")
spec.purp.calib$type <- "calBeta"
spec.purp.calib$by <- by
spec.purp.calib$partition <- partition
spec.purp.calib$drop.z <- drop.z
if (input.list) {
     # These must be built during the loop
     spec.purp.calib$benchmark <- list()
     spec.purp.calib$model <- list()
     spec.purp.calib$z.lin <- list()
    }
if (!is.null(by)) {
     # This must be built during the loop
     spec.purp.calib$benchmark <- list()
    }


# Initialize maximal z names vector, for partition = TRUE and in case of 
# differential aliasing
z.lin.max <- character(0)

## Here LOOP on model and Beta list components
for (i in 1:n.mod) {
     model <- modelS[[i]]
     Beta <- BetaS[[i]]
     group <- if (!is.null(by)) groups[[i]] else 1:NROW(design$variables)

     ## Prepare Calibration Variables z ## START -------------------------------

     ## Strings for output messages thrown after checks (non-positive weights,
     ## aliasing, consistency of mm and Beta...)
     if (input.list) {
         input.i.msg <- paste(" for input list element: ", i, sep = "")
        } else {
         if (!is.null(by)) {
             input.i.msg <- paste(" for input 'by' subpopulation: ", group.names[i], sep = "")
            } else {
             input.i.msg <- NULL
            }
        }

     ## NOTE: Same sanity checks as in function linB() + NO NAs allowed here!
     # Extract all symbols names from model
     varnames <- all.vars(model)
     # Extract names of design variables
     desnames <- names(design$variables)
     # Are there unknown symbols in model?
     unknown.vars <- varnames[!(varnames %in% desnames)]
     if (length(unknown.vars) > 0)
         stop("Unknown variables in regression model formula: ", paste(unknown.vars, collapse = ", "), input.i.msg)
     # Check model variable types
     typetest <- sapply(varnames, function(y) is.numeric(design$variables[, y]) |
                                              is.factor(design$variables[, y])   )
     if (!all(typetest))
         stop("Detected non-numeric or non-factor variables in regression model", input.i.msg)
     # No NAs allowed in any model variables!
     na.Fail(design$variables, varnames)
     # No NAs allowed in Beta!
     if (anyNA(Beta))
         stop("Missing values detected in Beta vector", input.i.msg)

     ## Get basic structures: 0. response variable name
     ##                       1. model frame
     ##                       2. model matrix
     ##                       3. model response
     ##                       4. weights
     # 0. Check that the model has a response variable and get its name
     respName <- get.respName(model)
     # 1. Model frame
     mf <- model.frame(model, data = design$variables, na.action = na.fail)[group, , drop = FALSE]
     # 2. Model matrix
     #    NOTE: No NAs are allowed in model variables!
     mm <- model.matrix(model, data = mf, na.action = na.fail)
     # 3. Model response
       # get the values
     resp <- model.response(mf)
       # check type
     if (!is.numeric(resp))
         stop("Detected non-numeric regression model response", input.i.msg)

     # 4. Sampling unit weights
     ww <- weights(design)[group]

     # Check for zero or negative weights
     # Negative weights produce an error, zero weights can be handled
     ww.neg <- (ww < 0)
     ww.0   <- (ww <= 0) & !ww.neg 
     N.ww.neg <- sum(ww.neg)
     if (N.ww.neg > 0) {
          stop("Detected ", N.ww.neg,
               " observations with weight < 0!\n  ",
               input.i.msg)
        }

     ## No NAs from now on: everything is ok.
     # Take weights sqrt (safe, as no negative weights are allowed here)
     whalf <- sqrt(ww)

     # Compute weighted model matrix
     mm.whalf <- mm * whalf

     ## Mandatory collinearity check: redundant variables will be aliased
     ## > Collinearity check STARTS <
         # Compute the QR decomposition of the weighted model matrix 
         tqr <- qr(mm.whalf)
         # Compute sample regression coefficients (weighted)
         Beta.in <- qr.coef(tqr, resp * whalf)

         # Handle variable aliasing
         alias <- any(aliased <- is.na(Beta.in))
         if (alias) {
             var.aliased <- names(Beta.in)[aliased]
             warning("Variables ", paste(var.aliased, collapse = ", "), " have been aliased due to collinearity!\n ",
                     input.i.msg)
             # Reduce the weighted model matrix and re-compute its QR decomposition
             mm <- mm[, !aliased, drop=FALSE]
             mm.whalf <- mm * whalf
             # tqr <- qr(mm.whalf)
             # Re-compute regression coefficients
             # Beta.in <- qr.coef(tqr, resp * whalf)
            }

         # Compute residuals
         # ee <- qr.resid(tqr, resp * whalf) / whalf
         # Take care of NaN arising from 0 whalf (i.e. zero weights, if any)
         # ee[ww.0] <- 0

     ## > Collinearity check ENDS <

     # Compute (generalized, if needed) inverse of t(X)%*%W%*%X
     T <- crossprod(mm.whalf)
     Tm1 <- try(solve(T), silent = TRUE)

     # Residual collinearity after aliasing (this should happen only
     # for numerical reasons...)
     if (collin.res <- inherits(Tm1, "try-error")) {
          warning("Model matrix is singular: switching to Moore-Penrose generalized inverse.")
          ## No longer needed: ReGenesees now IMPORTS MASS
          # require(MASS)
          Tm1 <- ginv(T)
        }

     # Consistency checks on mm and Beta

     # Check that model matrix and Beta coefficients have the same length
     if (!identical(nc <- NCOL(mm), nB <- length(Beta))) {
         stop.len <- paste("Model matrix columns (", nc,
                           ") do not match the number of specified regression coefficients (", nB,
                           ")", input.i.msg,
                           "!\n  Please double-check the consistency of 'model' and 'Beta' input values.", sep = "")
         stop(stop.len)
        }

     # Check that model matrix columns and Beta coefficients have the same names (if any)
     if (!is.null(B.names <- names(Beta)) && !is.null(mm.names <- colnames(mm)) && !identical(B.names, mm.names)) {
         if (setequal(mm.names, B.names)) {
             # Names are possibly permuted, try to reorder *Beta*
             Beta <- Beta[mm.names]
             warn.perm <- paste("Names of regression coefficients and model matrix columns were not in the same order",
                                input.i.msg,
                                "\n  Beta components have been re-ordered!", sep = "") 
             warning(warn.perm)
            } else {
             warn.names <- paste("Names of regression coefficients do not agree with model matrix columns",
                                 input.i.msg,
                                 "!\n  Please double-check the consistency of 'model' and 'Beta' input values.", sep = "")
             warning(warn.names)
            }
        }

     # Compute EXTERNAL residuals (i.e. y residuals w.r.t. an outside Beta)
     ee.out <- resp - mm %*% cbind(Beta)

     # Compute Beta Woodruff transform (see my notes on paper)
     z.Beta.out <- tcrossprod(mm * as.numeric(ee.out), Tm1)

     # Give z columns a name that reminds of the name of the response variable
     # NOTE: MUST GIVE A SUBSCRIPT FOR MODEL i, to avoid possible duplicates!
     z.lin <- make.names(paste(respName, colnames(z.Beta.out), sep = "_"))
     if (input.list) {
         mod.i <- paste("mod", i, sep = "")
         z.lin <- paste(z.lin, mod.i, sep = "_")
        }
     if (!is.null(by)) {
         mod.i <- group.names[i]
         if (!isTRUE(partition)) {
             z.lin <- paste(z.lin, mod.i, sep = "_")
            } else {
             # If partition = TRUE must update maximal z names, in case
             # differential aliasing in 'by' groups happens!
             if (length(z.lin) > length(z.lin.max)) z.lin.max <- z.lin
            }
        }
     colnames(z.Beta.out) <- z.lin

     # Add linearized variables to design
     if (is.null(by)) {
         # Then new z coming from *new models* are added as ADDITIONAL COLUMNS
         # as z.Beta.out has ALL ROWS (row-complete zetas)
         design$variables <- cbind(design$variables, z.Beta.out)
        } else {
                 if (!isTRUE(partition)) {
                     # Then, to solve calibration *globally*, new z coming from
                     # the *same model* fitted to *i-th group* subpopulation are
                     # added as ADDITIONAL COLUMNS...
                     design$variables[group, z.lin] <- z.Beta.out
                     # ...but z.Beta.out has *not* ALL ROWS, and values outside
                     # the *i-th group* subpopulation must be set to 0:
                     design$variables[-group, z.lin] <- 0
                    } else {
                     # Then new z values coming from new by subsets are added as
                     # ROW VALUES as z.Beta.out has only SOME ROWS
                     design$variables[group, z.lin] <- z.Beta.out
                     # NOTE: Beware that aliasing may delete different z variables
                     #       in different 'by' groups (i.e. rows), which will
                     #       generate NAs. These must be put to 0 (in accordance
                     #       to what is done for partition = FALSE.
                     #       This must be done at the end of the loop, as the
                     #       maximal number of z vars will be known only then!
                     #       See z.lin.max!
                    }
        }

     # Append info to "spec.purp.calib" attribute
     if (input.list || (!is.null(by) && !isTRUE(partition) )) {
         spec.purp.calib$benchmark <- c(spec.purp.calib$benchmark, list(Beta))
         names(spec.purp.calib$benchmark)[i] <- mod.i
         spec.purp.calib$model <- c(spec.purp.calib$model, list(model))
         names(spec.purp.calib$model)[i] <- mod.i
         spec.purp.calib$z.lin <- c(spec.purp.calib$z.lin, list(z.lin))
         names(spec.purp.calib$z.lin)[i] <- mod.i
        } else {
             if (!is.null(by)) {
                 spec.purp.calib$benchmark <- c(spec.purp.calib$benchmark, list(Beta))
                 names(spec.purp.calib$benchmark)[i] <- mod.i
                 spec.purp.calib$model <- model
                 # Recall z.lin.max!
                 # spec.purp.calib$z.lin <- z.lin
                 spec.purp.calib$z.lin <- z.lin.max
                } else {
                         spec.purp.calib$benchmark <- Beta
                         spec.purp.calib$model <- model
                         spec.purp.calib$z.lin <- z.lin
                }
        }
    }

# If !is.null(by) and partition = TRUE must set to 0 possible NAs in z variables
# that could have been generated due to differential aliasing in 'by' groups!
if (!is.null(by) & isTRUE(partition)) {
     ZETAS <- design$variables[, z.lin.max]
     ZETAS[is.na(ZETAS)] <- 0
     design$variables[, z.lin.max] <- ZETAS
    }

# Add to "spec.purp.calib" attribute the overall calmodel
if (input.list || (!is.null(by) && !isTRUE(partition)) ) {
     calmodel.chr <- paste(" ~ ",
                           paste(sapply(spec.purp.calib$z.lin, function(z.mod) paste(z.mod, collapse = " + ")),
                                 collapse = " + "
                                ),
                           sep = "")
     # Remove the intercept
     spec.purp.calib$calmodel <- as.formula(paste(calmodel.chr, " -1", sep = ""), env = .GlobalEnv)
    } else {
     calmodel.chr <- paste("~", paste(spec.purp.calib$z.lin, collapse = " + "), sep = "")
     # Remove the intercept
     spec.purp.calib$calmodel <- as.formula(paste(calmodel.chr, " -1", sep = ""), env = .GlobalEnv)
    }

# Add the "spec.purp.calib" attribute to the return design object
attr(design, "spec.purp.calib") <- spec.purp.calib

attr(design, "design") <- design.expr
return(design)
}


`is.spc` <- function(obj){
################################################################################
# Is object 'obj' prepared for a *special purpose calibration* task?           #
#                                                                              #
# NOTE: Argument 'obj' can be either a *survey design* or a *control totals    #
#       data frame*.                                                           #
#                                                                              #
# NOTE: If 'obj' is a control totals data frame, now that the new 'spc.pop'    #
#       class has been defined, this function could be refactored (or avoided) #
#       using inherits(obj, "spc.pop").                                        #
################################################################################
  !is.null(attr(obj, "spec.purp.calib"))
}


`pop.calBeta` <- function(design) {
################################################################################
# Given a survey design object already prepared to be calibrated on *Multiple  #
# Regression Coefficients*, generate the corresponding *control totals data    #
# frame*.                                                                      #
# NOTE: The return object will have a new class, 'spc.pop', which specializes  #
#       the 'pop.totals' class. This class is meant to handle also other cases #
#       of *Special Purpose Calibration* tasks.                                #
################################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Only most essential check on input type (more in inner functions...)
if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")

# Check that design has actually been prepared
if (!is.spc(design))
     stop("Input design object is not ready for this task yet!\n  Use function prep.calBeta to prepare it first.")

spec.purp.calib <- attr(design, "spec.purp.calib")

# Check that design has been prepared for the right task
if (spec.purp.calib$type != "calBeta")
     stop("Input design object has been prepared for a different calibration task: ", spec.purp.calib$type)

# Build the template and fill it with zeros
if (isFALSE(spec.purp.calib$partition)) {
     # Build
     pop <- pop.template(design, calmodel = spec.purp.calib$calmodel)
     # Fill
     pop[, ] <- 0
}

if (isTRUE(spec.purp.calib$partition)) {
     # Build
     pop <- pop.template(design, calmodel = spec.purp.calib$calmodel, spec.purp.calib$by)
     # Fill
     partition.vars <- all.vars(spec.purp.calib$by)
     n.part.vars <- length(partition.vars)
     pop[, -(1:n.part.vars)] <- 0
    }

# Copy the "spec.purp.calib" derived from design and attach it to the return
# pop.totals object
attr(pop, "spec.purp.calib") <- spec.purp.calib

attr(pop, "data") <- design.expr

# New class "spc.pop"
class(pop) <- c("spc.pop", class(pop))

return(pop)
}


`pop.fuse` <- function(spc.pop, pop, design){
################################################################################
# This function fuses population totals for a "special purpose calibration"    #
# task ('spc.pop) and for an "ordinary calibration" task ('pop'), so that      #
# survey data ('design') can be calibrated simultaneously on both.             #
#                                                                              #
# NOTE: All *special purpose calibration* tasks entail using as auxiliary      #
#       variables synthetic linearized variables whose calibration benchmarks  #
#       are zero. Calibration tasks of this kind always admit a false solution #
#       where calibration weights are *all zero*. This can be prevented by     #
#       fusing these zero calibration benchmarks with known totals for an      #
#       ordinary calibration task, as the latter cannot be matched with        #
#       calibration weights that are identically zero.                         #
################################################################################

# Check that spc.pop is actually a control totals dataframe for a
# "spec.purp.calib" task
if (!inherits(spc.pop, "spc.pop"))
     stop("Input object 'spc.pop' must be a control totals dataframe for a *special purpose calibration* task")

# Check that spc.pop has the right class and is filled
if (!inherits(spc.pop, "pop.totals"))
     stop("Object 'spc.pop' must be of class pop.totals") # Redundant, unless user is hacking
if (any(is.na(spc.pop)))
     stop("Object 'spc.pop' cannot contain any NAs!")

# Check that pop has the right class and is filled
if (data.class(pop) != "pop.totals")
     stop("Object 'pop' must be of class pop.totals, i.e. a known totals dataframe for an *ordinary calibration* task")
if (any(is.na(pop)))
     stop("Object 'pop' cannot contain any NAs!")

# Check that spc.pop and pop are either (i) both global or (ii) both partitioned
# on the same partition variables
if (!isTRUE(all.equal(attr(spc.pop, "partition"), attr(pop, "partition"))))
     stop("To fuse objects 'spc.pop' and 'pop' their 'partition' attributes must agree")

# Check that 'spc.pop' structure agrees with 'design' variables
spc.pop.test <- pop.calBeta(design) 
if (!isTRUE(all.equal(spc.pop, spc.pop.test, check.attributes = FALSE)))
     stop("Input object 'spc.pop' does not agree with 'design'")

# Check that 'pop' structure agrees with 'design' variables
# NOTE: Here try() is just to avoid printing on screen that the test has been
#       passed
pop.test <- try(population.check(pop, design, attr(pop, "calmodel"), attr(pop, "partition"))) 
if (inherits(pop.test, "try-error"))
     stop("Input object 'pop' does not agree with 'design'")

# Warn if names of data used to build spc.pop and pop do not agree (if any)?
# if ( !identical(attr(spc.pop, "data"), attr(pop, "data")) &&
     # !identical(spc.pop, eval(attr(pop, "data")))         &&
     # !is.null(attr(pop, "data"))                          &&
     # !is.null(attr(spc.pop, "data")) ) 
         # warning("Data frames used to build objects 'spc.pop' and 'pop' have different names",
                 # immediate. = TRUE)

# Get pop's partition attribute
partition <- attr(pop, "partition")

# Is it a global calibration task?
is.global <- identical(partition, FALSE)

## Build the *fused* data frame ( order: pop U spc.pop )
if (!is.global) {
     # If calibration task is partitioned, get partition variables and do not
     # duplicate them
     partition.vars <- all.vars(partition)
     n.part.vars <- length(partition.vars)
     # Drop spc.pop partition columns to restrict to auxiliary info columns
     spc.aux.totals <- spc.pop[, -(1:n.part.vars), drop = FALSE]
     # Fuse pop and spc.pop *data.frames*
     fused.df <- cbind(pop, spc.aux.totals)
    } else {
     # Fuse pop and spc.pop *data.frames*
     fused.df <- cbind(pop, spc.pop)
    }

## Build the *fused* calmodel
# Get pop's calmodel
pop.calmodel <- attr(pop, "calmodel")
# Does 'pop' include an intercept term?
pop.has.intercept <- any(attr(pop, "which.term.col") == 0)
# Get spc.pop's calmodel
spc.pop.calmodel <- attr(spc.pop, "spec.purp.calib")$calmodel
# Get pop's calmodel string
pop.calmodel.char <- form.to.char(pop.calmodel)
if (!pop.has.intercept){
     # Strip '- 1' from pop.calmodel
     pop.calmodel.char <- gsub("- 1", "", pop.calmodel.char)
    }
# Get pop's calmodel string
spc.pop.calmodel.char <- form.to.char(spc.pop.calmodel)
# Strip '~' from spc.pop.calmodel
spc.pop.calmodel.char <- gsub("~", "", spc.pop.calmodel.char)
# Fuse pop and spc.pop *calmodels* ( order: pop U spc.pop )
calmodel <- as.formula(paste(pop.calmodel.char, spc.pop.calmodel.char, sep = " + "), env = .GlobalEnv)

# Build a template for the fused calmodel
# NOTE: This is indeed *necessary*, as terms order in fused calmodel will *not* 
#       automatically match the columns order of fused.df
fused.pop <- pop.template(design, calmodel, partition)

# Copy fused.df values inside fused.pop *taking into account the columns order*
fused.pop[, names(fused.df)] <- fused.df[, names(fused.df)]

# Attach original calmodel attributes of 'pop' to fused.pop
attr(fused.pop, "pop.calmodel") <- pop.calmodel 

# Attach all the "spec.purp.calib" attributes of 'spc.pop' to fused.pop
attr(fused.pop, "spec.purp.calib") <- attr(spc.pop, "spec.purp.calib")

# New class "spc.pop"
class(fused.pop) <- c("spc.pop", class(fused.pop))

# Return the fused population totals
fused.pop
}


`is.fused.spc` <- function(pop.totals){
#############################################################################
# Is pop.totals a fused known totals data frame with both *special purpose* #
# and *ordinary* calibration benchmarks?                                    #
#############################################################################
  !is.null(attr(pop.totals, "pop.calmodel"))
}


`pop.desc.spc.pop` <- function(pop.totals, verbose = FALSE, ...){
#############################################################################
# This function provides a natural language description of a known totals   #
# data frame to be used for a *special purpose calibration* task.           #
#                                                                           #
# NOTE: Currently, *only* the case of calibration on *Multiple Regression   #
#       Coefficients* is handled.                                           #
#       Will have to *extend* it, once new cases of *special purpose        #
#       calibration* tasks are supported by ReGenesees.                     #
#                                                                           #
# NOTE: This is an S3 method for class 'spc.pop'. An S3 method for class    #
#       'pop.totals' do exist and is *directly* invoked on the input object #
#       'pop.totals' when verbose = TRUE.                                   #
#                                                                           #
# NOTE: Object 'pop.totals' can either be *simple* or *fused*.              #
#############################################################################

if (!is.spc(pop.totals))
    stop("Input object must be a known totals dataframe for a *special purpose calibration* task")

spec.purp.calib <- attr(pop.totals, "spec.purp.calib")
is.global <- identical(attr(pop.totals, "partition"), FALSE)
cal.task <- if (is.global) "Global" else "Partitioned"

cat("# Data frame of control totals for a *special purpose calibration* task", "\n", sep = "")

## -- Multiple Regression Coefficients - START -- ##
is.calBeta <- (spec.purp.calib$type == "calBeta")
spc.type <- if (is.calBeta) "Regression Coefficients" else "*UNKNOWN*"

cat("- Benchmark parameters:       ", spc.type, "\n", sep = "")

if (is.calBeta) {
     models.list <- attr(pop.totals, "spec.purp.calib")$model
     has.by <- !is.null(spec.purp.calib$by)
     if (has.by) by.vars <- paste(all.vars(spec.purp.calib$by), collapse = ":")
     if (!has.by) {

cat("- Benchmarks known at level:  ", "Overall Population", "\n", sep = "")
cat("- Calibration task:           ", cal.task, "\n", sep = "")
cat("\n")
cat("## Regression models", "\n", sep = "")
cat("\n")
print(models.list, showEnv = FALSE)
if (!is.list(models.list)) cat("\n")
cat("## Benchmark vectors (Beta)", "\n", sep = "")
cat("\n")
print(attr(pop.totals, "spec.purp.calib")$benchmark)
if (!is.list(models.list)) cat("\n")

                } else {

cat("- Benchmarks known at level:  ", "Domains [", by.vars, "]", "\n", sep = "")
cat("- Calibration task:           ", cal.task, "\n", sep = "")
cat("\n")
cat("## Regression models", "\n", sep = "")
cat("\n")
if (is.list(models.list)) print(unique(models.list)[[1]], showEnv = FALSE) else print(models.list, showEnv = FALSE)
cat("\n")
cat("## Domain benchmark vectors (Beta)", "\n", sep = "")
print(attr(pop.totals, "spec.purp.calib")$benchmark)
cat("\n")

                }
    }
## -- Multiple Regression Coefficients - END -- ##


if (is.fused.spc(pop.totals) & !isTRUE(verbose)) {
     # Notify if pop.totals is *fused*
cat("## NOTE: This is a *fused* known totals dataframe (see ?pop.fuse)", "\n", sep = "")
cat("##       Specify verbose = TRUE to learn about its *ordinary calibration* totals!", "\n", sep = "")
cat("\n")
    }

if (isTRUE(verbose)) {
     # Report the full structure of the auxiliary variables, both synthetic and
     # ordinary (for *fused* objects)
cat("# < Auxiliary Variables Structure >\n")
cat("\n")
     # Call pop.desc method for pop.totals class, temporarily stripping the
     # specialized class "spc.pop"
     temp <- pop.totals
     class(temp) <- class(temp)[class(temp) != "spc.pop"]
     temp <- pop.desc(temp)
    }

return(invisible(pop.totals))
}
