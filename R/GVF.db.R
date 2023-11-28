#######################################################################
# Archive of "registered" (i.e. built-in and user-defined) GVF models #
# supported by ReGenesees.                                            #
#                                                                     #
# The GVF.db data frame has columns:                                  #
# - Model.id:       an integer key ranging from 1 to the number of    #
#                   rows of GVF.db$db.                                #
# - GVF.model:      a character that can be transformed into a        #
#                   two-sided formula. Currently, the latter can      #
#                   involve only the following variables:             #
#                   <> 'Y', 'SE', 'CV', 'VAR', 'DEFF' <>              #
# - Estimator.kind: a character giving the kind of estimators for     #
#                   which the GVF model is deemed to be appropriate.  #
#                   Currently, valid kinds are the following:         #
#                   <> 'Total', 'Mean', 'Frequency',                  #
#                      'Absolute Frequency', 'Relative Frequency',    #
#                      'Ratio', 'Share', 'Share Ratio',               #
#                      'Regression Coefficient', 'Quantile',          #
#                      'Complex Estimator' <>                         #
# - Resp.to.CV:     a character that can be transformed into a valid  #
#                   R expression. The latter must map the response of #
#                   the GVF model (namely: variable 'resp') to 'CV',  #
#                   and may depend parametrically only on variable    #
#                   'Y'. !!! ARE YOU SURE? MUST REFLECT ON THIS !!!   #
#                                                                     #
# NOTE: Custom GVF models, beyond the ones available at startup, must #
#       be "registered" by using the dedicated accessor function      #
#       GVF.db$insert().                                              #
# NOTE: Models can be deleted from the archive via dedicated accessor #
#       function GVF.db$delete().                                     #
# NOTE: GVF archive can be get (e.g. in order to be copied outside    #
#       the Namespace) via GVF.db$get().                              #
# NOTE: GVF archive can be overwritten (e.g. to use a customized GVF  #
#       archive which was exported in a previous ReGenesees session)  #
#       via GVF.db$assign().                                          #
# NOTE: GVF archive can be reset to its default values by calling     #
#       GVF.db$reset().                                               #
# NOTE: All accessor functions unlock the package Namespace and use   #
#       lexical scoping to update the archive.                        #
#######################################################################

###############################################
# DB of builtin GVF models enabled at startup #
###############################################
GVF.db.ini <- data.frame(
                         Model.id = 1:5,
                         GVF.model = c("log(CV^2) ~ log(Y)",
                                       "CV^2 ~ I(1/Y)",
                                       "CV^2 ~ I(1/Y) + I(1/Y^2)",
                                       "SE ~ Y + I(Y^2)",
                                       "CV ~ I(1/Y) + Y"),
                         Estimator.kind = c("Frequency",
                                            "Frequency",
                                            "Frequency",
                                            "Total",
                                            "Total"),
                         Resp.to.CV = c("sqrt(exp(resp))",
                                        "sqrt(resp)",
                                        "sqrt(resp)",
                                        "resp/Y",
                                        "resp"),
                          stringsAsFactors = FALSE
                        )


GVF.db <- list(
     ######################
     # Set startup GVF db #
     ######################
     db = GVF.db.ini,
     ###################################################
     # Insert a *new*, *legal* GVF model inside the DB #
     ###################################################
     insert = function(GVF.model, Estimator.kind = NA, Resp.to.CV = NA, verbose = TRUE) {
                    #### Checks on GVF.model
                       # Handle formula and character GVF.model
                       this.env <- environment()
                       if (!inherits(GVF.model, "formula")) {
                           if (!is.character(GVF.model)) {
                               stop("GVF.model must be either a character or a formula")
                            }
                           GVF.model <- as.formula(GVF.model, env = this.env)
                        }
                       else {
                           environment(GVF.model) <- this.env
                        }
                       # Run all checks on GVF.model at once
                       GVF.moldel.check(GVF.model)
                       # Check if model to insert is actually new
                       mod.old <- lapply(GVF.db$db$GVF.model, FUN = as.formula, env = this.env)
                       mod.dup <- sapply(mod.old, function(mod) identical(mod, GVF.model))
                       if (any(mod.dup)) {
                           stop("GVF model already registered!")
                        }
                    #### Checks on Estimator.kind
                       if (!is.na(Estimator.kind)) {
                           # Check type
                           if (!is.character(Estimator.kind)) {
                               stop("Estimator.kind must be either a character or NA") 
                            }
                           # Check if passed Estimator.kind is valid
                           legal.kinds <- c("Total", "Mean",
                                            "Frequency", "Absolute Frequency", "Relative Frequency",
                                            "Ratio", "Share", "Share Ratio", "Regression Coefficient", "Quantile",
                                            "Complex Estimator")
                           if (!(Estimator.kind %in% legal.kinds)) {
                               stop("Only the following estimator kinds are allowed in a GVF model: '",
                                    paste(legal.kinds, collapse = "', '"), "'")
                            }
                        }
                    #### Checks on Resp.to.CV
                       if (!is.na(Resp.to.CV)) {
                           # Check type
                           if (!is.character(Resp.to.CV)) {
                               stop("Resp.to.CV must be either a character or NA") 
                            }
                           # Check if passed Resp.to.CV is valid
                           legal.sym <- c("resp", "Y")
                           fun.sym <- all.vars(parse(text=Resp.to.CV))
                           illegal.sym <- fun.sym[!(fun.sym %in% legal.sym)]
                           if (length(illegal.sym) > 0) {
                               stop("Only the following variables are allowed in Resp.to.CV: ",
                                    paste(legal.sym, collapse = ", "))
                            }
                        }
                    #### If nothing wrong, insert
                       n <- NROW(GVF.db[["db"]])
                       unlockBinding("GVF.db", getNamespace("ReGenesees"))
                       GVF.db$db[n + 1, "Model.id"] <<- n + 1
                       GVF.db$db[n + 1, "GVF.model"] <<- twosideform.to.char(GVF.model)
                       GVF.db$db[n + 1, "Estimator.kind"] <<- Estimator.kind
                       GVF.db$db[n + 1, "Resp.to.CV"] <<- Resp.to.CV
                       rownames(GVF.db$db) <<- NULL
                       lockBinding("GVF.db", getNamespace("ReGenesees"))
                       # Notify, if needed
                       if (isTRUE(verbose)) cat("\n# New GVF model has been registered\n\n")
                    },
     #########################################
     # Delete an *old* GVF model from the DB #
     #########################################
     delete = function(Model.id, verbose = TRUE) {
                       # Check that Model.id exists
                       if (!(Model.id %in% GVF.db$db$Model.id)) {
                             stop("Specified Model.id does not exist")
                        }
                       # Now delete
                       GVF.db.new <- GVF.db$db[GVF.db$db$Model.id != Model.id, ]
                       # Re-define Model.id column
                       n.new <- NROW(GVF.db.new)
                       if (n.new > 0) {
                           GVF.db.new$Model.id <- 1:n.new
                           rownames(GVF.db.new) <- NULL
                        }
                       unlockBinding("GVF.db", getNamespace("ReGenesees"))
                       GVF.db$db <<- GVF.db.new
                       lockBinding("GVF.db", getNamespace("ReGenesees"))
                       # Notify, if needed
                       if (isTRUE(verbose)) cat("\n# GVF model has been deleted\n\n")
                    },
     ###########################################################
     # Get the DB, so that it can be copied from the Namespace #
     ###########################################################
     get = function(verbose = TRUE) {
                     # Check if there are data, otherwise return NULL invisibly
                     if (NROW(GVF.db[["db"]]) < 1) {
                         # Notify, if needed
                         if (isTRUE(verbose)) cat("\n# GVF models db is empty\n\n")
                         return(invisible(NULL))
                        }
                     exported <- GVF.db[["db"]]
                     class(exported) <- c("GVF.db_exported", class(exported))
                     return(exported)
                    },
     ################################################################################
     # Overwrite the DB, so that one can use *its own* (previously exported) GVF DB #
     ################################################################################
     assign = function(value, verbose = TRUE) {
                       # Check class of value
                       if (!inherits(value, "GVF.db_exported")) {
                           stop("Cannot overwrite GVF models db: invalid value")
                        }
                       unlockBinding("GVF.db", getNamespace("ReGenesees"))
                       GVF.db$db <<- value
                       lockBinding("GVF.db", getNamespace("ReGenesees"))
                       # Notify, if needed
                       if (isTRUE(verbose)) cat("\n# GVF models db overwritten\n\n")
                    },
     ############################
     # Reset the startup GVF db #
     ############################
     reset = function(verbose = TRUE) {
                      unlockBinding("GVF.db", getNamespace("ReGenesees"))
                      GVF.db$db <<- GVF.db.ini
                      lockBinding("GVF.db", getNamespace("ReGenesees"))
                      # Notify, if needed
                      if (isTRUE(verbose)) cat("\n# Default GVF models db restored\n\n")
                    }
    )

class(GVF.db) <- c("GVF.db", class(GVF.db))


GVF.moldel.check <- function(GVF.model){
###################################################
# Check if a formula is elegible to be registered #
# into GVF.db.                                    #
###################################################

  # Check if GVF.model has a response term
  if (length(GVF.model) < 3) {
      stop("A GVF model must have a response term!")
    }
  # Check that GVF.model is legal (i.e. refers to meaningful variables)
  mod.vars <- all.vars(GVF.model)
  legal.vars <- c("Y", "SE", "CV", "VAR", "DEFF")
  illegal.vars <- mod.vars[!(mod.vars %in% legal.vars)]
  if (length(illegal.vars) > 0) {
      stop("Only the following variables are allowed in a GVF model: ",
           paste(legal.vars, collapse = ", "))
    }
  # Check that GVF.model response involves variances
  resp.vars <- all.vars(GVF.model[[2]])
  var.vars <- c("SE", "CV", "VAR")
  if (!any(resp.vars %in% var.vars)) {
      stop("GVF model response does not involve any of: ",
           paste(var.vars, collapse = ", "))
    }
  # Check that GVF.model explanatories include 'Y'
  pred.vars <- all.vars(GVF.model[[3]])
  if (!any(pred.vars %in% 'Y')) {
      stop("GVF model predictor must involve Y")
    }
}


print_GVF.db <- function() {
#########################
# (pseudo) print method #
#########################
  cat("\n# Registered GVF models currently available:\n\n")
  print.data.frame(GVF.db$db)
  cat("\n")
}

print.GVF.db <- function(x, ...) {
################
# print method #
################
  print_GVF.db()
}


str_GVF.db <- function() {
###################################
# (pseudo) str method             #
# NOTE: hides accessors functions #
###################################
  str(GVF.db$db)
}

str.GVF.db <- function(object, ...) {
###################################
# str method                      #
# NOTE: hides accessors functions #
###################################
  str_GVF.db()
}


dim_GVF.db <- function() {
###################################
# (pseudo) dim method             #
# NOTE: hides accessors functions #
###################################
  dim(GVF.db$db)
}

dim.GVF.db <- function(x) {
###################################
# dim method                      #
# NOTE: hides accessors functions #
###################################
  dim_GVF.db()
}


twosideform.to.char <- function(formula){
##################################################
# Turns a formula into a character, avoiding the #
# truncation problem of as.character for very    #
# long formulae.                                 #
##################################################
# remove initial white spaces generated while deparsing
char.form <- gsub("    ", "", paste(deparse(formula), collapse=""))
char.form
}

form.to.char <- function(formula){
##################################################
# Turns a formula into a character, avoiding the #
# truncation problem of as.character for very    #
# long formulae.                                 #
##################################################
# remove initial white spaces generated while deparsing
char.form <- gsub("    ", "", paste(deparse(formula), collapse=""))
# put a space after ~ (our standard)
char.form <- gsub("~", "~ ", char.form)
char.form
}


  estimator.kind <- function(stat, design) {
  #########################################################################
  # Given a survey statistic object and a design from which the former is #
  # supposed to have been derived, identify the 'precise kind' of the     #
  # estimator.                                                            #
  # Currently, possible return values (i.e. kinds) are the following:     #
  # (1)  'Total'                                                          #
  # (2)  'Absolute Frequency'                                             #
  # (3)  'Mix of Totals and Absolute Frequencies'                         #
  # (4)  'Mean'                                                           #
  # (5)  'Relative Frequency'                                             #
  # (6)  'Mix of Means and Relative Frequencies'                          #
  # (7)  'Ratio'                                                          #
  # (8)  'Share'                                                          #
  # (9)  'Share Ratio'                                                    #
  # (10) 'Regression Coefficient'                                         #
  # (11) 'Quantile'                                                       #
  # (12) 'Complex Estimator'                                              #
  # (13) 'Population Variance'                                            #
  # (14) 'Population Standard Deviation'                                  #
  #                                                                       #
  # NOTE: This is the exported version of private function Estimator.kind #
  #       (used by function gvf.input): it only has more error handling   #
  #       statements.                                                     #
  #########################################################################

    # Test design class
    if (!inherits(design, "analytic"))
        stop("Object 'design' must inherit from class analytic")

    # Test that input object stat is actually a survey statistic
    # 23/09/2015 added svystat.gr object generated by svystat
    #            when called with forGVF = FALSE
    svystat.classes <- c("svystatTM", "svystatR", "svystatS", "svystatSR", "svystatL", "svystatQ", "svystatB", "svySigma", "svySigma2", "svyby", "svystat.gr")
    ok.stat <- any(sapply(svystat.classes, function(class) inherits(stat, class)))
    if (!ok.stat){
        stop("Input object 'stat' is not a survey statistic!\n")
        }

    # Test that input statistic stat actually comes from design
    if (!inherits(stat, "svystat.gr")) {
         root.des <- as.character(attr(stat, "design"))
        }
    else {
         # 23/09/2015 if stat has class "svystat.gr", then it is a list of
         #            homogeneous survey statistics: hence I can obtain kind
         #            from its first slot
         if (!is.null(attr(stat, "group.vars"))){
             if (is.wrapped(stat)) {
                 root.des <- as.character(attr(stat[[1]], "design"))
                 stat <- stat[[1]][[1]]
                }
             else {
                 root.des <- as.character(attr(stat, "design"))
                 stat <- stat[[1]]
                }
            }
         # Code below added in 03/02/2016 for DEBUG: correctly handles the case
         # when 'stat' is generated by svystat with forGVF = FALSE *and group = NULL*
         else {
             root.des <- as.character(attr(stat[[1]], "design"))
             stat <- stat[[1]]
            }
        }
    # Here's the actual test on design
    if (!identical(root.des, character(0)) && !identical(design, get(root.des))){
        stop("Input estimates have not been computed from input design object!")
        }

    # First level assessment: what survey statistics function?
    kind <- switch(data.class(stat),
                   svystatTM = "TM", svystatTM.by = "TM",
                   svystatR  = "Ratio", svystatR.by = "Ratio",
                   svystatS  = "Share", svystatS.by = "Share",
                   svystatSR = "Share Ratio", svystatSR.by = "Share Ratio",
                   svystatB  = "Regression Coefficient",
                   svystatQ  = "Quantile", svystatQ.by = "Quantile",
                   svystatL  = "Complex Estimator", svystatL.by = "Complex Estimator",
                   svySigma2 = "Population Variance", svySigma2.by = "Population Variance",
                   svySigma  = "Population Standard Deviation", svySigma.by = "Population Standard Deviation")
    if (kind=="TM") {
        estfun <- if (!inherits(stat, "svyby")) {
                      attr(attr(stat, "origin"),"svyby")$statistic
                    }
                  else {
                      attr(stat,"svyby")$statistic
                    }

    # Second level assessment: what estimator? (i.e. split Total from Mean)
    if (estfun=="svytotal") kind <- "Total" else kind <- "Mean"

    # Third level assessment: numeric vs. factor variables, pure vs. mixed estimates
    y.vars <- attr(stat, "y.vars")

    # Test whether the variables exist in design
    vars.exist(y.vars, design)

    # Take into account variable types
    y.num  <- sapply(y.vars, function(var) is.numeric(design$variables[[var]]))

    # Treat 'ones' variables (i.e. columns of 1s) and
    # 'dummy' variables (i.e. columns of 0/1s) as if they were factors:
    # indeed they imply estimates of frequencies...
    y.od <- sapply(y.vars, function(var) is.ones(var, design) | is.dummy(var, design))

    # ...then identify REAL numeric variables
    y.Num <- y.num & !y.od

    if (kind=="Total") {
         if (all(!y.Num)) kind <- "Absolute Frequency"
         if (any(y.Num) && any(!y.Num)) kind <- "Mix of Totals and Absolute Frequencies"
        }

    if (kind=="Mean") {
         if (all(!y.Num)) kind <- "Relative Frequency"
         if (any(y.Num) && any(!y.Num)) kind <- "Mix of Means and Relative Frequencies"
        }
    }

    kind
  }
