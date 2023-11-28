gvf.input <- function(design, ..., stats = list(...)) {
####################################################################
# Given a set of survey statistic objects (or a list storing them) #
# and a design object from which the statics are supposed to have  #
# been derived, builds a dataframe storing all the data needed to  #
# fit a Generalized Variance Function (GVF) model.                 #
####################################################################

  # First verify if the function has been called inside another function:
  # this is needed to correctly manage metadata when e.g. the caller is a
  # GUI stratum
  directly <- !( length(sys.calls()) > 1 )
  design.expr <- if (directly) substitute(design)

  # Test design class
  if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")

  # Has 'stats' been passed?
  is.stats <- TRUE
  if (missing(stats)) {
      is.stats <- FALSE
    }

  # Check passed value of argument 'stats'
  if (data.class(stats)!="list") {
     warning("stats is not a list of survey statistics: will try to match the dot-dot-dot argument instead.")
     stats <- list(...)
     is.stats <- FALSE
    }
  else {
     # Strip off names of 'stats' list components
     names(stats) <- NULL
    }

  if (length(stats) < 1)
     stop("No survey estimates provided: cannot model sampling errors!")

  # Test that input objects are all survey statistics
  svystat.classes <- c("svystatTM", "svystatR", "svystatS", "svystatSR", "svystatL", "svystatQ", "svystatB", "svySigma", "svySigma2", "svyby")
  ok.stats <- sapply(stats, function(stat) any(sapply(svystat.classes, function(class) inherits(stat, class))))

  if (!all(ok.stats)) {
      where.txt <- if (is.stats) "in passed stats list" else "in dot-dot-dot argument"
      stop("All input objects must be survey statistics!\n (wrong objects positions ", where.txt,": ",
           paste(which(!ok.stats), collapse = ", "),")")
    }

  # Test that input statistics are all from the same design
  root.des <- sapply(stats, function(stat) {
                                     des <- attr(stat, "design")
                                     # if below returns NA when summary statistics
                                     # lacks a design attribute (e.g. when gvf.input
                                     # is invoked via example())
                                     if (!is.null(des)) as.character(des) else NA
                                    })
  # drop NAs (if any)
  root.des.noNA <- root.des[!is.na(root.des)]
  root.des.u <- unique(root.des.noNA)
  if (length(root.des.u) > 1) {
     root.des1 <- root.des[1]
     root.des2 <- root.des[2]
     stop("Input estimates come from different design objects!\n (e.g.: object ",
          match(root.des1, root.des)," from ", root.des1, ", ", match(root.des2, root.des)," from ", root.des2,")")
    }

  # Test that the common root design of all stats is actually what has been passed to 'design'
  # Note that first clause below skips cases where *all* stats lack the design attribute
  if ((length(root.des.u)==1) && !identical(design, get(root.des.u))){
      stop("Input estimates have not been computed from input design object!")
    }

  # Test that input statistics are elegible for a meaningful model
  stats.kind <- sapply(stats, Estimator.kind, design=design)
  stats.kind.u <- unique(stats.kind)
    # 1. Avoid mixed cases (i.e. svystatTM with both factor and numeric interest variables)
    mix <- stats.kind %in% c("Mix of Totals and Absolute Frequencies", "Mix of Means and Relative Frequencies")
    if (any(mix)) {
        stop("Some input estimates involve variables of different types, e.g. factor and numeric!\n (object positions: ",paste(which(mix), collapse = ", "),")")
    }
    # 2. Test that input statistics are all the same kind
    if (length(stats.kind.u) > 1) {
        kind1 <- stats.kind.u[1]
        kind2 <- stats.kind.u[2]
        # *Allow* (with a warning) two components mix "Complex Estimator" kind + "other kind"
        # The rationale is that svystatL can actually be used to compute exctly  "other kind"
        # of estimates! (i.e. no kind mismatch at all, actually)
        if ( (length(stats.kind.u) == 2) && any(stats.kind.u == "Complex Estimator")) {
            warning("Input estimators are not the same kind!\n (e.g.: object ",
                    match(kind1, stats.kind)," is a ", kind1, ", object ", match(kind2, stats.kind)," is a ",
                    kind2,")")
        } else {
            stop("Input estimators are not the same kind!\n (e.g.: object ",
                 match(kind1, stats.kind)," is a ", kind1, ", object ", match(kind2, stats.kind)," is a ", kind2,")")
        }
    }
    # 3. Warn if input statistics involve different variables
      # NOTE: CURRENTLY Works only for svystatTM objects (i.e. Totals, Means, Abs. and Rel. Freq.)
      #       as the rest generate NULL for ALL.y.vars below
      # NOTE: Check is meaningful for quantitative variables only: thus exludes frequencies
    ALL.y.vars <- unlist(lapply(stats, function(x) attr(x, "y.vars")))
    if (!is.null(ALL.y.vars)) {
         ALL.y.vars.u <- unique(ALL.y.vars)
         if ( (length(ALL.y.vars.u) > 1) && !(stats.kind.u %in% c("Absolute Frequency", "Relative Frequency")) ) {
             warning("Input estimates involve different interest variables: ", paste(ALL.y.vars.u, collapse = ", "))
         }
    }

    # Start building the estimates and errors dataframe
    Y    <- unlist(lapply(stats,coef))
    SE   <- unlist(lapply(stats,SE))
    CV   <- unlist(lapply(stats,cv))
    VAR  <- unlist(lapply(stats,VAR))
    name <- names(Y)
    data.to.fit <- data.frame(name = name, Y = Y, SE = SE, CV = CV, VAR = VAR)

    # Check for DEFFs
    getDEFF <- function(x) {
                 # NOTE: No deff() method available for quantiles estimators,
                 #       hence must avoid warnings that would be raised by
                 #       deff.default() 23/02/2016
                 if (data.class(x) == "svystatQ") {
                     DEFF <- NULL
                    }
                 else {
                     DEFF <- try(deff(x), silent=TRUE)
                    }
                 if (inherits(DEFF, "try-error") || is.null(DEFF)) {
                     DEFF.NA <- coef(x)
                     DEFF.NA[] <- NA
                     DEFF <- DEFF.NA
                    }
                 DEFF
                }
    DEFF <- unlist(lapply(stats, getDEFF))
    if (any(!is.na(DEFF))) {
         data.to.fit$DEFF <- DEFF
         has.Deff <- TRUE
        }
    else {
         has.Deff <- FALSE
        }

    # Drop duplicated rows that could arise from unintentionally
    # passing the same input stat more than once
    nobs.ini <- nrow(data.to.fit)
    data.to.fit <- unique(data.to.fit)
    nobs.fin <- nrow(data.to.fit)
    if (nobs.fin < nobs.ini) {
        warning("Duplicated estimates have been dropped: ", nobs.ini - nobs.fin)
    }

    # Renumber rownames
    rownames(data.to.fit) <- NULL

    # Attach useful attributes to the output dataframe
    attr(data.to.fit, "y.vars") <- if (!is.null(ALL.y.vars)) ALL.y.vars.u
    attr(data.to.fit, "stats.kind") <- stats.kind.u
    attr(data.to.fit, "has.Deff") <- has.Deff
    attr(data.to.fit, "design") <- design.expr
    class(data.to.fit) <- c("gvf.input", class(data.to.fit))

  # Return value
  data.to.fit
}


  Estimator.kind <- function(stat, design) {
  #####################################################
  # Given a survey statistic object and a design from #
  # which the former is supposed to have been derived #
  # identify the 'precise kind' of the estimator.     #
  # Currently, possible return values (i.e. kinds)    #
  # are the following:                                #
  # (1)  'Total'                                      #
  # (2)  'Absolute Frequency'                         #
  # (3)  'Mix of Totals and Absolute Frequencies'     #
  # (4)  'Mean'                                       #
  # (5)  'Relative Frequency'                         #
  # (6)  'Mix of Means and Relative Frequencies'      #
  # (7)  'Ratio'                                      #
  # (8)  'Share'                                      #
  # (9)  'Share Ratio'                                #
  # (10) 'Regression Coefficient'                     #
  # (11) 'Quantile'                                   #
  # (12) 'Complex Estimator'                          #
  # (13) 'Population Variance'                        #
  # (14) 'Population Standard Deviation'              #
  #                                                   #
  # NOTE: This is the unexported version of function  #
  #       estimator.kind, which has more error        #
  #       handling statements.                        #
  #####################################################

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


  kind.match <- function(kind, GVF.db.kind) {
  ###################################################################
  # Let 'kind' be the "stats.kind" attribute of a gvf.input object  #
  # and 'GVF.db.kind' the "Estimator.kind" of a registered model,   #
  # the function returns TRUE if the inputs are coherent, otherwise #
  # returns a diagnostic message.                                   #
  # NOTE: if 'GVF.db.kind' is NA (which is indeed a legal velue,    #
  #       see ?GVF.db), returns TRUE.                               #
  ###################################################################
    if (is.na(GVF.db.kind)) {
         return(TRUE)
        }
    if (kind == GVF.db.kind) {
         return(TRUE)
        }
    else {
         if ( (kind %in% c("Absolute Frequency", "Relative Frequency")) && (GVF.db.kind == "Frequency")) {
             return(TRUE)
            }
         else {
             warn.char <- paste("NOTE: fitted statistics have kind '", kind,
                                "'\n      whereas registered model has kind '", GVF.db.kind, "'!", sep = "")
             return(warn.char)
            }
        }
  }


  vars.exist <- function(vars, design) {
  #####################################
  # Check if variables named 'vars'   #
  # do actually exist inside 'design' #
  #####################################
    data <- design$variables
    vars <- unique(vars)
    absent.vars <- vars[!(vars %in% names(data))]
    if (length(absent.vars) > 0)
        stop("Variables not found: ", paste(absent.vars, collapse = ", "))
  }


  is.ones <- function(var, design, tol = 1E-6) {
  #######################################
  # Check if the variable named 'var'   #
  # is a column of 1s (within tolerance #
  # tol (which defaults to 1E-6).       #
  #######################################
    # vars.exist(var, design) #unnecessary: if var doesn't exist, returns FALSE!
    val <- design$variables[[var]]
    out <- FALSE
    if (is.numeric(val)) {
         test <- all.equal(sum(abs(range(val)-1)), 0, tol = tol)
         if (isTRUE(test)) {
             out <- TRUE
            }
        }
    out
  }


  is.dummy <- function(var, design, tol = 1E-6) {
  ########################################
  # Check if the variable named 'var' is #
  # a column of 1/0s (within tolerance   #
  # tol (which defaults to 1E-6).        #
  ########################################
    # vars.exist(var, design) #unnecessary: if var doesn't exist, returns FALSE!
    val <- design$variables[[var]]
    out <- FALSE
    if (is.numeric(val)) {
         test <- all.equal(sum(abs(range(val)-c(0,1))), 0, tol = tol)
         if (isTRUE(test)) {
             out <- TRUE
            }
        }
    out
  }


#`[.gvf.input` <- function(x, i, j)
################################################
# Subsetting method for 'gvf.input' objects.   #
# NOTE: The original attributes are preserved. #
################################################
#{
#    y <- `[.data.frame`(x, i, j, drop=FALSE)
#    attr(y, "y.vars") <- attr(x, "y.vars")
#    attr(y, "stats.kind") <- attr(x, "stats.kind")
#    attr(y, "has.Deff") <- attr(x, "has.Deff")
#    attr(y, "design") <- attr(y, "design")
#    class(y) <- class(x)
#    y
#}


plot.gvf.input <- function(x, ...){
##############################################
# Plot method on gvf.input objects: yields a #
# scatterplot matrix + polynomial smoother.  #
##############################################
  xx <- as.data.frame(x[,-1])
  plot(xx, panel=panel.smooth, ...)
# OLD CODE
# graphics:::plot.data.frame(x[,-1], panel=panel.smooth)
# This would avoid using ':::' but unfortunately can't manage to avoid plotting
# the first column...
# NextMethod("plot", object = x[,-1], panel=panel.smooth)
}


fit.gvf <- function(gvf.input, model = NULL, weights = NULL)
#################################################
# Generic function. Fit one or many GVF models. #
#################################################
{
UseMethod("fit.gvf")
}

fit.gvf.default <- function(gvf.input, model = NULL, weights = NULL)
####################################
# Default method: raises an error. #
####################################
{
stop("Object 'gvf.input' must inherit from class gvf.input or gvf.input.gr!")
}

fit.gvf.gvf.input <- function(gvf.input, model = NULL, weights = NULL) {
################################################################################
# This function fits one or more GVF models to a set of survey statistics.     #
# The primary purpose of fitting multiple models is for comparison: the user   #
# is expected to eventually choose his preferred model, in order to obtain     #
# sampling errors predictions.                                                 #
#                                                                              #
# Parameter 'model' can be either:                                             #
# 1) NULL (the default) meaning *all* registered models in GVF.db              #
# or:                                                                          #
# 2) any sub-vector of GVF.db$db$Model.id (i.e. an arbitrary selection of      #
#    registered models)                                                        #
# or:                                                                          #
# 3) an arbitrary (single) formula (i.e. any custom, user-defined model).      #
#                                                                              #
# NOTE: The class of the output is 'gvf.fit' if just a single model has been   #
#       fit, otherwise it is 'gvf.fits'. Objects of class 'gvf.fits'           #
#       essentially behave as lists of objects of class 'gvf.fit'.             #
################################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 2L )
gvf.input.expr <- if (directly) substitute(gvf.input)

# Test gvf.input class
if (!inherits(gvf.input, "gvf.input")) 
     stop("Object 'gvf.input' must inherit from class gvf.input")

# Check model type:
if (is.null(model)) {
     # i.e. *all* registered models
     model <- GVF.db$db$Model.id
    }
if (!inherits(model, "formula") && !is.numeric(model)) {
     stop("GVF models can be specified only as formulae or by passing the numeric ids of registered models! (see ?GVF.db)")
    }

#!!!!!! MANDATORY, otherwise cannot pass weights when calling lm
this.env <- environment()

  # 1) numeric case (i.e. registered models)
if (is.numeric(model)) {
     model.id <- sort(intersect(model, GVF.db$db$Model.id))
     if (length(model.id) < 1) {
         stop("Specified model doesn't match any available registered GVF! (see ?GVF.db)")
        }
     model.char <- GVF.db$db[GVF.db$db$Model.id %in% model.id, "GVF.model"]
     estkind.char <- GVF.db$db[GVF.db$db$Model.id %in% model.id, "Estimator.kind"]
     model.form <- lapply(model.char, FUN = as.formula, env = this.env)
     if (length(model.form) == 1) {
         model.form <- model.form[[1]]
        }
     # must check estimator kinds coherence? YES
     estkind.check <- TRUE
    }

  # 2) formula case (i.e. a custom model)
if (inherits(model, "formula")) {
     ## Set model.id to NULL (i.e. a custom model is being fit)
     model.id <- NULL

     # Run all checks on GVF.model at once
     GVF.moldel.check(model)

     model.form <- model
     #!!!!!! MANDATORY, otherwise cannot pass weights when calling lm
     environment(model.form) <- this.env
     # must check estimator kinds coherence? NO
     estkind.check <- FALSE
    }

# Check weights:
if (is.null(weights)) {
     w.char <- "NULL"
     weights <- NULL
    }
else {
     if (!inherits(weights,"formula"))
         stop("If supplied, weights must be passed as a formula")
     if (length(weights) >= 3)
         stop("weights formula must be one-sided!")
     w.char <- form.to.char(weights)
     w.form <- as.formula(paste(w.char, "- 1"))
     mf <- model.frame(w.form, gvf.input, na.action = na.pass)
     weights <- model.matrix(w.form, mf)
    }

# Fit the provided gvf model(s)
if (is.list(model.form)) {
     elm <- vector(mode = "list", length = length(model.form))
     i.lm <- 0
     for (form in model.form) {
          i.lm <- i.lm + 1
          elm[[i.lm]] <- lm(formula = form, data = gvf.input, weights = weights,
                            na.action = na.exclude, model = TRUE,
                            x = TRUE, y = TRUE)
          # Re-introduce NAs (if any) in response, as they have been dropped from $y slot
          if (!is.null(nas <- elm[[i.lm]]$na.action)) {
              y <- elm[[i.lm]]$y
              resp <- rep(NA, length(y) + length(nas))
              resp[-nas] <- y
              resp.names <- rep(NA, length(y) + length(nas))
              resp.names[-nas] <- names(y)
              resp.names[nas]  <- names(nas)
              names(resp) <- resp.names
              elm[[i.lm]]$resp <- resp
            }
          else {
              elm[[i.lm]]$resp <- elm[[i.lm]]$y
            }
          # Check estimator kind consistency
          chk <- kind.match(attr(gvf.input, "stats.kind"), estkind.char[i.lm])
          # Set attributes
          attr(elm[[i.lm]], "model") <- twosideform.to.char(form)
          attr(elm[[i.lm]], "model.id") <- model.id[i.lm]
          attr(elm[[i.lm]], "kind.mismatch") <- if (!isTRUE(chk)) chk
          attr(elm[[i.lm]], "weights") <- w.char
          attr(elm[[i.lm]], "gvf.input.expr") <- gvf.input.expr
          class(elm[[i.lm]]) <- c("gvf.fit", class(elm[[i.lm]]))
        }
    }
else {
      elm <- lm(formula = model.form, data = gvf.input, weights = weights,
                na.action = na.exclude, model = TRUE,
                x = TRUE, y = TRUE)
      # Re-introduce NAs (if any) in response, as they have been dropped from $y slot
      if (!is.null(nas <- elm$na.action)) {
          y <- elm$y
          resp <- rep(NA, length(y) + length(nas))
          resp[-nas] <- y
          resp.names <- rep(NA, length(y) + length(nas))
          resp.names[-nas] <- names(y)
          resp.names[nas]  <- names(nas)
          names(resp) <- resp.names
          elm$resp <- resp
        }
      else {
          elm$resp <- elm$y
        }
      # If needed, check estimator kind consistency
      if (estkind.check) {
          chk <- kind.match(attr(gvf.input, "stats.kind"), estkind.char)
        }
      # Set attributes
      attr(elm, "model") <- twosideform.to.char(model.form)
      attr(elm, "model.id") <- model.id
      attr(elm, "kind.mismatch") <- if (estkind.check && !isTRUE(chk)) chk
      attr(elm, "weights") <- w.char
      attr(elm, "gvf.input.expr") <- gvf.input.expr
      class(elm) <- c("gvf.fit", class(elm))
    }

# Set other attributes
attr(elm, "gvf.input") <- gvf.input

if (inherits(elm, "list")) {
     class(elm) <- c("gvf.fits", class(elm))
    }
# NOTE: if output inherits from list, then it is a collection of 'gvf.fit'
# objects, i.e. a 'gvf.fits' object
elm
}


`[.gvf.fits` <- function(x, ...)
#############################################
# Subsetting method for 'gvf.fits' objects. #
# NOTE: The original 'gvf.input' attributes #
#       are preserved.                      #
#############################################
{
    y <- NextMethod("[")
    attr(y, "gvf.input") <- attr(x, "gvf.input")
    attr(y, "gvf.input.expr") <- attr(x, "gvf.input.expr")
    class(y) <- class(x)
    y
}


`[[.gvf.fits` <- function(x, ...)
##################################################
# Extract method for 'gvf.fits' objects. Recall  #
# that components of x are actually 'gvf.fit'    #
# objects.                                       #
# NOTE: The original 'gvf.input' attributes is   #
#       preserved.                               #
##################################################
{
    y <- NextMethod("[[")
    attr(y, "gvf.input") <- attr(x, "gvf.input")
    y
}


print.gvf.fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
#######################################
# Print method for 'gvf.fit' objects. #
#######################################
{
    # Vertical separator
    vsep <- paste(rep("-", 0.5 * getOption("width")), collapse = "")

    cat(vsep, "\n")
    cat("GVF model:  ", attr(x, "model"),   "\n")
    if (!is.null(attr(x, "model.id"))) {
         cat("- model.id: ", attr(x, "model.id"), "\n")
        }
    if (!is.null(attr(x, "gvf.input.expr"))) {
         cat("- data:     ", deparse(attr(x, "gvf.input.expr")), "\n")
        }
    cat("- weights:  ", attr(x, "weights"), "\n")
    cat("\n")

    if (length(coef(x))) {
         cat("Coefficients:\n")
         print.default(format(coef(x), digits = digits), print.gap = 2L, 
                       quote = FALSE)
        }
    else {
         cat("No coefficients\n")
        }

    if (!is.null(attr(x, "kind.mismatch"))) {
         cat("\n")
         cat(attr(x, "kind.mismatch"), "\n")
        }

    cat("\n\n")
    invisible(x)
}


print.gvf.fits <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
########################################
# Print method for 'gvf.fits' objects. #
########################################
{
    for (e in x) {
        # print.gvf.fit(e, digits = digits, ...)
         print(e, digits = digits, ...) # PROVA!
        }
    invisible(x)
}


summary.gvf.fit <- function(object, correlation = FALSE,
                            symbolic.cor = FALSE, ...)
#########################################
# Summary method for 'gvf.fit' objects. #
#########################################
{
    s <- summary.lm(object, correlation = correlation,
                    symbolic.cor = symbolic.cor, ...)
    attr(s, "model")          <- attr(object, "model")
    attr(s, "model.id")       <- attr(object, "model.id")
    attr(s, "kind.mismatch")  <- attr(object, "kind.mismatch")
    attr(s, "weights")        <- attr(object, "weights")
    attr(s, "gvf.input.expr") <- attr(object, "gvf.input.expr")
    class(s) <- c("summary.gvf.fit", class(s))
    s
}


summary.gvf.fits <- function(object, correlation = FALSE,
                             symbolic.cor = FALSE, ...)
##########################################
# Summary method for 'gvf.fits' objects. #
##########################################
{
    ss <- vector(mode = "list", length = length(object))
    i.ss <- 0
    for (e in object) {
         i.ss <- i.ss + 1
        # ss[[i.ss]] <- summary.gvf.fit(e, correlation = correlation,
        #                               symbolic.cor = symbolic.cor, ...)
         ss[[i.ss]] <- summary(e, correlation = correlation,
                               symbolic.cor = symbolic.cor, ...) # PROVA!
        }
    class(ss) <- c("summary.gvf.fits", class(ss))
    ss
}


print.summary.gvf.fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  symbolic.cor = x$symbolic.cor, 
                                  signif.stars = getOption("show.signif.stars"),
                                  ...)
##################################################
# Print method for summary on 'gvf.fit' objects. #
##################################################
{
    # Vertical separator
    vsep <- paste(rep("-", 0.6 * getOption("width")), collapse = "")

    cat(vsep, "\n")
    cat("GVF model:  ", attr(x, "model"),   "\n")
    if (!is.null(attr(x, "model.id"))) {
         cat("- model.id: ", attr(x, "model.id"), "\n")
        }
    if (!is.null(attr(x, "gvf.input.expr"))) {
         cat("- data:     ", deparse(attr(x, "gvf.input.expr")), "\n")
        }
    cat("- weights:  ", attr(x, "weights"), "\n")
    cat("\n")

    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if (!is.null(x$weights) && diff(range(x$weights))) 
        "Weighted ", "Residuals:\n", sep = "")
    if (rdf > 5L) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2L) 
            structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                dimnames(resid)[[2L]]))
        else {
            zz <- zapsmall(quantile(resid), digits + 1L)
            structure(zz, names = nam)
        }
        print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
        print(resid, digits = digits, ...)
    }
    else {
        cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
        cat("\n")
    }
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    }
    else {
        if (nsingular <- df[3L] - df[1L]) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\nResidual standard error:", format(signif(x$sigma, 
        digits)), "on", rdf, "degrees of freedom")
    cat("\n")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
        cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, 
            digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L], 
            digits = digits), "on", x$fstatistic[2L], "and", 
            x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L], 
                x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), 
                digits = digits))
        cat("\n")
    }
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }

    if (!is.null(attr(x, "kind.mismatch"))) {
         cat("\n")
         cat(attr(x, "kind.mismatch"), "\n")
        }

    cat("\n\n")
    invisible(x)
}


print.summary.gvf.fits <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   symbolic.cor = x$symbolic.cor, 
                                   signif.stars = getOption("show.signif.stars"),
                                   ...)
###################################################
# Print method for summary on 'gvf.fits' objects. #
###################################################
{
    for (e in x) {
        # print.summary.gvf.fit(e, digits = digits, symbolic.cor = symbolic.cor,
        #                       signif.stars = signif.stars, ...)
         print(e, digits = digits, symbolic.cor = symbolic.cor,
               signif.stars = signif.stars, ...) # PROVA!
        }
    invisible(x)
}


plot.gvf.fit <- function(x, which.more = 1:3, id.n = 3,
                         labels.id = names(residuals(x)), cex.id = 0.75,
                         label.pos = c(4, 2), cex.caption = 1, Main = NULL, ...)
#####################################################
# Plot method for 'gvf.fit' objects.                #
# NOTE: which.more = NULL means that only the       #
#       'Observed vs Fitted' plot will be provided. #
#####################################################
{
    if (!is.null(which.more)) {
         if (!is.numeric(which.more) || any(which.more < 1) || any(which.more > 6) || (length(which.more) > 3)) {
             stop("If specified, 'which.more' must be in 1:6 and of length <= 3!")
            }
         # How many plots?
         np <- 1 + length(unique(which.more))
         # Set plot region according to the number of plots requested
         old.par <- if (np == 4) par(mfrow = c(2, 2)) else par(mfrow = c(1, np))
        }

    ### Plot1: Observed vs Fitted - START ###
    r <- residuals(x)
    yh <- predict(x)
    yobs <- x$resp
    w <- weights(x)
    if (!is.null(w)) {
         wind <- w != 0
         r <- r[wind]
         yh <- yh[wind]
         yobs <- yobs[wind]
         w <- w[wind]
         labels.id <- labels.id[wind]
        }
    n <- length(r)
    if (is.null(id.n)) {
         id.n <- 0
        }
    else {
         id.n <- as.integer(id.n)
         if (id.n < 0L || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
        }
    if (id.n > 0L) {
         if (is.null(labels.id)) 
             labels.id <- paste(1L:n)
         iid <- 1L:id.n
         show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
         text.id <- function(x, y, ind, adj.x = TRUE) {
                             labpos <- if (adj.x) {
                                         label.pos[1 + as.numeric(x > mean(range(x)))]
                                        }
                                       else 3
                             text(x, y, labels.id[ind], cex = cex.id,
                                  xpd = TRUE, pos = labpos, offset = 0.25)
                            }
        }
    ylim <- range(yobs, na.rm = TRUE)
    if (id.n > 0) {
         ylim <- extendrange(r = ylim, f = 0.08)
        }
    dev.hold()
    plot(yh, yobs, xlab = "Fitted values", ylab = "Observed values",
         ylim = ylim, main = Main, ...)
    mtext("Observed vs Fitted", 3, 0.25, cex = cex.caption)
    if (id.n > 0) {
         y.id <- yobs[show.r]
         y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
         text.id(yh[show.r], y.id, show.r)
        }
    abline(0:1, col = "red")
    dev.flush()
    ### Plot1: Observed vs Fitted - END ###

    if (!is.null(which.more)) {
         stats_plot.lm(x, which = which.more, id.n = id.n, labels.id = labels.id,
                       cex.id = cex.id, label.pos = label.pos,
                       cex.caption = cex.caption, main = Main, ...)
         par(old.par)
        }
}


plot.gvf.fits <- function(x, which.more = NULL, id.n = 3,
                          labels.id = names(residuals(x)), cex.id = 0.75,
                          label.pos = c(4, 2), cex.caption = 1, Main = NULL, ...)
#########################################################
# Plot method for 'gvf.fits' objects.                   #
# NOTE: The default value of which.more differ from the #
#       one adopted for class gvf.fit and provides just #
#       the 'Observed vs Fitted' plot.                  #
#########################################################
{
    # How many models?
    nm <- length(x)
    # PROVA! Must finalize this function...
    if (nm > 1) {
         oask <- devAskNewPage(TRUE)
         on.exit(devAskNewPage(oask))
        }
    for (m in x) {
         Model.lab <- paste("Model: ", attr(m, "model"), sep = "")
         main <- if (is.null(Main)) Model.lab else paste(Main, " - ", Model.lab, sep = "")

         plot(m, which.more = which.more, id.n = id.n, labels.id = labels.id,
              cex.id = cex.id, label.pos = label.pos, cex.caption = cex.caption,
              Main = main, ...)
        }
}


###############
# AIC methods #
###############
# No need of AIC.gvf.fit
AIC.gvf.fits <- function(object, ...)
{
    AA <- vector(mode = "numeric", length = length(object))
    i.AA <- 0
    for (e in object) {
         i.AA <- i.AA + 1
         AA[i.AA] <- AIC(e, ...)
        }
    AA
}


###############
# BIC methods #
###############
# No need of BIC.gvf.fit
BIC.gvf.fits <- function(object, ...)
{
    BB <- vector(mode = "numeric", length = length(object))
    i.BB <- 0
    for (e in object) {
         i.BB <- i.BB + 1
         BB[i.BB] <- BIC(e, ...)
        }
    BB
}


getR2 <- function (object, adjusted = FALSE, ...){
##############################################
# Generic function. Get R2 (adjusted or not) #
# of fitted GVF models.                      #
##############################################
UseMethod("getR2")
}

getR2.default <- function(object, adjusted = FALSE, ...)
####################################
# Default method: raises an error. #
####################################
{
stop("Input object is not a fitted GVF model!")
}

getR2.gvf.fit <- function(object,
                          adjusted = FALSE, ...)
###################
# gvf.fit method. #
###################
{
s <- summary(object)
if (isTRUE(adjusted)) {
     r <- s$adj.r.squared
    }
else {
     r <- s$r.squared
    }
r
}

getR2.gvf.fits <- function(object,
                           adjusted = FALSE, ...)
####################
# gvf.fits method. #
####################
{
    rr <- vector(mode = "numeric", length = length(object))
    i.rr <- 0
    for (e in object) {
         i.rr <- i.rr + 1
         rr[i.rr] <- getR2(e, adjusted = adjusted, ...)
        }
    rr
}


drop.gvf.points <- function(x, method = c("pick", "cut"), which.plot = 1:2,
                            res.type = c("standard", "student"), res.cut = 3,
                            id.n = 3, labels.id = NULL,
                            cex.id = 0.75, label.pos = c(4, 2),
                            cex.caption = 1, col = NULL, drop.col = "red",
                            ...)
#############################################
# Generic function. Drops observations from #
# a fitted GVF model and refits the model.  #
#############################################
{
UseMethod("drop.gvf.points")
}

drop.gvf.points.default <- function(x, method = c("pick", "cut"), which.plot = 1:2,
                                    res.type = c("standard", "student"), res.cut = 3,
                                    id.n = 3, labels.id = NULL,
                                    cex.id = 0.75, label.pos = c(4, 2),
                                    cex.caption = 1, col = NULL, drop.col = "red",
                                    Main = NULL, ...)
####################################
# Default method: raises an error. #
####################################
{
stop("Input object must be a single fitted GVF model! (i.e. of class gvf.fit or gvf.fit.gr)")
}

drop.gvf.points.gvf.fit <- function(x, method = c("pick", "cut"), which.plot = 1:2,
                                    res.type = c("standard", "student"), res.cut = 3,
                                    id.n = 3, labels.id = NULL,
                                    cex.id = 0.75, label.pos = c(4, 2),
                                    cex.caption = 1, col = NULL, drop.col = "red",
                                    Main = NULL, ...)
##################################################################################
# This function drops observations from a fitted GVF model and refits the model. #
# Method for class gvf.fit.                                                      #
#                                                                                #
# If method = "pick", observations to be dropped are identified by clicking on   #
# points of a plot. Argumemt 'which.plot' determines the nature of the plot:     #
# value 1 is for 'Observed vs Fitted', value 2 is for 'Residuals vs Fitted'. In  #
# the latter case, argument 'res.type' specifies what kind of residuals have to  #
# be plotted. Argument 'id.n' specifies how many points have to be labelled      #
# initially, starting with the most extreme in terms of the selected residuals:  #
# this applies to both kinds of plots.                                           #
#                                                                                #
# If method = "cut", observations to be dropped are those with residuals whose   #
# absolute value exceeds the value of argument 'res.cut'. Again, argument        #
# 'res.type' specifies what kind of residuals have to be used. The points which  #
# have been cut are highlighted on a plot, whose nature is again specified by    #
# argument 'which.plot'. If which.plot = 1:2, dropped points are visualized on   #
# both the 'Observed vs Fitted' and the 'Residuals vs Fitted' graphs             #
# simultaneously.                                                                #
#                                                                                #
# All the other arguments have the same meaning and function as in plot.lm       #
##################################################################################
{
    # Check on x: must be a single fitted GVF model
    if (!inherits(x, "gvf.fit")) stop("Input object must be a single fitted GVF model!")

    method <- match.arg(method)
    if (method == "pick") {
         if (!missing(which.plot)) {
             if (!is.numeric(which.plot) || any(which.plot < 1) || any(which.plot > 2) || (length(which.plot) > 1)) {
                 stop("For method = 'pick', 'which.plot' must be either 1 or 2! (1 means 'Observed vs Fitted', 2 means 'Residuals vs Fitted')")
                }
            }
         else {
             which.plot <- which.plot[1]
            }
         # Trigger graphic device to prompt user:
         oask <- devAskNewPage(TRUE)
         on.exit(devAskNewPage(oask))
    }
    if (method == "cut") {
         if (!is.numeric(which.plot) || any(which.plot < 1) || any(which.plot > 2) || (length(which.plot) > 2)) {
             stop("For method = 'cut', 'which.plot' must be in 1:2! (1 means 'Observed vs Fitted', 2 means 'Residuals vs Fitted')")
            }
         which.plot <- unique(which.plot)
    }

    res.type <- match.arg(res.type)
    get.residuals <- switch(res.type,  "standard" = rstandard, "student" = rstudent)
    if (!is.numeric(res.cut) || (res.cut < 0)) stop("Argument 'res.cut' must be positive!")

    # How many and which plots?
    if (length(which.plot) == 2) {
         old.par <- par(mfrow = c(1, 2))
        }
    show <- rep(FALSE, 2)
    show[which.plot] <- TRUE

    # Get initial color for points
    if (is.null(col)) col <- par("col")

    if (show[1L]) {
    ### Plot1: Observed vs Fitted - START ###
#         r <- residuals(x) !OLD CODE USED rough RESIDUALS!
         r <- get.residuals(x)
         yh <- predict(x)
         yobs <- x$resp
         w <- weights(x)
         if (!is.null(w)) {
             wind <- w != 0
             r <- r[wind]
             yh <- yh[wind]
             yobs <- yobs[wind]
             w <- w[wind]
             labels.id <- labels.id[wind]
            }
         n <- length(r)

         if (method == "cut") {
             # Observations to cut
             id.n <- sum(abs(r) > res.cut, na.rm = TRUE)
            }
         if (is.null(id.n)) {
             id.n <- 0
            }
         else {
             id.n <- as.integer(id.n)
             if (id.n < 0L || id.n > n) 
                 stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
            }
         if (id.n > 0L) {
             if (is.null(labels.id)) 
                 labels.id <- paste(1L:n)
                 iid <- 1L:id.n
                 show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
                 text.id <- function(x, y, ind, adj.x = TRUE, ...) {
                                     labpos <- if (adj.x) {
                                                     label.pos[1 + as.numeric(x > mean(range(x)))]
                                                    }
                                               else 3
                                     text(x, y, labels.id[ind], cex = cex.id,
                                          xpd = TRUE, pos = labpos, offset = 0.25, ...)
                                    }
            }
         ylim <- range(yobs, na.rm = TRUE)
         if (id.n > 0) {
             ylim <- extendrange(r = ylim, f = 0.08)
            }
         dev.hold()
         # Set new color for cut points (if any)
         plot.col <- rep_len(col, n)
         if ( (method == "cut") && (id.n > 0) ) {
             plot.col[show.r] <- drop.col
            }
         plot(yh, yobs, xlab = "Fitted values", ylab = "Observed values",
              ylim = ylim, col = plot.col, main = Main, ...)
         mtext("Observed vs Fitted", 3, 0.25, cex = cex.caption)
         if (id.n > 0) {
             y.id <- yobs[show.r]
             y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
             text.col <- if (method == "cut") drop.col
             text.id(yh[show.r], y.id, show.r, col = text.col)
            }
         abline(0:1, col = "red")
         dev.flush()
    ### Plot1: Observed vs Fitted - END ###
    }

    if (show[2L]) {
    ### Plot2: Residuals vs Fitted - START ###
         res.string <- switch(res.type,  "standard" = "Standardized Residuals", "student" = "Studentized Residuals")
         r <- get.residuals(x)
         yh <- predict(x)
         w <- weights(x)
         if (!is.null(w)) {
             wind <- w != 0
             r <- r[wind]
             yh <- yh[wind]
             w <- w[wind]
             labels.id <- labels.id[wind]
            }
         n <- length(r)

         if (method == "cut") {
             # Observations to cut
             id.n <- sum(abs(r) > res.cut, na.rm = TRUE)
            }
         if (is.null(id.n)) {
             id.n <- 0
            }
         else {
             id.n <- as.integer(id.n)
             if (id.n < 0L || id.n > n) 
                 stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
            }
         if (id.n > 0L) {
             if (is.null(labels.id)) 
                 labels.id <- paste(1L:n)
                 iid <- 1L:id.n
                 show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
                 text.id <- function(x, y, ind, adj.x = TRUE, ...) {
                                     labpos <- if (adj.x) {
                                                     label.pos[1 + as.numeric(x > mean(range(x)))]
                                                    }
                                               else 3
                                     text(x, y, labels.id[ind], cex = cex.id,
                                          xpd = TRUE, pos = labpos, offset = 0.25, ...)
                                    }
            }
         ylim <- range(r, na.rm = TRUE)
         if (id.n > 0) {
             ylim <- extendrange(r = ylim, f = 0.08)
            }
         dev.hold()
         # Set new color for cut points (if any)
         plot.col <- rep_len(col, n)
         if ( (method == "cut") && (id.n > 0) ) {
             plot.col[show.r] <- drop.col
            }
         plot(yh, r, xlab = "Fitted values", ylab = res.string,
              ylim = ylim, col = plot.col, main = Main, ...)
         panel.smooth(yh, r, col = plot.col, ...)
         mtext(paste(res.string, " vs Fitted", sep = ""), 3, 0.25, cex = cex.caption)
         if (id.n > 0) {
             y.id <- r[show.r]
             y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
             text.col <- if (method == "cut") drop.col
             text.id(yh[show.r], y.id, show.r, col = text.col)
            }
         abline(h = 0, lty = 3, col = "gray")
         if (method == "cut") {
             abline(h =   res.cut, lty = 2, col = drop.col)
             abline(h = - res.cut, lty = 2, col = drop.col)
            }
         dev.flush()
    ### Plot2: Residuals vs Fitted - END ###
    }


    if (method == "pick") {
    # Identify points by mouse click...
         # What variable was there on the y axis?
         if (which.plot == 1) yy <- yobs
         if (which.plot == 2) yy <- r
         ii <- identify(yh, yy, labels = labels.id, 
                        cex = cex.id, col = drop.col)
        }
    if (method == "cut") {
    # Identify points by residuals threshold...
         if (id.n > 0L) { 
             ii <- show.r
            }
         else {
             ii <- NULL
            }
         # If needed, re-set graphical parameters
         if (length(which.plot) == 2) par(old.par)        
        }


    if ( (n.ii <- length(ii)) > 0) {
         # Now subset gvf.input by dropping identified points...
         gvf.input <- attr(x, "gvf.input")
         gvf.input <- gvf.input[-ii, ]

         # Lastly re-fit the obtained gvf.input with the same model...
         model <- attr(x, "model")
         model.form <- as.formula(model)
         w.char <- attr(x, "weights")
         if (w.char != "NULL"){
             weights <- as.formula(w.char)
            }
         else {
             weights <- NULL
            }

         # Was x a registered GVF model?
         old.model.id <- attr(x, "model.id")
         if (!is.null(old.model.id)) {
             # Is there still in GVF.db an entry with the same Model.id?
             if (any(GVF.db$db$Model.id == old.model.id)) {
                 now.model <- GVF.db$db$GVF.model[GVF.db$db$Model.id == old.model.id]
                 # Check if GVF.db entry with the old model id has been changed
                 # since the time the old model x was fitted
                 if (!identical(model, now.model)) {
                     warning("Input GVF fitted model does not match the corresponding entry in GVF models db!\n (perhaps you changed GVF.db after fitting?)")
                     # drop useless Model.id info
                     old.model.id <- NULL
                    }
                }
             else {
                 # No GVF.db entry with same Model.id: drop useless Model.id info
                 old.model.id <- NULL
                }
            }

         # Use GVF Model.id info when meaningful, else re-fit using the model formula
         if (!is.null(old.model.id)) {
             out <- fit.gvf(gvf.input = gvf.input, model = old.model.id, weights = weights)
            }
         else {
             out <- fit.gvf(gvf.input = gvf.input, model = model.form, weights = weights)
            }
         # Copy 'gvf.input.expr' attribute from x, as out would loose it being fit.gvf *not* invoked directly
         attr(out, "gvf.input.expr") <- attr(x, "gvf.input.expr")

           # Notify how many cases have been dropped...
           cat("\n# GVF model has been re-fitted after having dropped ", n.ii," observations", sep="")
           cat("\n#")
           # ...the corresponding *R2* variation
           cat("\n# - R^2             passed from ", getR2(x)," to ", getR2(out), sep="")
           # ...the corresponding *adjusted R2* variation
           cat("\n# - Adjusted R^2    passed from ", getR2(x, adjusted = TRUE)," to ", getR2(out, adjusted = TRUE),
               sep="")
           # ...the corresponding *AIC* variation
           cat("\n# - AIC             passed from ", AIC(x)," to ", AIC(out), sep="")
           # ...the corresponding adjusted *BIC* variation
           cat("\n# - BIC             passed from ", BIC(x)," to ", BIC(out), "\n\n", sep="")
         out
        }
     else {
         x
        }
}


predictCV <- function (object, new.Y = NULL,
                       scale = NULL, df = Inf,
                       interval = c("none", "confidence", "prediction"),
                       level = 0.95, na.action = na.pass,
                       pred.var = NULL, weights = 1)
#######################################
# Generic function. Predict CV values #
# on the basis of a fitted GVF model. #
#######################################
{
UseMethod("predictCV")
}

predictCV.default <- function(object, new.Y = NULL,
                              scale = NULL, df = Inf,
                              interval = c("none", "confidence", "prediction"),
                              level = 0.95, na.action = na.pass,
                              pred.var = NULL, weights = 1)
####################################
# Default method: raises an error. #
####################################
{
stop("Input object is not a fitted GVF model!")
}

predictCV.gvf.fit <- function(object, new.Y = NULL,
                              scale = NULL, df = Inf,
                              interval = c("none", "confidence", "prediction"),
                              level = 0.95, na.action = na.pass,
                              pred.var = NULL, weights = 1)
###################################
# gvf.fit method.                 #
###################################
{
# Checks on input object (which is a *single* fitted GVF model)
  # Get GVF model id
  model.id <- attr(object, "model.id")
  ## 1. must be a *registered* fitted GVF model
  if (is.null(model.id)){
      stop("Input fitted model not registered in GVF models db: cannot predict CVs!")
    }
  # Get appropriate resp.to.cv transformation
  resp.to.cv <- GVF.db$db$Resp.to.CV[GVF.db$db$Model.id == model.id]
  ## 2. must have a *non-missing* Resp.to.CV entry in GVF archive
  if (is.na(resp.to.cv)){
      stop("Input fitted model has missing Resp.to.CV in GVF models db: cannot predict CVs!")
    }

# By default new.Y is the fitted gvf.input object. This will return the CVs
# obtained by transforming the fitted response values
if (is.null(new.Y)){
     new.Y <- attr(object, "gvf.input")
    }

# Build resp.to.cv expression
expr <- parse(text = resp.to.cv)

# Compute predictions for the GVF model response:
if (!is.null(pred.var)){
     # Passed variance(s) for prediction intervals around future observations
     p <- predict.lm(object = object, newdata = new.Y,
                     scale = scale, df = df,
                     interval = interval, level = level,
                     type = "response", terms = NULL, na.action = na.action,
                     pred.var = pred.var, weights = weights)
    }
else {
     # This assumes the same default for 'pred.var' as in predict.lm()
     p <- predict.lm(object = object, newdata = new.Y,
                     scale = scale, df = df,
                     interval = interval, level = level,
                     type = "response", terms = NULL, na.action = na.action,
                     weights = weights)
    }

# Handle possible p types (which depend on interval)
if (is.matrix(p)){
     p <- as.data.frame(p)
    }
else {
     p <- data.frame(fit = p)
    }

# Start building actual CV predictions
p.CV <- p
colnames(p.CV) <- paste("CV", colnames(p), sep=".")
# Evaluate resp.to.cv expression (expr)
i.col <- 0
for (cols in p){
     i.col <- i.col + 1
     new.Y.cols <- data.frame(new.Y, resp = cols)
     p.CV[, i.col] <- eval(expr, envir = new.Y.cols)
}
# Bind CV predictions to the new.Y input data
p.CV <- data.frame(new.Y, p.CV)
# Return output object
return(p.CV)
}

predictCV.gvf.fits <- function(object, new.Y = NULL,
                               scale = NULL, df = Inf,
                               interval = c("none", "confidence", "prediction"),
                               level = 0.95, na.action = na.pass,
                               pred.var = NULL, weights = 1)
###################################
# gvf.fits method.                #
###################################
{
    pp <- vector(mode = "list", length = length(object))
    for (i.pp in 1:length(object)) {
         pp[[i.pp]] <- predictCV(object[[i.pp]], new.Y = new.Y,
                                 scale = scale, df = df,
                                 interval = interval,
                                 level = level, na.action = na.action,
                                 pred.var = pred.var, weights = weights)
        }
    pp
}


###################
# predict methods #
###################
predict.gvf.fit <- function(object, ...)
{
stats::predict.lm(object, ...)
}

predict.gvf.fits <- function(object, ...)
{
    pp <- vector(mode = "list", length = length(object))
    i.pp <- 0
    for (e in object) {
         i.pp <- i.pp + 1
         pp[[i.pp]] <- predict(e, ...)
        }
    pp
}


#####################
# rstandard methods #
#####################
rstandard.gvf.fit <- function(model, ...)
{
stats_rstandard.lm(model, ...)
}

rstandard.gvf.fits <- function(model, ...)
{
    pp <- vector(mode = "list", length = length(model))
    i.pp <- 0
    for (e in model) {
         i.pp <- i.pp + 1
         pp[[i.pp]] <- rstandard(e, ...)
        }
    pp
}


####################
# rstudent methods #
####################
rstudent.gvf.fit <- function(model, ...)
{
stats_rstudent.lm(model, ...)
}

rstudent.gvf.fits <- function(model, ...)
{
    pp <- vector(mode = "list", length = length(model))
    i.pp <- 0
    for (e in model) {
         i.pp <- i.pp + 1
         pp[[i.pp]] <- rstudent(e, ...)
        }
    pp
}


#################
# anova methods #
#################
anova.gvf.fit <- function(object, ...)
{
stats_anova.lm(object, ...)
}

anova.gvf.fits <- function(object, ...)
{
    aa <- vector(mode = "list", length = length(object))
    i.aa <- 0
    for (e in object) {
         i.aa <- i.aa + 1
         aa[[i.aa]] <- anova(e, ...)
        }
    aa
}


################
# coef methods #
################
coef.gvf.fit <- function(object, ...)
{
stats_coef.default(object, ...)
}

coef.gvf.fits <- function(object, ...)
{
    cc <- vector(mode = "list", length = length(object))
    i.cc <- 0
    for (e in object) {
         i.cc <- i.cc + 1
         cc[[i.cc]] <- coef(e, ...)
        }
    cc
}


###################
# effects methods #
###################
effects.gvf.fit <- function(object, ...)
{
stats_effects.lm(object, ...)
}

effects.gvf.fits <- function(object, ...)
{
    ee <- vector(mode = "list", length = length(object))
    i.ee <- 0
    for (e in object) {
         i.ee <- i.ee + 1
         ee[[i.ee]] <- effects(e, ...)
        }
    ee
}


#####################
# residuals methods #
#####################
residuals.gvf.fit <- function(object, ...)
{
residuals.lm(object, ...)
}

residuals.gvf.fits <- function(object, ...)
{
    rr <- vector(mode = "list", length = length(object))
    i.rr <- 0
    for (e in object) {
         i.rr <- i.rr + 1
         rr[[i.rr]] <- residuals(e, ...)
        }
    rr
}


##################
# fitted methods #
##################
fitted.gvf.fit <- function(object, ...)
{
stats_fitted.default(object, ...)
}

fitted.gvf.fits <- function(object, ...)
{
    ff <- vector(mode = "list", length = length(object))
    i.ff <- 0
    for (e in object) {
         i.ff <- i.ff + 1
         ff[[i.ff]] <- fitted(e, ...)
        }
    ff
}


################
# vcov methods #
################
vcov.gvf.fit <- function(object, ...)
{
stats_vcov.lm(object, ...)
}

vcov.gvf.fits <- function(object, ...)
{
    vv <- vector(mode = "list", length = length(object))
    i.vv <- 0
    for (e in object) {
         i.vv <- i.vv + 1
         vv[[i.vv]] <- vcov(e, ...)
        }
    vv
}


getBest <- function(object,
                    criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
###############################################
# Generic function. Get best GVF fitted model #
# on the basis of a given criterion.          #
###############################################
{
UseMethod("getBest")
}

getBest.default <- function(object,
                            criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
####################################
# Default method: raises an error. #
####################################
{
stop("Input object is not a fitted GVF model!")
}

getBest.gvf.fit <- function(object,
                            criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
#####################################
# gvf.fit method: no actual choice, #
# returns object.                   #
#####################################
{
criterion <- match.arg(criterion)
object
}

getBest.gvf.fits <- function(object,
                             criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
################################################
# gvf.fits method: returns the fitted GVF with #
# highest score in the given criterion.        #
################################################
{
criterion <- match.arg(criterion)
score <- switch(criterion, "R2"     = getR2(object, ...),
                           "adj.R2" = getR2(object, adjusted = TRUE, ...),
                           "AIC"    = -AIC(object, ...),
                           "BIC"    = -BIC(object, ...)
                )

object[[which.max(score)]]
}


############################
# Possible methods TO ADD! #
############################

# influence.measures(model)

# dffits(model, infl = , res = )

# dfbeta(model, ...)

# dfbetas(model, ...)

# covratio(model, infl = lm.influence(model, do.coef = FALSE), res = weighted.residuals(model))

# cooks.distance(model, ...)

# hatvalues(model, ...)
