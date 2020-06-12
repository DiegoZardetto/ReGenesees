svystatSR <- function(design, y, classes, by = NULL,
                      vartype = c("se", "cv", "cvpct", "var"),
                      conf.int = FALSE, conf.lev = 0.95, deff = FALSE,
                      na.rm = FALSE){

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Only most essential check on input type (more in inner functions...)
if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")
if (!inherits(y, "formula")) 
    stop("Interest variable must be supplied as a formula")
y.char <- all.vars(y)
if (length(y.char) > 1)
    stop("Can specify only one variable of interest")

# Check for missing values in interest variable
NA.estvars(design = design, estvars = y.char, na.rm = na.rm)
# Type check for interest variable
if (!is.numeric(design$variables[, y.char]))
    stop("Interest variable must be numeric")

# Checks on new argument classes
if (!inherits(classes, "formula")) 
    stop("Argument 'classes' must be supplied as a formula")
classes.char <- all.vars(classes)
if (length(classes.char) > 1)
    stop("Argument 'classes' can reference only one variable")
# Check for missing values in classes variables
NA.estvars(design = design, estvars = classes.char, na.rm = na.rm)
# Type check for interest variables
classes.typetest <- sapply(classes.char, function(y) is.factor(design$variables[, y]))
if (!all(classes.typetest)) 
    stop("The variable referenced by 'classes' must be factor")

# Checks on by variables
if (!is.null(by)) {
     if (!inherits(by, "formula")) 
         stop("If specified, 'by' must be supplied as a formula")
# No common variables in 'classes' and 'by'
# EDIT: could make sense, actually. Thus, statements below have been commented
#     by.vars <- all.vars(by)
#     if (any(by.vars %in% classes.char))
#         stop("No common variables allowed in formulae 'classes' and 'by'")
    }

if (missing(vartype)) vartype <- "se"
vartype <- match.arg(vartype, several.ok = TRUE)
vartype <- unique(vartype)
if (conf.int) {
    vartype <- c(vartype, "ci")
    l.conf.tag <- paste("CI.l(", round(100 * conf.lev, 1), "%)", sep = "")
    u.conf.tag <- paste("CI.u(", round(100 * conf.lev, 1), "%)", sep = "")
    if (!is.numeric(conf.lev))
        stop("conf.lev must be numeric")
    if (conf.lev < 0 || conf.lev > 1) 
        stop("conf.lev must fall inside [0,1]")
    }

# Build formula to generate share variables
f.y.classes <- as.formula(paste("~", paste(c(y.char, classes.char), collapse = ":"), " - 1", sep=""), env = .GlobalEnv)

# Build shares model matrix
smf <- model.frame(f.y.classes, data = design$variables, na.action = na.pass)
smm <- model.matrix(f.y.classes, data = smf)

# If by == NULL, change smm colnames dropping 'y:'
# if by != NULL, map ":" to deliberately weird string "THIS_WAS_COLON"
if (is.null(by)) {
     colnames(smm) <- gsub(paste(y.char, ":", sep = ""), "", colnames(smm))
    }
else {
     colnames(smm) <- gsub(":", "THIS_WAS_COLON", colnames(smm))
    }

# Compute possible ORDERED pairs of shares
snames <- colnames(smm)
scomb <- combn(snames, m = 2)
snum <- c(scomb[1, ], scomb[2, ])
sden <- c(scomb[2, ], scomb[1, ])
s.n <- length(snum)

# EXTEND smm (because svyratio does NOT ALLOW DUPLICATED variables in numerator
# and denominator formulae)
# NOTE: here chack.names = TRUE ensures that variable names are DE-DUPLICATED
smm <- data.frame(smm[, snum, drop = FALSE], smm[, sden, drop = FALSE], check.names = TRUE)

# Prepare a map linking names changed for DE-DUPLICATION to old names (for later
# post processing) 
old.names <- c(snum, sden)
new.names <- make.names(old.names, unique = TRUE)
new.changed <- new.names[new.names != old.names]
old.changed <- old.names[new.names != old.names]

# Build numerator and denominator formulae to compute share ratios
srnames <- colnames(smm)
srnum <- as.formula(paste("~`", paste(srnames[1:s.n],             collapse="`+`"), "`", sep = ""), env = .GlobalEnv)
srden <- as.formula(paste("~`", paste(srnames[(s.n + 1):(2*s.n)], collapse="`+`"), "`", sep = ""), env = .GlobalEnv)

# Fuse share ratios model matrix with survey data
design$variables <- data.frame(design$variables, smm, check.names = FALSE)

# Compute shares (which are "special" ratios)

if (!is.null(by)) {
    stat <- svyby(formula = srnum, by = by, design = design, FUN = svyratio, denominator = srden,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = FALSE)

      ### ONLY for shares, post-process output names
        # 1) RE-MAP names changed for DE-DUPLICATION to old names
        # 2) RESURGE ":"
        # 3) change variables names
      # 1)
      for (new in new.changed) {
             attr(stat, "svyby")$variables <- sub(new, old.changed[new.changed == new], attr(stat, "svyby")$variables)
            }
      # 2)
      attr(stat, "svyby")$variables <- gsub("THIS_WAS_COLON", ":", attr(stat, "svyby")$variables)
      # 3)
      attr(stat, "svyby")$variables <- gsub(paste("/", y.char, ":", sep = ""), "/", attr(stat, "svyby")$variables)

      # 1)
      for (new in new.changed) {
             names(stat) <- sub(new, old.changed[new.changed == new], names(stat))
            }
      # 2)
      names(stat) <- gsub("THIS_WAS_COLON", ":", names(stat))
      # 3)
      names(stat) <- gsub(paste("/", y.char, ":", sep = ""), "/", names(stat))
      ### DONE

    class(stat) <- c("svystatSR.by", class(stat))
    attr(stat,"design") <- design.expr
    return(stat)
    }
else {
    stat <- svyby(formula = srnum, by = rep(1, nrow(design)), design = design, FUN = svyratio, denominator = srden,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = FALSE)

  # Now change stat's format (recall stat is a dataframe)  
    # Drop fake 'by' column
    stat.f <- stat[, -which(names(stat)=="by"), drop=FALSE]

      # NOTE: If by == NULL, already changed variables names
      ### ONLY for shares, post-process output names
        # 1) RE-MAP names changed for DE-DUPLICATION to old names
      # 1)
      for (new in new.changed) {
             attr(stat, "svyby")$variables <- sub(new, old.changed[new.changed == new], attr(stat, "svyby")$variables)
            }
      # 2)
      for (new in new.changed) {
             names(stat) <- sub(new, old.changed[new.changed == new], names(stat))
            }
      ### DONE

    # Count how many stats (numbers) inside
    nstats <- ncol(stat.f)
    # Count how many variability measures inside (recall "ci" counts twice but already in attr(stat, "svyby")$vars)
    nvariances <- attr(stat, "svyby")$vars + attr(stat, "svyby")$deffs
    # Count how many variables inside
    nvariables <- attr(stat, "svyby")$nstats
    # Consistency check
      # DEBUG 12/06/2020: BUG: if clauses have to use !isTRUE(all.equal()). Fixed
    if (!isTRUE(all.equal(nstats, nvariables*(1 + nvariances)))) stop("Format modification failed")
    # Store properly the relevant stats into a matrix
    stat.mat <- matrix(nrow = nvariables, ncol = (1 + nvariances))
    stat.mat[, ] <- as.matrix(stat.f[, ])
    # Rownames are the variables names
    rownames(stat.mat) <- attr(stat, "svyby")$variables
    # Colnames are Estimator type + Variability measures
    o.vartypes <- c("se","ci","ci","cv","cvpct","var")
    varia.name  <- o.vartypes[o.vartypes %in% attr(stat, "svyby")$vartype]
    colnames(stat.mat) <- c(paste(y.char, "ShareRatio", sep = "."), varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatSR", class(stat.df))
    # Better column names
    tmp.names <- names(stat.df)
    tmp.names[match("se", tmp.names)]    <- "SE"
    tmp.names[match("cv", tmp.names)]    <- "CV"
    tmp.names[match("cvpct", tmp.names)] <- "CV%"
    tmp.names[match("var", tmp.names)]   <- "VAR"
    if (conf.int) {
        tmp.names[match("ci", tmp.names)] <- l.conf.tag
        tmp.names[match("ci", tmp.names)] <- u.conf.tag
    }
    names(stat.df) <- tmp.names
    attr(stat.df,"design") <- design.expr
    stat.df
    }
}


coef.svystatSR <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatSR <- function(object, ...){
    by.object <- attr(object, "origin")
    rval <- attr(by.object,"var")
    if(is.null(rval)){
       # warning("Only diagonal elements of vcov() available")
       se <- SE(by.object)
    if (is.data.frame(se)) se <- as.vector(as.matrix(se))
    if (length(se) > 1)
        rval<-diag(se^2)
    else
        rval<-as.matrix(se^2)
  }
  nms <- names(coef(object))
  dimnames(rval)<-list(nms,nms)
  rval
}

SE.svystatSR <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatSR <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatSR <- function(object, ...) cv(attr(object, "origin"))

deff.svystatSR <- function(object, ...) deff(attr(object, "origin"))

confint.svystatSR <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
