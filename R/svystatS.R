svystatS <- function(design, y, classes, by = NULL,
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
# Check for missing values in classes variables
NA.estvars(design = design, estvars = classes.char, na.rm = na.rm)
# Type check for interest variables
classes.typetest <- sapply(classes.char, function(y) is.factor(design$variables[, y]))
if (!all(classes.typetest)) 
    stop("Variables referenced by 'classes' must be factor")

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
if (is.null(by)) 
    colnames(smm) <- gsub(paste(y.char, ":", sep = ""), "", colnames(smm))

# Fuse shares model matrix with survey data
design$variables <- data.frame(design$variables, smm, check.names = FALSE)

# Build formula to compute shares
snum <- as.formula(paste("~`", paste(colnames(smm), collapse="`+`"), "`", sep = ""), env = .GlobalEnv)

# Compute shares (which are "special" ratios)
if (!is.null(by)) {
    stat <- svyby(formula = snum, by = by, design = design, FUN = svyratio, denominator = y,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = FALSE)

      # ONLY for shares, change variables names
      attr(stat, "svyby")$variables <- gsub(paste("/", y.char, sep = ""), "", attr(stat, "svyby")$variables)
      names(stat) <- gsub(paste("/", y.char, sep = ""), "", names(stat))

    class(stat) <- c("svystatS.by", class(stat))
    attr(stat,"design") <- design.expr
    return(stat)
    }
else {
    stat <- svyby(formula = snum, by = rep(1, nrow(design)), design = design, FUN = svyratio, denominator = y,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = FALSE)

  # Now change stat's format (recall stat is a dataframe)  
    # Drop fake 'by' column
    stat.f <- stat[, -which(names(stat)=="by"), drop=FALSE]

      # ONLY for shares, change variables names
      attr(stat, "svyby")$variables <- gsub(paste("/", y.char, sep = ""), "", attr(stat, "svyby")$variables)
      names(stat) <- gsub(paste("/", y.char, sep = ""), "", names(stat))

    # Count how many stats (numbers) inside
    nstats <- ncol(stat.f)
    # Count how many variability measures inside (recall "ci" counts twice but already in attr(stat, "svyby")$vars)
    nvariances <- attr(stat, "svyby")$vars + attr(stat, "svyby")$deffs
    # Count how many variables inside
    nvariables <- attr(stat, "svyby")$nstats
    # Consistency check
    if (!all.equal(nstats, nvariables*(1 + nvariances))) stop("Format modification failed")
    # Store properly the relevant stats into a matrix
    stat.mat <- matrix(nrow = nvariables, ncol = (1 + nvariances))
    stat.mat[, ] <- as.matrix(stat.f[, ])
    # Rownames are the variables names
    rownames(stat.mat) <- attr(stat, "svyby")$variables
    # Colnames are Estimator type + Variability measures
    o.vartypes <- c("se","ci","ci","cv","cvpct","var")
    varia.name  <- o.vartypes[o.vartypes %in% attr(stat, "svyby")$vartype]
    colnames(stat.mat) <- c(paste(y.char, "Share", sep = "."), varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatS", class(stat.df))
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


coef.svystatS <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatS <- function(object, ...){
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

SE.svystatS <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatS <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatS <- function(object, ...) cv(attr(object, "origin"))

deff.svystatS <- function(object, ...) deff(attr(object, "origin"))

confint.svystatS <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
