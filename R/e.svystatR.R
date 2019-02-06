svystatR <- function(design, num, den, by = NULL, cross = FALSE,
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
if (!inherits(num, "formula")) 
    stop("Numerator variables must be supplied as a formula")
num.char <- all.vars(num)
# Check for missing values in interest variables
NA.estvars(design = design, estvars = num.char, na.rm = na.rm)
# Type check for interest variables
num.typetest <- sapply(num.char, function(y) is.numeric(design$variables[, y]))
if (!all(num.typetest)) 
    stop("Numerator variables must be numeric")

if (!inherits(den, "formula")) 
    stop("Denominator variables must be supplied as a formula")
den.char <- all.vars(den)
# Check for missing values in interest variables
NA.estvars(design = design, estvars = den.char, na.rm = na.rm)
# Type check for interest variables
den.typetest <- sapply(den.char, function(y) is.numeric(design$variables[, y]))
if (!all(den.typetest)) 
    stop("Denominator variables must be numeric")

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
if (!is.null(by)) {
    if (!inherits(by, "formula")) 
        stop("If specified, 'by' must be supplied as a formula")
    stat <- svyby(formula = num, by = by, design = design, FUN = svyratio, denominator = den,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = cross)
    class(stat) <- c("svystatR.by", class(stat))
    attr(stat,"design") <- design.expr
    return(stat)
    }
else {
    stat <- svyby(formula = num, by = rep(1, nrow(design)), design = design, FUN = svyratio, denominator = den,
                  deff = deff, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE,
                  cross = cross)

  # Now change stat's format (recall stat is a dataframe)  
    # Drop fake 'by' column
    stat.f <- stat[, -which(names(stat)=="by"), drop=FALSE]
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
    colnames(stat.mat) <- c("Ratio", varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatR", class(stat.df))
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


coef.svystatR <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatR <- function(object, ...){
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

SE.svystatR <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatR <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatR <- function(object, ...) cv(attr(object, "origin"))

deff.svystatR <- function(object, ...) deff(attr(object, "origin"))

confint.svystatR <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
