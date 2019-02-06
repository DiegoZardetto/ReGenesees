svystatQ <- function(design, y, probs = c(0.25, 0.5, 0.75), by = NULL,
                     vartype = c("se", "cv", "cvpct", "var"), conf.lev = 0.95,
                     na.rm = FALSE, ties=c("discrete", "rounded")){

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Only most essential check on input type (more in inner functions...)
if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")
if (!inherits(y, "formula")) 
    stop("Variable of interest must be supplied as a formula")
y.char <- all.vars(y)
if (length(y.char) > 1) 
    stop("Can specify only one variable of interest")
# Check for missing values in interest variables
NA.estvars(design = design, estvars = y.char, na.rm = na.rm)
# Type check for interest variables
y.typetest <- is.numeric(design$variables[, y.char])
if (!y.typetest) 
    stop("Variables of interest must be numeric")

# round probs to permille
probs <- round(probs, 3)
probs <- sort(unique(probs))
if (any(probs < 0.001 | probs > 0.999)) 
    stop("'probs' values must fall inside the interval [0.001,0.999]")
if (missing(vartype)) vartype <- "se"
vartype <- match.arg(vartype, several.ok = TRUE)
vartype <- unique(vartype)
# Confidence Intervals are always needed for quantiles (otherwise no variance estiamation...)
vartype <- c(vartype, "ci")
l.conf.tag <- paste("CI.l(", round(100 * conf.lev, 1), "%)", sep = "")
u.conf.tag <- paste("CI.u(", round(100 * conf.lev, 1), "%)", sep = "")
if (!is.numeric(conf.lev))
    stop("conf.lev must be numeric")
if (conf.lev < 0 || conf.lev > 1) 
    stop("conf.lev must fall inside [0,1]")
ties <- match.arg(ties)

if (!is.null(by)) {
    if (!inherits(by, "formula")) 
        stop("If specified, 'by' must be supplied as a formula")
    stat <- svyby(formula = y, by = by, design = design, FUN = svyquantile, quantiles = probs,
                  ties = ties, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE)
    # Better column names
    # Remark: old.names trick is needed to avoid multiple updates for names
    # were multiple matches are found (see svystatTM)
    tmp.names <- old.names <- names(stat)
    # 1) More than a single prob
    if (length(probs) > 1) {
        for( p in  rev(attr(stat, "svyby")$variables) ){
            tmp.names[tmp.names==old.names] <- sub(p, paste(y.char, ".Q[", p, "]", sep=""),
                                                    tmp.names[tmp.names==old.names])
            }
        }
    # 1) A single prob
    else{
        tmp.names[tmp.names==old.names] <- sub(y.char, paste(y.char, ".Q[", probs, "]", sep=""),
                                                tmp.names[tmp.names==old.names])
        # must treat SE and VAR separately, since in svyby output
        # they happen to be not pasted with the interest variable
        tmp.names[(tmp.names == old.names) & (tmp.names != y.char) & (tmp.names == "SE")] <-
        paste("SE.", y.char, ".Q[", probs, "]", sep="")
        tmp.names[(tmp.names == old.names) & (tmp.names != y.char) & (tmp.names == "VAR")] <-
        paste("VAR.", y.char, ".Q[", probs, "]", sep="")
        }
    names(stat) <- tmp.names
    class(stat) <- c("svystatQ.by", class(stat))
    attr(stat,"design") <- design.expr
    return(stat)
    }
else {
    stat <- svyby(formula = y, by = rep(1, nrow(design)), design = design, FUN = svyquantile, quantiles = probs,
                  ties = ties, keep.var = TRUE, keep.names = TRUE, verbose = FALSE,
                  vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE)

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
    # Rownames are the probability breaks
    rownames(stat.mat) <- paste("p = ", format(round(probs, 3), digits = 3, nsmall=3), sep="")
    # Colnames are Estimator type + Variability measures
    o.vartypes <- c("se","ci","ci","cv","cvpct","var")
    varia.name  <- o.vartypes[o.vartypes %in% attr(stat, "svyby")$vartype]
    colnames(stat.mat) <- c( paste(y.char, "Q[p]", sep="."),
                             varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatQ", class(stat.df))
    # Better column names
    tmp.names <- names(stat.df)
    tmp.names[match("se", tmp.names)]    <- "SE"
    tmp.names[match("cv", tmp.names)]    <- "CV"
    tmp.names[match("cvpct", tmp.names)] <- "CV%"
    tmp.names[match("var", tmp.names)]   <- "VAR"
    tmp.names[match("ci", tmp.names)] <- l.conf.tag
    tmp.names[match("ci", tmp.names)] <- u.conf.tag
    names(stat.df) <- tmp.names
    attr(stat.df,"design") <- design.expr
    stat.df
    }
}


coef.svystatQ <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatQ <- function(object, ...){
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

SE.svystatQ <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatQ <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatQ <- function(object, ...) cv(attr(object, "origin"))

# No deff for quantiles estimators, currently...
# deff.svystatQ <- function(object, ...) deff(attr(object, "origin"))

confint.svystatQ <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
