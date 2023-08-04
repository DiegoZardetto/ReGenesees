svystatTM <- function(design, y, by = NULL, estimator = c("Total", "Mean"),
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
    stop("Interest variables must be supplied as a formula")
y.vars <- all.vars(y)
# Check for missing values in interest variables
NA.estvars(design = design, estvars = y.vars, na.rm = na.rm)
#if (isTRUE(na.rm)) {
#    if (length(y.vars) > 1) stop("na.rm = TRUE allowed only when a single interest variable is supplied!")
#    }
# Type check for interest variables
typetest <- sapply(y.vars, function(y) is.numeric(design$variables[, y]) |
                                       is.factor(design$variables[, y])   )
if (!all(typetest))
    stop("Interest variables must be numeric or factor")
estimator <- match.arg(estimator)
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
    stat <- switch(estimator, 
                  Total = svyby(formula = y, by = by, design = design, FUN = svytotal,
                                deff = deff, keep.var = TRUE, keep.names = TRUE,
                                verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
                                drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                                return.replicates = FALSE),
                  Mean  = svyby(formula = y, by = by, design = design, FUN = svymean,
                                deff = deff, keep.var = TRUE, keep.names = TRUE,
                                verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
                                drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                                return.replicates = FALSE)
                  )
    # Better column names -> embed estimator type (Total or Mean)
    stat.names <- names(stat)
    by.vars <- stat.names[attr(stat,"svyby")$margins]
    names(stat) <- better.names(y.vars, by.vars, stat.names, estimator)
    class(stat) <- c("svystatTM.by", class(stat))
    attr(stat,"design") <- design.expr
    attr(stat,"y.vars") <- y.vars
    attr(stat,"by.vars") <- by.vars
    return(stat)
    }
else {
    stat <- switch(estimator, 
                   Total = svyby(formula = y, by = rep(1, nrow(design)), design = design, FUN = svytotal,
                                 deff = deff, keep.var = TRUE, keep.names = TRUE,
                                 verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
                                 drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                                 return.replicates = FALSE),
                   Mean  = svyby(formula = y, by = rep(1, nrow(design)), design = design, FUN = svymean,
                                 deff = deff, keep.var = TRUE, keep.names = TRUE,
                                 verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
                                 drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                                 return.replicates = FALSE)
                  )
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
    colnames(stat.mat) <- c(estimator, varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatTM", class(stat.df))
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
    attr(stat.df,"y.vars") <- y.vars
    stat.df
    }
}


better.names <- function(y.vars, by.vars, stat.names, estimator="Total"){
#################################################
# Build better column names for svyby output by #
# embedding estimator type ("Total" or "Mean"). #
#################################################

# Identify stat.names referring to y.vars statistics
is.ystat <- !(stat.names %in% by.vars)
ystat <- stat.names[is.ystat]

# Order y.vars by decreasing length (to avoid double updates for partial matches)
y.vars <- y.vars[order(sapply(y.vars, nchar), decreasing = TRUE)]

tmp.names <- old.names <- ystat
# First process names of variance measures (i.e. containing string ".yvars[i]")
for (y in y.vars){
     dot.y <- paste(".", y, sep="")
     tmp.names[tmp.names==old.names] <- sub(dot.y, paste(".", estimator, dot.y, sep=""),
                                            tmp.names[tmp.names==old.names], fixed = TRUE)
    }
# Then process names of estimates for quantitative vars (i.e. exact matches to "yvars[i]")
for (y in y.vars){
     tmp.names[tmp.names==y] <- paste(estimator, y, sep=".")
     }
# Lastly process names of estimates for factor vars (i.e. partial matches to "yvars[i]")
for (y in y.vars){
     tmp.names[tmp.names==old.names] <- sub(y, paste(estimator, y, sep="."),
                                            tmp.names[tmp.names==old.names], fixed = TRUE)
    }

# Check if all names have been processed
if (any(tmp.names==old.names))
    warning("Check output column names: some could be wrong/non-standard!")
stat.names[is.ystat] <- tmp.names
stat.names
}


coef.svystatTM <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatTM <- function(object, ...){
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

SE.svystatTM <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatTM <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatTM <- function(object, ...) cv(attr(object, "origin"))

deff.svystatTM <- function(object, ...) deff(attr(object, "origin"))

confint.svystatTM <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
