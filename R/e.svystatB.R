svystatB <- function(design, model,
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
if (!inherits(model, "formula"))
    stop("Linear regression model must be specified as a formula")
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

# Compute estimates and errors of multiple regression coefficients
    stat <- svyby(formula = model, by = rep(1, nrow(design)), design = design, FUN = svylinB,
                  deff = deff, keep.var = TRUE, keep.names = TRUE,
                  verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
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
    # Rownames are the variables names
    rownames(stat.mat) <- attr(stat, "svyby")$variables
    # Colnames are Estimator type + Variability measures
    o.vartypes <- c("se","ci","ci","cv","cvpct","var")
    varia.name  <- o.vartypes[o.vartypes %in% attr(stat, "svyby")$vartype]
    resp.name <- gsub(" ", "", deparse(model[[2]]))
    est.name <- paste("RegCoef", resp.name, sep=".")
    colnames(stat.mat) <- c(est.name, varia.name, if (attr(stat, "svyby")$deffs) "DEff")
    # Cast into a data.frame
    stat.df <- as.data.frame(stat.mat)
    attr(stat.df, "origin") <- stat
    class(stat.df) <- c("svystatB", class(stat.df))
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
    attr(stat.df,"model") <- model
    attr(stat.df,"design") <- design.expr
    stat.df
}


coef.svystatB <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svystatB <- function(object, ...){
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

SE.svystatB <- function(object, ...) SE(attr(object, "origin"))

VAR.svystatB <- function(object, ...) VAR(attr(object, "origin"))

cv.svystatB <- function(object, ...) cv(attr(object, "origin"))

deff.svystatB <- function(object, ...) deff(attr(object, "origin"))

confint.svystatB <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}

summary.svystatB <- function(object, ...){
# Compute Z statistics for estimated regression coefficients
object[["z value"]] <- coef(object)/unlist(SE(object))
# Compute p values for the test b=0
object[["Pr(>|z|)"]] <- 2*(1 - pnorm(abs(object[["z value"]])))
printCoefmat(object, signif.stars = TRUE, signif.legend = TRUE,
                     P.values = TRUE, has.Pvalue = TRUE, ...)
cat("\nNOTE: the distribution of regression coefficients estimators has been assumed normal\n\n")
invisible(object)
}
