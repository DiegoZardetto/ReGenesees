svystatB <- function(design, model, by = NULL,
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
# Get response variable name
respName <- get.respName(model)

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
    model.vars <- all.vars(model)
    by.vars <- all.vars(by)
    if (any(by.vars %in% model.vars))
         stop("Variables referenced by argument 'by' must not appear in the input 'model' formula!")
    stat <- svyby(formula = model, by = by, design = design, FUN = svylinB,
                  deff = deff, keep.var = TRUE, keep.names = TRUE,
                  verbose = FALSE, vartype = vartype, ci.lev = conf.lev,
                  drop.empty.groups = TRUE, covmat = FALSE, na.rm = na.rm,
                  return.replicates = FALSE)

    # Handle unfortunate list output cases caused by different aliasing in 'by'
    # subpops
    if (inherits(stat, "svyby.list")) {
         class(stat) <- c("svystatB.by.list", class(stat))
        } else {
         # Better column names -> embed response variable name separated by "_"
         # See below for function better.namesB()
         stat.names <- names(stat)
         regcoef.vars <- attr(stat,"svyby")$variables
         by.vars <- stat.names[attr(stat,"svyby")$margins]
         names(stat) <- better.namesB(regcoef.vars, by.vars, stat.names, respName)
         class(stat) <- c("svystatB.by", class(stat))
        }
    attr(stat,"model") <- model
    attr(stat,"respName") <- respName
    attr(stat,"design") <- design.expr
    return(stat)
    }
else {
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
    # resp.name <- gsub(" ", "", deparse(model[[2]]))
    est.name <- paste("RegCoef", respName, sep=".")
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
    attr(stat.df,"respName") <- respName
    attr(stat.df,"design") <- design.expr
    stat.df
    }
}

get.respName <- function(formula){
###################################################
# Get the name of the response term of a formula. #
# NOTE: The function checks that the formula:     #
#       1) do have a response term                #
#       2) the response term is a variable with   #
#          a name (not, e.g., an expression)      #
###################################################
  if (length(formula) < 3)
     stop(paste("The specified model does not have a response term: ", form.to.char(formula), sep = ""))
     respName <- formula[[2]]

  if (!is.name(respName))
     stop("Response term is not a named variable!")

  respName <- as.character(respName)

  return(respName)
}


## Accessor functions for overall population and *ordinary* domain estimates ##

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


## summary method for overall population estimates ##

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


## summary method for *ordinary* domain estimates ##

summary.svystatB.by <- function(object, ...){
# Build a lm-like coefficients matrix...
RegCoefs <- cbind(coef(object))
colnames(RegCoefs) <- paste("RegCoef", attr(object, "respName"), sep = ".")
cmat <- cbind(RegCoefs, "SE" = unlist(SE(object)))
# Compute Z statistics for estimated regression coefficients
z_val <- coef(object)/unlist(SE(object))
# Compute p values for the test b=0
Pr.gt.z <- 2*(1 - pnorm(abs(z_val)))
# Add columns to cmat
cmat <- cbind(cmat, "z value" = z_val, "Pr(>|z|)" = Pr.gt.z)
printCoefmat(cmat, signif.stars = TRUE, signif.legend = TRUE,
                     P.values = TRUE, has.Pvalue = TRUE, ...)
cat("\nNOTE: the distribution of regression coefficients estimators has been assumed normal\n\n")
invisible(cmat)
}


## Non-ordinary domain estimates:
## Print method that drops attributes
print.svystatB.by.list <- function(x, ...){
out <- x
attributes(out) <- NULL
names(out) <- names(x)
print(out)
return(invisible(x))
}


## Non-ordinary domain estimates:
## No accessor functions for svyby.list objects (yet)

coef.svystatB.by.list <- function(object, ...){
cat("\nNo 'coef' method available for this class, sorry.\n\n")
invisible(NULL)
}

vcov.svystatB.by.list <- function(object, ...){
cat("\nNo 'vcov' method available for this class, sorry.\n\n")
invisible(NULL)
}

SE.svystatB.by.list <- function(object, ...){
cat("\nNo 'SE' method available for this class, sorry.\n\n")
invisible(NULL)
}

VAR.svystatB.by.list <- function(object, ...){
cat("\nNo 'VAR' method available for this class, sorry.\n\n")
invisible(NULL)
}

cv.svystatB.by.list <- function(object, ...){
cat("\nNo 'cv' method available for this class, sorry.\n\n")
invisible(NULL)
}

deff.svystatB.by.list <- function(object, ...){
cat("\nNo 'deff' method available for this class, sorry.\n\n")
invisible(NULL)
}

confint.svystatB.by.list <- function(object, ...){
cat("\nNo 'confint' method available for this class, sorry.\n\n")
invisible(NULL)
}


## Non-ordinary domain estimates:
## No summary method for svyby.list objects (yet)

summary.svystatB.by.list <- function(object, ...){
cat("\nNo 'summary' method available for this class, sorry.\n\n")
invisible(NULL)
}



better.namesB <- function(y.vars, by.vars, stat.names, prefix, sepchar = "_"){
###########################################################
# Build better column names for svyby output by embedding #
# a customizable prefix (e.g. the model response name).   #
###########################################################

# Identify stat.names referring to y.vars statistics
is.ystat <- !(stat.names %in% by.vars)
ystat <- stat.names[is.ystat]

# Order y.vars by decreasing length (to avoid double updates for partial matches)
y.vars <- y.vars[order(sapply(y.vars, nchar), decreasing = TRUE)]

tmp.names <- old.names <- ystat
# First process names of variance measures (i.e. containing string ".yvars[i]")
for (y in y.vars){
     dot.y <- paste(sepchar, y, sep="")
     tmp.names[tmp.names==old.names] <- sub(dot.y, paste(sepchar, prefix, dot.y, sep=""),
                                            tmp.names[tmp.names==old.names], fixed = TRUE)
    }
# Then process names of estimates for quantitative vars (i.e. exact matches to "yvars[i]")
for (y in y.vars){
     tmp.names[tmp.names==y] <- paste(prefix, y, sep=sepchar)
     }
# Lastly process names of estimates for factor vars (i.e. partial matches to "yvars[i]")
for (y in y.vars){
     tmp.names[tmp.names==old.names] <- sub(y, paste(prefix, y, sep=sepchar),
                                            tmp.names[tmp.names==old.names], fixed = TRUE)
    }

# Check if all names have been processed
if (any(tmp.names==old.names))
    warning("Check output column names: some could be wrong/non-standard!")
stat.names[is.ystat] <- tmp.names
stat.names
}
