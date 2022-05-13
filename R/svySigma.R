svySigma <- function(design, y, by = NULL,
                     fin.pop = TRUE,
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
    stop("Variable of interest must be supplied as a formula")
# Vartype
if (missing(vartype)) vartype <- "se"
vartype <- match.arg(vartype, several.ok = TRUE)
vartype <- unique(vartype)
# Interest variable
y.char <- all.vars(y)
if (length(y.char) > 1) 
    stop("Can specify only one variable of interest")
# Check for missing values in interest variables
NA.estvars(design = design, estvars = y.char, na.rm = na.rm)
# Type check for interest variables
y.typetest <- is.numeric(design$variables[, y.char])
if (!y.typetest) 
    stop("Variables of interest must be numeric")

design$variables$ones <- 1
design$variables$this.y  <- (design$variables[[y.char]])
design$variables$this.y2 <- (design$variables$this.y)^2

if (fin.pop) {
     sigma <- svystatL(design, expression( sqrt( (ones / (ones - 1)) * 
                                                 ( (this.y2 / ones) - (this.y / ones)^2 )
                                                )
                                        ),
                               by = by, vartype = vartype, conf.int = conf.int, conf.lev = conf.lev,
                               deff = deff, na.rm = na.rm)
    } else {
     sigma <- svystatL(design, expression( sqrt(  
                                                 ( (this.y2 / ones) - (this.y / ones)^2 )
                                                )
                                        ),
                               by = by, vartype = vartype, conf.int = conf.int, conf.lev = conf.lev,
                               deff = deff, na.rm = na.rm)
    }

is.by <- inherits(sigma, "svyby")

if (!is.by) {
     rownames(sigma) <- y.char
     names(sigma)[names(sigma) == "Complex"] <- "Sigma"
    } else {
     complex.expr <- attr(sigma, "svyby")$variables
     names(sigma) <- gsub(paste(".", complex.expr, sep = ""), "", names(sigma), fixed = TRUE)
     names(sigma)[names(sigma) == complex.expr] <- paste("Sigma", y.char, sep = ".")
    }

# Better internal names (instead of expr)
sigma <- sub.expr(sigma, "Sigma", y.char)

attr(sigma, "design") <- design.expr
class(sigma) <- c(ifelse(is.by, "svySigma.by", "svySigma") , class(sigma)[class(sigma) != ifelse(is.by, "svystatL.by", "svystatL")])
sigma
}


sub.expr <- function(object, stat, y.char) {
###############################################################################
# Given an object returned by svystatL, make the names of rows / cols of its  #
# inner svyby objects more standard and concise (i.e. avoid the reporting the #
# input expression to the original call.                                      #
# NOTE: This is meant to yield better output of extractors like coef, SE, cv, #
# VAR, confint, ...                                                           #
###############################################################################
     if (inherits(object, "svyby")) {
         attr(object, "svyby")$variables <- paste(stat, y.char, sep = ".")
        } else {
         expr <- attr(attr(object, "origin"),"svyby")$variables
         attr(attr(object, "origin"),"svyby")$variables <- paste(stat, y.char, sep = ".")
         names(attr(object, "origin")) <- gsub(paste(".", expr, sep = ""), "", names(attr(object, "origin")), fixed = TRUE)
         names(attr(object, "origin"))[names(attr(object, "origin")) == expr] <- paste(stat, y.char, sep = ".")
        }
     object
    }


coef.svySigma <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svySigma <- function(object, ...){
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

SE.svySigma <- function(object, ...) SE(attr(object, "origin"))

VAR.svySigma <- function(object, ...) VAR(attr(object, "origin"))

cv.svySigma <- function(object, ...) cv(attr(object, "origin"))

deff.svySigma <- function(object, ...) deff(attr(object, "origin"))

confint.svySigma <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
