svySigma2 <- function(design, y, by = NULL,
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
     sigma2 <- svystatL(design, expression( (ones / (ones - 1)) * 
                                            ( (this.y2 / ones) - (this.y / ones)^2 )
                                        ),
                               by = by, vartype = vartype, conf.int = conf.int, conf.lev = conf.lev,
                               deff = deff, na.rm = na.rm)
    } else {
     sigma2 <- svystatL(design, expression(   
                                            ( (this.y2 / ones) - (this.y / ones)^2 )
                                        ),
                               by = by, vartype = vartype, conf.int = conf.int, conf.lev = conf.lev,
                               deff = deff, na.rm = na.rm)
    }

is.by <- inherits(sigma2, "svyby")

if (!is.by) {
     rownames(sigma2) <- y.char
     names(sigma2)[names(sigma2) == "Complex"] <- "Sigma2"
    } else {
     complex.expr <- attr(sigma2, "svyby")$variables
     names(sigma2) <- gsub(paste(".", complex.expr, sep = ""), "", names(sigma2), fixed = TRUE)
     names(sigma2)[names(sigma2) == complex.expr] <- paste("Sigma2", y.char, sep = ".")
    }

# Better internal names (instead of expr)
sigma2 <- sub.expr(sigma2, "Sigma2", y.char)

attr(sigma2, "design") <- design.expr
class(sigma2) <- c(ifelse(is.by, "svySigma2.by", "svySigma2") , class(sigma2)[class(sigma2) != ifelse(is.by, "svystatL.by", "svystatL")])
sigma2
}


coef.svySigma2 <- function(object, ...){
    by.object <- attr(object, "origin")
    aa <- attr(by.object, "svyby")
    rval <- by.object[, max(aa$margins) + (1:aa$nstats)]
    if (!is.null(dim(rval))){
        rval<-as.vector(as.matrix(rval))
        }
    names(rval) <- gsub("statistics\\.", "", aa$variables)
    rval
}

vcov.svySigma2 <- function(object, ...){
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

SE.svySigma2 <- function(object, ...) SE(attr(object, "origin"))

VAR.svySigma2 <- function(object, ...) VAR(attr(object, "origin"))

cv.svySigma2 <- function(object, ...) cv(attr(object, "origin"))

deff.svySigma2 <- function(object, ...) deff(attr(object, "origin"))

confint.svySigma2 <- function(object, ...){
    ci <- confint(attr(object, "origin"), ...)
    rownames(ci) <- sub("1:", "", x = rownames(ci), fixed = TRUE)
    ci
}
