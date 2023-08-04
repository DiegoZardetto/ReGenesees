Zapsmall <- function(x, digits = getOption("digits"), ...){
UseMethod("Zapsmall", x)
}

Zapsmall.default <- function(x, digits = getOption("digits"), ...){
######################################################
# Default method: zapsmall function in package base. #
######################################################
base::zapsmall(x, digits = digits)
}

Zapsmall.data.frame <- function (x, digits = getOption("digits"), except = NULL, ...) {
    if (!inherits(x, "data.frame")) 
        stop("object 'x' is not of class data.frame")
    # Prevent havoc caused by tibbles:
    if (inherits(x, c("tbl_df", "tbl")))
        x <- as.data.frame(x)
    ix <- vapply(x, is.numeric, NA)
    if (!is.null(except))
        x[except] <- FALSE
    x[ix] <- lapply(x[ix], Zapsmall, ...)
    x
}
