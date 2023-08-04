`des.merge` <-
function(design, data, key)
#########################################################
# Survey data contained in 'design' are merged with the #
# 'data' dataframe on a common key column 'key'.        #
#########################################################
{
    # Check on 'design' class
    if (!inherits(design, "analytic"))
         stop("Object 'design' must inherit from class analytic")

    # Check on 'data' class
    if (!inherits(data, "data.frame"))
         stop("Object 'data' must inherit from class data.frame")
    # Prevent havoc caused by tibbles:
    if (inherits(data, c("tbl_df", "tbl")))
        data <- as.data.frame(data)

    # Check on 'key' class and length
    if (!inherits(key,"formula"))
         stop("'key' variable must be passed as a formula")
    key.char <- all.vars(key)
    if (length(key.char) < 1) 
         stop("'key' formula must reference a survey data variable")
    if (length(key.char) > 1) 
         stop("'key' formula must reference only one variable")

    # Check 'key' existence, NAs and type ('design' side)
    fail.key.design <- try(na.Fail(design$variables, key.char), silent = TRUE)
    if (inherits(fail.key.design, "try-error")){
         cond <- attr(fail.key.design, "condition")$message
         stop("Error on the 'design' side. ", cond)
        }
    if (!is.numeric(design$variables[, key.char]) && !is.character(design$variables[, key.char])) 
        stop("Error on the 'design' side. 'key' variable must be numeric or character")

    # Check 'key' existence, NAs and type ('data' side)
    fail.key.data <- try(na.Fail(data, key.char), silent = TRUE)
    if (inherits(fail.key.data, "try-error")){
         cond <- attr(fail.key.data, "condition")$message
         stop("Error on the 'data' side. ", cond)
        }
    if (!is.numeric(data[, key.char]) && !is.character(data[, key.char])) 
        stop("Error on the 'data' side. 'key' variable must be numeric or character")

    # Check nrow mismatch
    if ( (nr.des <- nrow(design$variables)) != (nr.data <- nrow(data)) )
         stop("'design' and 'data' must have the same number of rows")

    # Check that 'key' is actually a key ('design' side)
    des.key <- design$variables[[key.char]]
    if (length(unique(des.key)) != nr.des)
         stop("Error on the 'design' side. 'key' variable is not an actual key")

    # Check that 'key' is actually a key ('data' side)
    data.key <- data[[key.char]]
    if (length(unique(data.key)) != nr.data)
         stop("Error on the 'data' side. 'key' variable is not an actual key")

    # Check that 'key' variables in 'design' and 'data' are
    # actually in *1:1 correspondence*
    if (length(unique(c(des.key, data.key))) != nr.des){
         stop("'key' values in 'design' and 'data' are not in a 1:1 correspondence")
        }

    # If common variables exist in 'design' and 'data' (besides the 'key'),
    # then retain their 'design' version only
    datanames <- names(data)
    newnames <- datanames[datanames != key.char]
    ## If only the 'key' in 'data', do nothing
         if (length(newnames) < 1) {
             return(design)
            }
    desnames <- names(design$variables)
    desvars <- desnames[desnames != key.char]
    which.common <- newnames %in% desvars
    common <- newnames[which.common]
    not.common <- newnames[!which.common]

    if (any(which.common)){
         warning("Common variables found in 'design' and 'data' (besides the 'key'): ",
                 paste(common, collapse = ", "), ".\nOnly their 'design' version will be retained")
         # If no non-common variables in 'data', do nothing else
         if (length(not.common) < 1) {
             return(design)
            }
        }

    # Perform the merge
    ## Store original row ordering of 'design'
    design$variables$old.order <- 1:nr.des
    ## Order 'design' on 'key'
    design$variables <- design$variables[order(des.key), ]
    ## Order 'data' on 'key'
    data <- data[order(data.key), ]
    ## Add new (i.e. non-common) columns in 'data' to 'design'
    for (var in not.common) {
         design$variables[, var] <- data[[var]]
        }
    ## Re-order 'design' rows in the original way
    design$variables <- design$variables[order(design$variables$old.order), ]
    ## Drop the original ordering column
    design$variables$old.order <- NULL

    # Done
    design
}
