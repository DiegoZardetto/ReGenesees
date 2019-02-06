`ez.design` <-
function (ids, variables = NULL, data = NULL, weights = NULL)
#######################################################################
#  Versione modificata della funzione svydesign del package survey    #
#  NOTA: Per motivi di efficienza sono state eliminate TUTTE le       #
#        funzionalita' originali NON NECESSARIE per e.calibrate       #
#######################################################################
{
    interaction <- function(..., drop = TRUE) {
        args <- list(...)
        narg <- length(args)
        if (narg == 1 && is.list(args[[1]])) {
            args <- args[[1]]
            narg <- length(args)
        }
        ans <- do.call("paste", c(lapply(args, as.character), 
            sep = "."))
        ans <- factor(ans)
        return(ans)
    }
    ids <- eval.parent(substitute(model.frame(ids, data = data)))
    weights <- eval.parent(substitute(model.frame(weights, data = data)))
    dir.weights <- as.numeric(as.matrix(weights))
    if (is.character(variables)) {
        variables <- data[, variables]
    }
    else {
        variables <- data
    }
    if (!is.factor(ids[, 1]))
        ids[, 1] <- as.factor(ids[, 1])
    if (NCOL(ids) > 1) {
        N <- ncol(ids)
        for (i in 2:N) {
            ids[, i] <- do.call("interaction", ids[, 1:i, drop = FALSE])
        }
    }
    design <- list(cluster = ids)
    design$dir.weights <- dir.weights
    design$variables <- variables
    design
}
