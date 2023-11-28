`write.svystat` <-
function(x, ...)
{
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",    "svystatB",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by", "svystatB.by",
                "svySigma",       "svySigma2",
                "svySigma.by", "svySigma2.by",
                "CoV",                 "Corr",       "svyby")
check.classes <- sapply(OK.classes, function(class) inherits(x, class))
if (!any(check.classes))
    stop("Object has not the right class (should be ",paste(OK.classes, collapse=" or "),")")
Call <- match.call(expand.dots = TRUE)
for (argname in c("col.names", "row.names", "quote", "qmethod")) {
     if (!is.null(Call[[argname]]))
        warning(gettextf("attempt to set '%s' ignored", argname), 
                domain = NA)
    }
rows <- as.character(rownames(x))
stat.w <- cbind(rows, x)
rownames(stat.w) <- NULL
if (inherits(x, "svyby")) {
    names(stat.w)[1] <- "Domain"
    }
else {
    names(stat.w)[1] <- "Variable"
    if (data.class(x)=="svystatQ") {
         names(stat.w)[1] <- "Probability"
         # Strip the 'p = ' prefix
         stat.w[, "Probability"] <- sub("p = ", "", stat.w[, "Probability"])
        }
    if (data.class(x) %in% c("svystatS", "svystatSR")) {
         names(stat.w)[1] <- "Class"
        }
    if (data.class(x)=="svystatB") names(stat.w)[1] <- "Predictor"
    if (data.class(x) %in% c("CoV", "Corr")) stat.w[[1]] <- NULL
    # Classes Corr.by and CoV.by do NOT inherit from svyby
    if (data.class(x) %in% c("CoV.by", "Corr.by")) names(stat.w)[1] <- "Domain"
    }
Call$x <- stat.w
Call$col.names <- TRUE
Call$row.names <- FALSE
Call$quote <- FALSE
Call$qmethod <- "double"
Call[[1L]] <- as.name("write.table")
eval.parent(Call)
}
