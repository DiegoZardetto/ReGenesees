# BarPlot with confidence intervals visualized via error bars
BarPlotCI <- function(stat, level = 0.95, arr.len = 0.05, arr.col = "black", arr.lwd = 1,
                      pch = "", xlab = "Label", ylab = "Estimate", names.arg = NULL, ...){
#
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",    "svystatB",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by", "svystatB.by", 
                "svyby")
check.classes <- sapply(OK.classes, function(class) inherits(stat, class))
if (!any(check.classes))
    stop("Object has not the right class (should be ",paste(OK.classes, collapse=" or "),")")
#
e <- coef(stat)
ci <- as.numeric(as.matrix(confint(stat, level = level)))
ci.l <- ci[1:(length(ci)/2)]
ci.u <- ci[(length(ci)/2 + 1):length(ci)]
yylim <- range(0, ci.l, ci.u)
if (is.null(names.arg)) names.arg <- names(e)
x <- barplot(e, ylim = yylim, xlab = xlab, ylab = ylab, names.arg = names.arg, ...)
points(x, e, pch = pch, col = arr.col)
arrows(x, ci.l, x, ci.u, length = arr.len, angle = 90, code = 3, col = arr.col, lwd = arr.lwd)
}

# Plot with confidence intervals visualized via error bars
PlotCI <- function(stat, level = 0.95, arr.len = 0.05, arr.col = "black", arr.lwd = 1,
                   xlab = "Label", ylab = "Estimate", labels = NULL, ...){
#
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",    "svystatB",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by", "svystatB.by", 
                "svyby")
check.classes <- sapply(OK.classes, function(class) inherits(stat, class))
if (!any(check.classes))
    stop("Object has not the right class (should be ",paste(OK.classes, collapse=" or "),")")
#
e <- coef(stat)
ci <- as.numeric(as.matrix(confint(stat, level = level)))
ci.l <- ci[1:(length(ci)/2)]
ci.u <- ci[(length(ci)/2 + 1):length(ci)]
plot(e, ylim = range(ci.l, ci.u), xaxt = "n", col = arr.col, xlab = xlab, ylab = ylab, ...)
x <- 1:length(e)
if (is.null(labels)) labels <- names(e)
suppressWarnings(axis(1, at = x, labels = labels, ...))
arrows(x, ci.l, x, ci.u, length = arr.len, angle = 90, code = 3, col = arr.col, lwd = arr.lwd)
}
