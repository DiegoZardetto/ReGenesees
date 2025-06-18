###############################################################
# BarPlot with confidence intervals visualized via error bars #
###############################################################
BarPlotCI <- function(stat, level = 0.95, eb.len = 0.05, eb.col = "black", eb.lwd = 1,
                      xlab = "Domain", ylab = "Estimate", names.arg = NULL, pch = "", ...){

# Allowed summary statistics
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",    "svystatB",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by", "svystatB.by",
                "svySigma",       "svySigma2",
                "svySigma.by", "svySigma2.by",
                "svyby")
check.classes <- sapply(OK.classes, function(class) inherits(stat, class))
if (!any(check.classes))
    stop("Object has not the right class (should be ",paste(OK.classes, collapse=" or "),")")

# No support for horizontal bars
Call <- match.call(expand.dots = TRUE)
if (!is.null(Call[["horiz"]])) {
     stop("Horizontal bars are not supported!")
    }

e <- coef(stat)
ci <- as.numeric(as.matrix(confint(stat, level = level)))
ci.l <- ci[1:(length(ci)/2)]
ci.u <- ci[(length(ci)/2 + 1):length(ci)]
yylim <- range(0, ci.l, ci.u)
if (is.null(names.arg)) names.arg <- names(e)

x <- barplot(e, ylim = yylim, xlab = xlab, ylab = ylab, names.arg = names.arg, horiz = FALSE, ...)
points(x, e, pch = pch, col = eb.col)
arrows(x, ci.l, x, ci.u, length = eb.len, angle = 90, code = 3, col = eb.col, lwd = eb.lwd)
}


###################################################################
# ScatterPlot with confidence intervals visualized via error bars #
###################################################################
PlotCI <- function(stat, level = 0.95, eb.len = 0.05, eb.col = "black", eb.lwd = 1,
                   xlab = "Domain", ylab = "Estimate", labels = NULL, ...){

# Allowed summary statistics
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",    "svystatB",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by", "svystatB.by",
                "svySigma",       "svySigma2",
                "svySigma.by", "svySigma2.by",
                "svyby")
check.classes <- sapply(OK.classes, function(class) inherits(stat, class))
if (!any(check.classes))
    stop("Object has not the right class (should be ",paste(OK.classes, collapse=" or "),")")

e <- coef(stat)
ci <- as.numeric(as.matrix(confint(stat, level = level)))
ci.l <- ci[1:(length(ci)/2)]
ci.u <- ci[(length(ci)/2 + 1):length(ci)]

plot(e, ylim = range(ci.l, ci.u), xaxt = "n", col = eb.col, xlab = xlab, ylab = ylab, ...)
x <- 1:length(e)
if (is.null(labels)) labels <- names(e)
suppressWarnings(axis(1, at = x, labels = labels))
arrows(x, ci.l, x, ci.u, length = eb.len, angle = 90, code = 3, col = eb.col, lwd = eb.lwd)
}
