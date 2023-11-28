svystat <- function(design, kind = c("TM", "R", "S", "SR", "B", "Q", "L", "Sigma", "Sigma2"),
                    by = NULL, group = NULL, forGVF = TRUE, combo = -1, ...) {
################################################################################
# This function can compute all the summary statistics provided by ReGenesees, #
# and is principally meant to return a lot of them in just a single shot.      #
# If 'forGVF' is TRUE the output will be ready to feed ReGenesees GVF fitting  #
# infrastructure, otherwise the output will be a set of summary statistic      #
# objects.                                                                     #
# The 'group' formula (if any) specifies a way of partitioning the population  #
# into groups: the output will be reported separately for each group.          #
# Parameter 'combo' is only meaningful if 'by' is passed; if the 'by' formula  #
# involves n variables, specifying combo = m requests to compute outputs for   #
# all the domains determined by all the interactions of 'by' variables up to   #
# order m (-1 <= m <= n):                                                      #
# - m = -1 means "no combo", i.e. treat by formula as usual;                   #
# - m =  0 means "order zero" combination, i.e. just a single domain, which is #
#          the whole population;                                               #
# - m =  1 means "order zero" plus "order one" combinations, the latter being  #
#          all the marginal domains defined by 'by' variables;                 #
# - m =  n means combinations of any order, the maximum being the one with all #
#          'by' variables interacting simultaneously.                          #
#                                                                              #
# NOTE: Class "svystat", which already exists in ReGenesees (and survey), has  #
#       *NOTHING* to do with output objects of function svystat().             #
#       Object created by svystat have the following data.classes:             #
#       - gvf.input    if forGVF = TRUE  and group =  NULL                     #
#       - gvf.input.gr if forGVF = TRUE  and group != NULL                     #
#       - svystat.gr   if forGVF = FALSE and the output is a *LIST* of summary #
#                                            statistics                        #
#       - other        if forGVF = FALSE and the output is a summary statistic #
#                                            object of class 'other'           #
#                                                                              #
################################################################################
# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Test design class
if (!inherits(design, "analytic")) 
     stop("Object 'design' must inherit from class analytic")

kind <- match.arg(kind)

by.is.NULL <- is.null(by)
if (!by.is.NULL && !inherits(by, "formula")) 
     stop("If specified, 'by' must be supplied as a formula")

if (!is.null(group) && !inherits(group, "formula")) 
     stop("If specified, 'group' must be supplied as a formula")

# DEBUG 06/11/2023
# Code below commented because, from ReGenesees 2.1, we DO have a 'by' argument for
# svystatB!
# # DEBUG 14/11/2017
# # Handle special case: kind = 'B' with either of 'by' and 'group' not NULL.
# # This is because function svystatB does NOT have a 'by' argument.
# if ( (kind == "B") && ( !by.is.NULL || !is.null(group) ) ) {
#      stop("When kind == 'B' cannot specify 'by' nor 'group' (because of possible aliasing issues, see ?svystatB)")
#     }

by.var <- all.vars(by)
group.var <- all.vars(group)

if (any(group.var %in% by.var)) {
     stop("'Group' variables cannot appear in 'by' formula!")
    }

combo <- as.integer(round(combo))

if ((combo != -1) && !by.is.NULL) {
     by <- combf(by, m = combo)
    }
else {
     by <- list(by)
    }

stats <- lapply(by, FUN = function(bb) svystat1(design = design, kind = kind,
                                                by = bb, group = group, forGVF = forGVF,
                                                design.expr = design.expr, ...)
                )
if (length(stats) == 1L) stats <- stats[[1L]]

if ( isTRUE(forGVF) && (combo >= 1) ) {
     # post-production: re-structures the output list in order to get...
     if ( !by.is.NULL && !is.null(group) ) {
         # ...EITHER a LIST of gvf.input objects, one for each group
         # NOTE: this deals only with sub-cases: !by.is.NULL .and. !is.null(group)
         outnames <- lapply(stats, names)[[1L]]
         outlen <- length(outnames)
         stats.new <- vector(mode = "list", length = outlen)
         for (j in 1:outlen) {
              dd.g <- lapply(stats, function(el) el[[j]])
              stats.new[[j]] <- Reduce(f = rbind, x = dd.g)
            }
         names(stats.new) <- outnames
         stats <- stats.new
        }
     if ( !by.is.NULL && is.null(group) ) {
         # ...OR a SINGLE gvf.input object
         # NOTE: this deals only with sub-cases: !by.is.NULL .and. is.null(group)
         stats <- Reduce(f = rbind, x = stats)
        }
    }

attr(stats, "group.vars") <- if (!is.null(group)) group.var

if (data.class(stats) == "list") {
     if ( isTRUE(forGVF) ) {
         class(stats) <- c("gvf.input.gr", class(stats))
        }
     else {
         class(stats) <- c("svystat.gr", class(stats))
        }
    }

stats
}

is.wrapped <- function(svystat.gr) {
#################################################
# Is systat.gr a 'nested' list where statistics #
# are wrapped at an inner level?                #
#################################################
# Get top level object (i.e. first level slot)
top <- svystat.gr[[1]]

# Summary statistics classes 
OK.classes <- c("svystatTM",       "svystatR",    "svystatL",    "svystatQ",    "svystatS",    "svystatSR",
                "svystatTM.by", "svystatR.by", "svystatL.by", "svystatQ.by", "svystatS.by", "svystatSR.by",
                "svySigma",       "svySigma2",
                "svySigma.by", "svySigma2.by",
                "svystatB",     "svyby")

# Is top level object a summary statistic?
check.classes <- sapply(OK.classes, function(class) inherits(top, class))

# Thus, is the object wrapped?
return(!any(check.classes))
}

svystat1 <- function(design, kind = c("TM", "R", "S", "SR", "B", "Q", "L", "Sigma", "Sigma2"),
                     by = NULL, group = NULL, forGVF = TRUE, design.expr = NULL, ...) {
##############################################################################
# Workhorse function serving exported function svystat: does the brute force #
# job and get called for each combination when 'combo' is used in svystat.   #
##############################################################################
kind <- match.arg(kind)
if (kind %in% c("Sigma", "Sigma2")) {
     FUN <- paste("svy", kind, sep = "")
    } else {
     FUN <- paste("svystat", kind, sep = "")
    }
by.var <- all.vars(by)
group.var <- all.vars(group)
group.by <- paste(c(group.var, by.var), collapse = ":")
if (nchar(group.by) > 0) {
     group.by <- as.formula(paste("~", group.by, sep=""), env = .GlobalEnv)
    }
else {
     group.by <- NULL
    }

stats <- do.call(FUN, c(list(design = design, by = group.by), list(...)))

if (!is.null(group)) {
     stats <- split(stats, stats[group.var], drop = TRUE)
    }

if (isTRUE(forGVF)) {
     if (!is.null(group)) {
         stats <- lapply(stats, FUN = function(ee) {
                                         gi <- gvf.input(design = design, ee)
                                         if (!is.null(design.expr)) {
                                             attr(gi, "design") <- design.expr
                                            }
                                         gi
                                        }
                        )
        }
     else {
         stats <- gvf.input(design = design, stats)
         if (!is.null(design.expr)) attr(stats, "design") <- design.expr
        }
    }

if ( !isTRUE(forGVF) && !is.null(design.expr) ) {
     attr(stats, "design") <- design.expr
    }
stats
}


combf <- function(formula, m) {
############################################################################
# Given a formula referencing n variables, returns a list storing all the  #
# formulae containing interactions of those variables up to order m, with  #
# 0 <= m <= n.                                                             #
# NOTE: m = 0 yields ~1.                                                   #
############################################################################
vars <- all.vars(formula)
nvars <- length(vars)
m <- as.integer(round(m))

out0 <- list(NULL)
names(out0) <- "population"
if (m == 0) {
     return(out0)
    }

if ( !(m %in% 1:nvars) ) stop("For the current formula, parameter 'combo' must be in -1:", nvars)
out <- vector(mode = "list", length = m)

for (i in 1:m) {
     out[[i]] <- combn(vars, m = i, simplify = FALSE,
                       FUN = function(v) as.formula(paste("~", paste(v, collapse = ":", sep = ""), sep = ""),
                                                    env = .GlobalEnv)
                    )
    }
out <- unlist(out)
fnames <- sapply(out, function(f) gsub("~ ", "", form.to.char(f)))
names(out) <- fnames
out <- c(out0, out)
out
}


################################
# Methods for class svystat.gr #
################################
## coef
coef.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, coef, ...))
    }
else {
     lapply(object, coef, ...)
    }
}

## vcov
vcov.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, vcov, ...))
    }
else {
     lapply(object, vcov, ...)
    }
}

## SE
SE.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, SE, ...))
    }
else {
     lapply(object, SE, ...)
    }
}

## VAR
VAR.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, VAR, ...))
    }
else {
     lapply(object, VAR, ...)
    }
}

## cv
cv.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, cv, ...))
    }
else {
     lapply(object, cv, ...)
    }
}

## deff
deff.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, deff, ...))
    }
else {
     lapply(object, deff, ...)
    }
}

## confint
confint.svystat.gr <- function(object, ...)
{
if (is.wrapped(object)) {
     lapply(object, function(el) lapply(el, confint, ...))
    }
else {
     lapply(object, confint, ...)
    }
}


##################################
# Methods for class gvf.input.gr #
##################################
## Plot
plot.gvf.input.gr <- function(x, ...) 
{
    # How many input?
    ng.i <- length(x)
    groups <- names(x)

    if (ng.i > 1) {
         oask <- devAskNewPage(TRUE)
         on.exit(devAskNewPage(oask))
        }

    jg.i <- 0
    for (jg.i in 1:ng.i) {
         plot(x[[jg.i]], main = paste("Group: ", groups[[jg.i]], sep=""), ...)
        }
}

## Fit
fit.gvf.gvf.input.gr <- function(gvf.input, model = NULL, weights = NULL)
{
    # First verify if the function has been called inside another function:
    # this is needed to correctly manage metadata when e.g. the caller is a
    # GUI stratum
    directly <- !( length(sys.calls()) > 2L )
    gvf.input.expr <- if (directly) substitute(gvf.input)

    mm <- vector(mode = "list", length = length(gvf.input))

    for (i.mm in 1:length(gvf.input)) {
         mm[[i.mm]] <- fit.gvf(gvf.input[[i.mm]], model = model, weights = weights)

         if (is.gvf.fits <- inherits(mm[[i.mm]], "gvf.fits")) {
             for (j.i.mm in 1:length(mm[[i.mm]])) {
                  attr(mm[[i.mm]][[j.i.mm]], "gvf.input.expr") <- gvf.input.expr
                }
            }
         else {
             attr(mm[[i.mm]], "gvf.input.expr") <- gvf.input.expr
            }
        }

    names(mm) <- names(gvf.input)

    attr(mm, "group.vars") <- attr(gvf.input, "group.vars") 

    class(mm) <- if (is.gvf.fits) c("gvf.fits.gr", class(mm)) else c("gvf.fit.gr", class(mm))
    mm
}


##################################################
# Methods for classes gvf.fit.gr and gvf.fits.gr #
##################################################
# print
print.gvf.fit.gr <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    # Strip group.vars attribute, otherwise it would be printed...
    x.strip <- x
    attr(x.strip, "group.vars") <- NULL

    print(unclass(x.strip), digits = digits, ...)
    invisible(x)
}
print.gvf.fits.gr <- print.gvf.fit.gr

# summary
summary.gvf.fit.gr <- function(object, correlation = FALSE,
                               symbolic.cor = FALSE, ...)
{
    lapply(object, FUN = summary, correlation = correlation,
           symbolic.cor = symbolic.cor, ...)
}
summary.gvf.fits.gr <- summary.gvf.fit.gr

# plot
plot.gvf.fit.gr <- function(x, which.more = 1:3, id.n = 3,
                            labels.id = NULL, cex.id = 0.75,
                            label.pos = c(4, 2), cex.caption = 1, ...)
{
    # How many input?
    ng.i <- length(x)
    groups <- names(x)

    if (ng.i > 1) {
         oask <- devAskNewPage(TRUE)
         on.exit(devAskNewPage(oask))
        }

    jg.i <- 0
    for (jg.i in 1:ng.i) {
         plot(x[[jg.i]], which.more = which.more, id.n = id.n, labels.id = labels.id,
              cex.id = cex.id, label.pos = label.pos, cex.caption = cex.caption,
              Main = paste("Group: ", groups[[jg.i]], sep=""), ...)
        }
}
plot.gvf.fits.gr <- function(x, which.more = NULL, id.n = 3,
                            labels.id = NULL, cex.id = 0.75,
                            label.pos = c(4, 2), cex.caption = 1, ...)
{
    # How many input?
    ng.i <- length(x)
    groups <- names(x)

    if (ng.i > 1) {
         oask <- devAskNewPage(TRUE)
         on.exit(devAskNewPage(oask))
        }

    jg.i <- 0
    for (jg.i in 1:ng.i) {
         plot(x[[jg.i]], which.more = which.more, id.n = id.n, labels.id = labels.id,
              cex.id = cex.id, label.pos = label.pos, cex.caption = cex.caption,
              Main = paste("Group: ", groups[[jg.i]], sep=""), ...)
        }
}

# drop.gvf.points
drop.gvf.points.gvf.fit.gr <- function(x, method = c("pick", "cut"), which.plot = 1:2,
                                       res.type = c("standard", "student"), res.cut = 3,
                                       id.n = 3, labels.id = NULL,
                                       cex.id = 0.75, label.pos = c(4, 2),
                                       cex.caption = 1, col = NULL, drop.col = "red", ...)
{
  method <- match.arg(method)
  xd <- x
  groups <- names(x)

  jm <- 0
  for (m in x) {
     jm <- jm + 1

     if (method == "pick") {
         cat("\n## Processing Group: ", groups[[jm]], sep="")
         cat("\n## ")
        }
     if (method == "cut") {
         cat("\n## Processing Group: ", groups[[jm]], sep="")
         cat("\n##\n")
         dev.new()
        }

     xd[[jm]] <- drop.gvf.points(x = m, method = method, which.plot = which.plot,
                                 res.type = res.type, res.cut = res.cut,
                                 id.n = id.n, labels.id = labels.id,
                                 cex.id = cex.id, label.pos = label.pos,
                                 cex.caption = cex.caption, col = col, drop.col = drop.col,
                                 Main = paste("Group: ", groups[[jm]], sep=""), ...)
    }
  cat("\n")

  xd
}

# AIC
AIC.gvf.fit.gr <- function(object, ...)
{
lapply(object, AIC, ...)
}
AIC.gvf.fits.gr <- AIC.gvf.fit.gr

# BIC
BIC.gvf.fit.gr <- function(object, ...)
{
lapply(object, BIC, ...)
}
BIC.gvf.fits.gr <- BIC.gvf.fit.gr

# getR2
getR2.gvf.fit.gr <- function(object, adjusted = FALSE, ...)
{
lapply(object, FUN = getR2, adjusted = adjusted, ...) 
}
getR2.gvf.fits.gr <- getR2.gvf.fit.gr

# getBest
getBest.gvf.fit.gr <- function(object,
                               criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
#####################################
# gvf.fit method: no actual choice, #
# returns object.                   #
#####################################
{
criterion <- match.arg(criterion)
object
}
getBest.gvf.fits.gr <- function(object,
                                criterion = c("R2", "adj.R2", "AIC", "BIC"), ...)
###################################################
# gvf.fits.gr method: returns the fitted GVF with #
# highest average score in the given criterion    #
# over the groups.                                #
###################################################
{
criterion <- match.arg(criterion)
scoreFUN <- switch(criterion, "R2"     = function (object, ...) getR2(object, ...),
                              "adj.R2" = function (object, ...) getR2(object, adjusted = TRUE, ...),
                              "AIC"    = function (object, ...) -AIC(object, ...),
                              "BIC"    = function (object, ...) -BIC(object, ...)
                )
score <- rowMeans(sapply(object, scoreFUN, ...))
m <- which.max(score)
for (i.g in 1:length(object)) {
     object[[i.g]] <- object[[i.g]][[m]]
     }
class(object)[1] <- "gvf.fit.gr"
object
}

# predictCV
predictCV.gvf.fit.gr <- function(object, new.Y = NULL,
                                 scale = NULL, df = Inf,
                                 interval = c("none", "confidence", "prediction"),
                                 level = 0.95, na.action = na.pass,
                                 pred.var = NULL, weights = 1)
{
gvf.input.list <- lapply(object, function(x) attr(x, "gvf.input"))
group.vars <- attr(object, "group.vars")

if (!is.null(new.Y)){

     # If has only 'Y' column, predict for *all* groups
     ncol.new.Y <- ncol(new.Y)
     has.Y <- any(names(new.Y) == "Y")
     if (has.Y && (ncol.new.Y == 1)) {
         group.df <- Reduce(x = strsplit(names(gvf.input.list), ".", fixed = TRUE), f = rbind)
         rownames(group.df) <- NULL
         colnames(group.df) <- group.vars
         group.df <- as.data.frame(group.df)
         # Expand group.df rows so as to match new.Y values
         group.df <- group.df[rep(1:nrow(group.df), each = nrow(new.Y)), , drop = FALSE]
         new.Y <- data.frame(new.Y, group.df)
         rownames(new.Y) <- NULL
        }

     absent.vars <- group.vars[!(group.vars %in% names(new.Y))]
     if (length(absent.vars) > 0) {
         stop("Group variables not found in 'new.Y': ", paste(absent.vars, collapse = ", "))
        }
    }

# By default new.Y is the fitted gvf.input.gr object. This will return the CVs
# obtained by transforming the fitted response values
if (is.null(new.Y)){
     new.Y.list <- gvf.input.list
    }
else {
     new.Y.list <- split(new.Y, new.Y[group.vars], drop = TRUE)
    }
# Returning it for test, must finalize!!!!!
# new.Y.list

# 1. Subset 'object' so as to retain only models which are mapped to
#    realized groups in 'new.Y'
# 2. Subset new.Y.list so as to retain only components (groups) which are
#    mapped to retained models in 'object'
obj.retain <- (names(object) %in% names(new.Y.list))
new.Y.list.retain <- (names(new.Y.list) %in% names(object))
# subset
object <- object[obj.retain]
new.Y.list <- new.Y.list[new.Y.list.retain]

# Check that input estimates actually match at least one fitted group
if (length(new.Y.list) == 0) {
     stop("Input groups in 'new.Y' do not match any fitted group in 'object'") 
    }

# Compute CV predictions
pp <- mapply(FUN = predictCV, object, new.Y.list, SIMPLIFY = FALSE,
             MoreArgs = list(scale = scale, df = df, interval = interval,
                             level = level, na.action = na.action,
                             pred.var = pred.var, weights = weights)
            )

# Reduce the group list to a single data.frame
p.CV <- Reduce(f = rbind, x = pp)
rownames(p.CV) <- NULL
return(p.CV)
}
predictCV.gvf.fits.gr <- function(object, new.Y = NULL,
                                  scale = NULL, df = Inf,
                                  interval = c("none", "confidence", "prediction"),
                                  level = 0.95, na.action = na.pass,
                                  pred.var = NULL, weights = 1)
{
gvf.input.list <- lapply(object, function(x) attr(x, "gvf.input"))
group.vars <- attr(object, "group.vars")

if (!is.null(new.Y)){

     # If has only 'Y' column, predict for *all* groups
     ncol.new.Y <- ncol(new.Y)
     has.Y <- any(names(new.Y) == "Y")
     if (has.Y && (ncol.new.Y == 1)) {
         group.df <- Reduce(x = strsplit(names(gvf.input.list), ".", fixed = TRUE), f = rbind)
         rownames(group.df) <- NULL
         colnames(group.df) <- group.vars
         group.df <- as.data.frame(group.df)
         # Expand group.df rows so as to match new.Y values
         group.df <- group.df[rep(1:nrow(group.df), each = nrow(new.Y)), , drop = FALSE]
         new.Y <- data.frame(new.Y, group.df)
         rownames(new.Y) <- NULL
        }

     absent.vars <- group.vars[!(group.vars %in% names(new.Y))]
     if (length(absent.vars) > 0) {
         stop("Group variables not found in 'new.Y': ", paste(absent.vars, collapse = ", "))
        }
    }

# By default new.Y is the fitted gvf.input.gr object. This will return the CVs
# obtained by transforming the fitted response values
if (is.null(new.Y)){
     new.Y.list <- gvf.input.list
    }
else {
     new.Y.list <- split(new.Y, new.Y[group.vars], drop = TRUE)
    }
# Returning it for test, must finalize!!!!!
# new.Y.list

# Subset 'object' so as to retain only models which are mapped to
# realized groups in 'new.Y' 
object <- object[names(object) %in% names(new.Y.list)]

# Compute CV predictions
pp <- mapply(FUN = predictCV, object, new.Y.list, SIMPLIFY = FALSE,
             MoreArgs = list(scale = scale, df = df, interval = interval,
                             level = level, na.action = na.action,
                             pred.var = pred.var, weights = weights)
            )

# Here pp structure is [[group]][[model]]: must put all groups in each model
# into a single data.frame:
nmod <- length(pp[[1L]])
p.CV <- vector(mode = "list", length = nmod)
for (j in 1:nmod) {
     pp.j <- lapply(pp, function(el) el[[j]])
     p.CV[[j]] <- Reduce(f = rbind, x = pp.j)
     rownames(p.CV[[j]]) <- NULL
    }

return(p.CV)
}

# predict
predict.gvf.fit.gr <- function(object, ...)
{
lapply(object, predict, ...)
}
predict.gvf.fits.gr <- predict.gvf.fit.gr

# rstandard
rstandard.gvf.fit.gr <- function(model, ...)
{
lapply(model, rstandard, ...)
}
rstandard.gvf.fits.gr <- rstandard.gvf.fit.gr

# rstudent
rstudent.gvf.fit.gr <- function(model, ...)
{
lapply(model, rstudent, ...)
}
rstudent.gvf.fits.gr <- rstudent.gvf.fit.gr

# anova
anova.gvf.fit.gr <- function(object, ...)
{
lapply(object, anova, ...)
}
anova.gvf.fits.gr <- anova.gvf.fit.gr

# coef
coef.gvf.fit.gr <- function(object, ...)
{
lapply(object, coef, ...)
}
coef.gvf.fits.gr <- coef.gvf.fit.gr

# effects
effects.gvf.fit.gr <- function(object, ...)
{
lapply(object, effects, ...)
}
effects.gvf.fits.gr <- effects.gvf.fit.gr

# residuals
residuals.gvf.fit.gr <- function(object, ...)
{
lapply(object, residuals, ...)
}
residuals.gvf.fits.gr <- residuals.gvf.fit.gr

# fitted
fitted.gvf.fit.gr <- function(object, ...)
{
lapply(object, fitted, ...)
}
fitted.gvf.fits.gr <- fitted.gvf.fit.gr

# vcov
vcov.gvf.fit.gr <- function(object, ...)
{
lapply(object, vcov, ...)
}
vcov.gvf.fits.gr <- vcov.gvf.fit.gr
