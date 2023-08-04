# References for method = "Beaumont"
# - [Beaumont 08]
# Beaumont, J. F. (2008). A new approach to weighting and inference in sample surveys. Biometrika, 95(3), 539-553.
# - [Beaumont, Rivest 09]
# Beaumont, J. F., Rivest, L. P. (2009). Dealing with outliers in survey data. In Handbook of statistics (Vol. 29, pp. 247-279). Elsevier.

smooth.strat.jump <- function(design, curr.strata, method = c("MinChange", "Beaumont")) {
#########################################################################
# Given a stratified one-stage unit sampling design object, smooth its  #
# weights to possibly mitigate the issue of unstable estimates that may #
# arise from stratum jumpers.                                           #
#                                                                       #
# Formula 'curr.strata' must identify the *current strata* variable,    #
# namely the strata variable measured at survey time by answers to      #
# questionnaire items. This is different from the *design strata*       #
# variable used to build 'design', which was measured at sampling time  #
# by means of the sampling frame.                                       #
#                                                                       #
# Stratum jumpers are, by definition, units whose design and current    #
# strata do differ. If the current value is reliable, stratum jumpers   #
# are evidence of frame imperfections (typically, the frame was not     #
# up-to-date). In enterprise surveys, stratum jumpers are often units   #
# that underwent a fast growth in size from sampling to survey time.    #
# Now, in most enterprise surveys the sampling design is such that      #
# smaller firms receive a smaller inclusion probability (and hence a    #
# larger weight). Therefore, stratum jumpers often have a "too large"   #
# weight, in the sense that they would have received a smaller weight   #
# had their actual size been known at sampling time. When these units   #
# also happen to have a large value of interest variable y, they may    #
# become influential in estimation. Weight smoothing tries to cope with #
# the problem of unstable (or even biased) estimators as a result of    #
# influential stratum jumpers.                                          #
#                                                                       #
# The smoothing process entails two steps:                              #
# 1) The weights are smoothed according to a given method (see below)   #
#    w -> w'                                                            #
# 2) The weights w' of all units are scaled by a global factor so as to #
#    preserve the initial overall sum of weights                        #
#    w' -> w'' = scale * w'  with  scale = sum(w)/sum(w')               #
#    so that sum(w'') = sum(w)                                          #
#                                                                       #
# Argument 'method' controls the smoothing algorithm in step 1).        #
# - MinChange (my proposal):                                            #
#   > Only the weights of stratum jumpers are smoothed, by setting      #
#     their value to the average weight of units belonging to the same  #
#     *current* stratum. Weights of all other units are left unchanged. #
# - Beaumont:                                                           #
#   > The weights of all units are smoothed, by setting their value to  #
#     the average weight of units belonging to the same *current*       #
#     stratum, with the exception of *minimum weight* units, whose      #
#     weights are left unchanged. Therefore all weights, exluding only  #
#     minimum weights within each *current* stratum, are smoothed. Note #
#     that this smoothing affects all *current* strata, even those that #
#     do NOT include any stratum jumper.                                #
#                                                                       #
# Put briefly: both methods lead to very similar smoothed weights for   #
# units that are stratum jumpers; however method "Beaumont" smooths the #
# weights of all other (i.e. non stratum jumpers) units much more       #
# aggressively than method "MinChange" (which alters them only in step  #
# 2 to preserve the overall sum of weights).                            #
#                                                                       #
#########################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum. (affects ONLY diagnostics written in the .GlobalEnv)
directly <- !( length(sys.calls()) > 1 )

# Restrict to un-calibrated designs
if (data.class(design) != "analytic") {
     stop("Design object must be of class analytic!")
    }

# Smoothing an already smoothed object must be avoided, in order to prevent
# 'smoothing cycles'
if (is.smoothed(design)) 
     stop("Weights of object 'design' have already been smoothed!")

# Check if design is actually stratified
if (!isTRUE(design$has.strata)) {
     stop("Design object must be stratified!")
    }

# Check if design is actually one-stage
ids.charvect <- all.vars(attr(design, "ids"))
if(length(ids.charvect) > 1){
     stop("Design object must be a one-stage sampling!")
    }

# Check if design is a unit sampling
n <- NROW(design$cluster)
un <- length(unique(design$cluster[,1]))
if(n != un){
     stop("Design object must be a unit sampling!")
    }

# Get the data
data <- design$variables

# Get *design strata* variable (i.e. measured at sampling time by the sampling 
# frame)... 
# ... paying attention to handle well *collapsed strata* design (DEBUG 9/12/2021)
is.clps <- !is.null(attr(design, "collapse.strata"))
# get strata variable name
des.strata.char <- all.vars(attr(design, "strata"))
# get strata variable values
if (!is.clps) {
     des.strata <- design$strata[, 1]
    } else {
     des.strata <- design$variables[[des.strata.char]]
    }

# Get the *current strata* variable (i.e. measured at survey time by answers to
# questionnaire items)
if (!inherits(curr.strata,"formula")) {
     stop("Current strata must be passed as a formula")
    }
curr.strata.char  <- all.vars(curr.strata)
if (length(curr.strata.char)<1) { 
     stop("Formula for current strata must reference a survey data variable")
    }
if (length(curr.strata.char)>1) {
     stop("Formula for current strata must reference only one variable")
    }

na.Fail(data, curr.strata.char)

if (!is.factor(data[, curr.strata.char])) { 
     stop("Current strata variable ", curr.strata.char, " is not a factor")
    }
curr.strata <- data[, curr.strata.char]

# Warn if design and current strata have different levels
des.strl <- levels(des.strata)
curr.strl <- levels(curr.strata)
if (!identical(des.strl, curr.strl)) {
     # get unmatched design strata
     des.unml <- des.strl[!(des.strl %in% curr.strl)]
     n.des.unml <- length(des.unml)
     # get unmatched current strata
     curr.unml <- curr.strl[!(curr.strl %in% des.strl)]
     n.curr.unml <- length(curr.unml)
     # build warning message
     warn <- ""
     if (n.des.unml > 0) warn <- paste(warn, "\n# Unmatched design strata found (",n.des.unml,"): ", paste(des.unml, collapse = ", "), "\n", sep = "")
     if (n.curr.unml > 0) warn <- paste(warn, "\n# Unmatched current strata found (",n.curr.unml,"): ", paste(curr.unml, collapse = ", "), "\n", sep = "")
     warning("Design strata and current strata factors have different levels. Please double-check!\n", warn)
    }

# Identify stratum jumpers...
# ... paying attention to handle well cases where des.strata and curr.strata have
#     different *number of levels* (DEBUG 9/12/2021), i.e.
#     - some des.strata disappeared at survey-time
#     - some curr.strata are brand new and did not exist at sampling-time
is.jumper <- ( as.character(curr.strata) != as.character(des.strata) )

# How many stratum jumpers?
n.jumpers <- sum(is.jumper)

# If no stratum jumpers found, raise an error. Else, print their number
if (n.jumpers > 0) {
     msg <- paste("Found ", n.jumpers, " stratum jumpers (out of ", NROW(design), " units), see strat.jump.status", sep="")
     cat(paste("\n# ", msg, "\n\n", sep=""))
    }
else {
     stop("No point in weight smoothing: no stratum jumpers found!")
    }

# Compute number of jumpers from design- to current-strata (for each pair)
data$nn.jumpers <- ave(is.jumper, des.strata, curr.strata, FUN = sum)

# Compute design- and current-strata averages of weights
weights.char <- all.vars(attr(design, "weights"))
w <- data[, weights.char]
data$des.str.w.avg <- ave(w, des.strata, FUN = mean)
data$curr.str.w.avg <- ave(w, curr.strata, FUN = mean)

# Build smoothed weights
w.smooth.char <- paste(weights.char, "smooth", sep = ".")
w.smooth.unsc.char <- paste(weights.char, "smooth.unsc", sep = ".")
data[[w.smooth.char]] <- data[[w.smooth.unsc.char]] <- w
method <- match.arg(method)
## SETP 1
# My minimum-change way:
if (method == "MinChange"){
     data[[w.smooth.char]][is.jumper] <- data$curr.str.w.avg[is.jumper]
    }

# Beaumont's way (in my view, a too aggressive and unnecessary smoothing):
if (method == "Beaumont"){
groups <- split(1:nrow(data), curr.strata)
n.groups <- length(groups)
group.names <- levels(curr.strata)
# Loop on current strata cells
# NOTE: A mitigation would be, unlike in Beaumont's approach, *to restrict to
#       current strata that are *actually affected* by the stratum jumpers
#       phenomenon (i.e. do contain at least one stratum jumper). In any case,
#       even with this modification, the approach would still seem to me way too
#       aggressive...
for (i in 1:n.groups) {
     group <- groups[[i]]
     w.avg <- mean(w[group])
     is.w.min <- w[group] <= min(w[group]) 
     data[group, w.smooth.char][!is.w.min] <- w.avg
    }
}

# Up to now, unscaled smoothed weights
# IMPORTANT REMARK: It could seem that unscaled smoothed weights for *stratum
#                   jumpers* must be *the same* for methods 'MinChange' and
#                   'Beaumont'. This is almost always true. However, it is false
#                   if any stratum jumper *happen to have the minimum weight*
#                   in its current stratum!
data[[w.smooth.unsc.char]] <- data[[w.smooth.char]]


## STEP 2
# Scale smoothed weights globally
scale <- sum(data[[weights.char]]) / sum(data[[w.smooth.char]])
data[[w.smooth.char]] <- scale * data[[w.smooth.char]]


# Build a diagnostic data frame to assess the impact of the smoothing process on
# stratum jumpers (these are the most relevant units to screen)
df.jumpers <- data[is.jumper, c(ids.charvect, weights.char, des.strata.char, "des.str.w.avg", curr.strata.char, "curr.str.w.avg", "nn.jumpers", w.smooth.unsc.char, w.smooth.char)]
colnames(df.jumpers) <- c("IDS", "W", "DES_STR", "DES_STR_W_AVG", "CURR_STR", "CURR_STR_W_AVG", "N_JUMP_DES_CURR_STR","W_SMOOTH_UNSC", "W_SMOOTH")

## Build output design object
# Add a jumper flag to design$variables:
design$variables$is.jumper <- is.jumper

# Add the smoothed weights to design:
## design$variables slot:
design$variables[[w.smooth.char]] <- data[[w.smooth.char]]
## design$prob slot:
design$prob <- 1/data[[w.smooth.char]]
## design attribute "weights":
attr(design, "weights") <- as.formula(paste("~", w.smooth.char, sep = ""), env = .GlobalEnv)
# NOTE: All ReGenesees function that *change* the weights of design objects (e.g.
#       e.calibrate, trimcal, ...) do *not* update the *design$allprob* slot.
#       Therefore, *design$allprob* acts as a persistent memory of the *initial*
#       weights in arbitrary multi-step weights adjustment procedures!

# Add a token to testify smoothing
# NOTE: THIS TOKEN COULD (AND MUST) BE REMOVED BY ANY SUBSEQUENT CALL OF
#       e.calibrate
attr(design, "smoothed") <- TRUE

# Store the current strata variable
attr(design, "curr.strata") <- curr.strata.char

# Store the smoothing method
attr(design, "smoothing.method") <- method

# Get the current call
design$call <- sys.call() # MUST check the right way, for possible GUI development

# Save the diagnostic data frame of stratum jumpers
assign2GE("strat.jump.status", df.jumpers)

design
}


is.smoothed <- function(design){
  isTRUE(attr(design, "smoothed"))
}
