`pop.plot` <- function(pop.totals, ...)
#################################################
# Generic function. Describe population totals. #
#################################################
{
UseMethod("pop.plot")
}

`pop.plot.pop.totals` <- function(pop.totals, design,
                                 xlab = "Current Estimates", ylab = "Calibration Control Totals",
                                 lcol = c("red", "green", "blue"), lwd = c(1, 1, 1), lty = c(2, 1, 2), verbose = TRUE, ...){
################################################################################
# Scatter plot of X.pop vs. X.hat, where X.pop are population benchmarks for a #
# calibration task and X.hat are the corresponding HT estimates using current  #
# design weights.                                                              #
# This is intended as a diagnostic tool to have a quick and dirty assessment   #
# of the expected optimization complexity of the calibration task, and to spot #
# potential issues (e.g. bias) in the underlying auxiliary sample data.        #
################################################################################
# Get HT estimates (and delegate consistency checks on arguments)
HT.totals <- aux.estimates(design, template = pop.totals)

# Is it a global calibration task?
partition <- attr(pop.totals, "partition")
is.global <- identical(partition, FALSE)
n.part.vars <- 0
aux.totals <- pop.totals
aux.HT <- HT.totals
# If calibration task is partitioned, get rid of partition variables
if (!is.global) {
     partition.vars <- all.vars(partition)
     n.part.vars <- length(partition.vars)
     # Drop partition columns to restrict to auxiliary info columns
     aux.totals <- pop.totals[, -(1:n.part.vars), drop = FALSE]
     aux.HT <- aux.HT[, -(1:n.part.vars), drop = FALSE]
    }

# Coerce to numeric vectors
aux.totals.v <- as.numeric(unlist(aux.totals))
aux.HT.v <- as.numeric(unlist(aux.HT))

# xlim <- ylim <- range(aux.totals.v, aux.HT.v)

# Take care of simultaneous zeros (which would be skipped in calibration) 
zero.over.zero <- (aux.totals.v == 0) & (aux.HT.v == 0)
aux.totals.v <- aux.totals.v[!zero.over.zero]
aux.HT.v <- aux.HT.v[!zero.over.zero]

# Get g-style ratios and slopes of envelope bounds
g <- aux.totals.v / aux.HT.v
g.range <- range(g)
min.g <- g.range[1]
max.g <- g.range[2]

# If the calibration task is *not* special purpose
# if (data.class(pop.totals) == "pop.totals") {
     # Draw the scatter plot
     plot(aux.HT.v, aux.totals.v, xlab = xlab, ylab = ylab, ...)
     # Draw the 3 lines
     if (is.finite(min.g)){
         abline(c(0, min.g), col = lcol[1], lwd = lwd[1], lty = lty[1], untf = TRUE)
        }
     if (is.finite(max.g)){
         abline(c(0, max.g), col = lcol[3], lwd = lwd[3], lty = lty[3], untf = TRUE)
        }
     abline(0:1, col = lcol[2], lwd = lwd[2], lty = lty[2], untf = TRUE)
     # Print on screen
     if (isTRUE(verbose)) {
         cat("\nMinimum and maximum slopes: [", round(min.g, 3), ", ", round(max.g, 3), "]\n\n", sep = "")
        }
     # Return
     return(invisible(g.range))
#     }

# # This should not happen
# return(invisible(NULL))
}

`pop.plot.spc.pop` <- function(pop.totals, design,
                               xlab = "Current Estimates", ylab = "Calibration Control Totals",
                               lcol = c("red", "green", "blue"), lwd = c(1, 1, 1), lty = c(2, 1, 2), verbose = TRUE, ...){
################################################################################
# Scatter plot of X.pop vs. X.hat, where X.pop are population benchmarks for a #
# calibration task and X.hat are the corresponding HT estimates using current  #
# design weights.                                                              #
# This is intended as a diagnostic tool to have a quick and dirty assessment   #
# of the expected optimization complexity of the calibration task, and to spot #
# potential issues (e.g. bias) in the underlying auxiliary sample data.        #
################################################################################
# Get HT estimates (and delegate consistency checks on arguments)
HT.totals <- aux.estimates(design, template = pop.totals)

# Is it a global calibration task?
partition <- attr(pop.totals, "partition")
is.global <- identical(partition, FALSE)
n.part.vars <- 0
aux.totals <- pop.totals
aux.HT <- HT.totals
# If calibration task is partitioned, get rid of partition variables
if (!is.global) {
     partition.vars <- all.vars(partition)
     n.part.vars <- length(partition.vars)
     # Drop partition columns to restrict to auxiliary info columns
     aux.totals <- pop.totals[, -(1:n.part.vars), drop = FALSE]
     aux.HT <- aux.HT[, -(1:n.part.vars), drop = FALSE]
    }

# Coerce to numeric vectors
aux.totals.v <- as.numeric(unlist(aux.totals))
aux.HT.v <- as.numeric(unlist(aux.HT))

# xlim <- ylim <- range(aux.totals.v, aux.HT.v)

# Take care of simultaneous zeros (which would be skipped in calibration) 
zero.over.zero <- (aux.totals.v == 0) & (aux.HT.v == 0)
aux.totals.v <- aux.totals.v[!zero.over.zero]
aux.HT.v <- aux.HT.v[!zero.over.zero]

# Get g-style ratios and slopes of envelope bounds
g <- aux.totals.v / aux.HT.v
g.range <- range(g)
min.g <- g.range[1]
max.g <- g.range[2]

# If the calibration task is a special purpose one
# if (data.class(pop.totals) == "spc.pop") {
     # If it is *not* fused
     if (!is.fused.spc(pop.totals)) {
         # Draw the scatter plot
         plot(aux.HT.v, aux.totals.v, xlab = xlab, ylab = ylab, ...)
         # plot one vertical green line
         abline(v = 0, col = lcol[2], lwd = lwd[2], lty = lty[2], untf = TRUE)
         # Return NULL invisibly (no meaningful g.range)
         return(invisible(NULL))
        }
     # If it is fused
     else {
         # Identify sub-dataframes of ordinary and special purpose benchmarks
         aux.origin <- attr(pop.totals, "aux.origin")
         # If calibration task is partitioned, get rid of partition variables
         if (!is.global) {
                         aux.origin <- attr(pop.totals, "aux.origin")[-(1:n.part.vars)]
                        }
         ord.totals <- aux.totals[, aux.origin == "pop", drop = FALSE]
         # ord.totals <- ord.totals[, -(1:n.part.vars), drop = FALSE]
         spc.totals <- aux.totals[, aux.origin == "spc", drop = FALSE]
         # Identify sub-dataframes of ordinary and special purpose HT estimates
         ord.HT <- aux.HT[, aux.origin == "pop", drop = FALSE]
         # ord.HT <- ord.HT[, -(1:n.part.vars), drop = FALSE]
         spc.HT <- aux.HT[, aux.origin == "spc", drop = FALSE]
         # Coerce to numeric vectors
         ## Ordinary ##
         ord.totals.v <- as.numeric(unlist(ord.totals))
         ord.HT.v <- as.numeric(unlist(ord.HT))
         # Take care of simultaneous zeros (which would be skipped in calibration) 
         zero.over.zero <- (ord.totals.v == 0) & (ord.HT.v == 0)
         ord.totals.v <- ord.totals.v[!zero.over.zero]
         ord.HT.v <- ord.HT.v[!zero.over.zero]
         # Get g-style ratios and slopes of envelope bounds
         g <- ord.totals.v / ord.HT.v
         g.range <- range(g)
         min.g <- g.range[1]
         max.g <- g.range[2]
         ## Spc ##
         spc.totals.v <- as.numeric(unlist(spc.totals))
         spc.HT.v <- as.numeric(unlist(spc.HT))
         # Plot
         opar <- par(mfrow = c(1, 2))
         # Ordinary plot
         # Draw the scatter plot
         plot(ord.HT.v, ord.totals.v, xlab = xlab, ylab = ylab, main = "Ordinary Calibration", ...)
         # Draw the 3 lines
         if (is.finite(min.g)){
             abline(c(0, min.g), col = lcol[1], lwd = lwd[1], lty = lty[1], untf = TRUE)
            }
         if (is.finite(max.g)){
             abline(c(0, max.g), col = lcol[3], lwd = lwd[3], lty = lty[3], untf = TRUE)
            }
         abline(0:1, col = lcol[2], lwd = lwd[2], lty = lty[2], untf = TRUE)
         # Spc plot
         # Draw the scatter plot
         plot(spc.HT.v, spc.totals.v, xlab = xlab, ylab = ylab, main = "Special Purpose Calibration", ...)
         # plot one vertical green line
         abline(v = 0, col = lcol[2], lwd = lwd[2], lty = lty[2], untf = TRUE)
         # Reset par
         par(opar)
         # Print on screen
         if (isTRUE(verbose)) {
             cat("\nMinimum and maximum slopes: [", round(min.g, 3), ", ", round(max.g, 3), "]\n\n", sep = "")
            }
         # Return
         return(invisible(g.range))
        }
#     }

# # This should not happen
# return(invisible(NULL))
}
