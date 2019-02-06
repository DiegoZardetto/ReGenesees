#############################################
# Possible metrics for calibration trimming #
#############################################

# LINEAR -> the ONLY enabled, to date
`trim.linear` <-
list(
        Fm1 = function(u, bounds, sigma2){
              u <- u / sigma2
              pmin(pmax(u + 1, bounds[, 1]), bounds[, 2]) - 1
            },
        dF = function(u, bounds, sigma2){
             u <- u / sigma2
             as.numeric(u < bounds[, 2] - 1 & u > bounds[, 1] - 1) * (1 / sigma2)
            },
        name = "Generalized Euclidean Trimming"
    )
class(trim.linear) <- c("trimfun", "calfun")


# RAKING -> NOT enabled, to date
`trim.raking` <-
list(
        Fm1 = function (u, bounds, sigma2){
              u <- u / sigma2
              pmin(pmax(exp(u), bounds[, 1]), bounds[, 2]) - 1
            },
        dF = function (u, bounds, sigma2){
             u <- u / sigma2
             ifelse(u < bounds[, 2] - 1 & u > bounds[, 1] - 1, exp(u), 0) * (1 / sigma2)
            },
        name = "generalized raking trimming"
    )
class(trim.raking) <- c("trimfun", "calfun")


# LOGIT -> NOT enabled, to date
`trim.logit` <-
###############################################################
# If eAu goes to Inf, dF goes to 0: this was badly programmed #
# and has been fixed.                                         #
# NOTE: Currently DOES NOT WORK PROPERLY FOR TRIMMING...      #
#       ...MUST STILL INVESTIGATE!                            #
###############################################################
list(
        Fm1 = function (u, bounds, sigma2){
              u <- u / sigma2
              L <- bounds[, 1]
              U <- bounds[, 2]
              A <- (U - L)/((U - 1) * (1 - L))
              eAu <- exp(A * u)
              Fm1 <- (L * (U - 1) + U * (1 - L) * eAu)/(U - 1 + (1 - L) * eAu) - 1
              ifelse(is.finite(Fm1), Fm1, U-1)
            },
        dF = function (u, bounds, sigma2){
             u <- u / sigma2
             L <- bounds[, 1]
             U <- bounds[, 2]
             A <- (U - L)/((U - 1) * (1 - L))
             eAu <- exp(A * u)
             dF <- U*(1-L)*eAu*A/(U-1+(1-L)*eAu)-((L*(U-1)+U*(1-L)*eAu)*((1-L)*eAu*A))/(U-1+(1-L)*eAu)^2
             dF <- dF * (1 / sigma2)
             ifelse(is.finite(dF), dF, 0) 
            },
        name = "generalized exponential trimming"
    )
class(trim.logit) <- c("trimfun", "calfun")


# Print method
print.trimfun <- function(x,...) cat("Trimming Metric: ", x$name, "\n")
