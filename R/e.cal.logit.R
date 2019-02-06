`cal.logit` <-
###############################################################
# If eAu goes to Inf, dF goes to 0: this was badly programmed #
# and has been fixed.                                         #
###############################################################
list(
        Fm1=function (u, bounds, sigma2){
            u <- u / sigma2
            L <- bounds[1]
            U <- bounds[2]
            A <- (U - L)/((U - 1) * (1 - L))
            eAu <- exp(A * u)
            Fm1 <- (L * (U - 1) + U * (1 - L) * eAu)/(U - 1 + (1 - L) * eAu) - 1
            ifelse(is.finite(Fm1), Fm1, U-1)
        },
        dF= function (u, bounds, sigma2){
            u <- u / sigma2
            L <- bounds[1]
            U <- bounds[2]
            A <- (U - L)/((U - 1) * (1 - L))
            eAu <- exp(A * u)
            dF <- U*(1 - L)*eAu*A/(U-1+(1-L)*eAu)-((L*(U-1)+U*(1-L)*eAu)*((1-L)*eAu*A))/(U-1+(1-L)*eAu)^2
            dF <- dF * (1 / sigma2)
            ifelse(is.finite(dF), dF, 0) 
        },
        name="logit distance"
    )
class(cal.logit) <- "calfun"
