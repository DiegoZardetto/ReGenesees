`cal.raking` <-
list(
        Fm1=function (u, bounds, sigma2){
            u <- u / sigma2
            pmin(pmax(exp(u), bounds[1]), bounds[2]) - 1
        },
        dF= function (u, bounds, sigma2){
            u <- u / sigma2
            ifelse(u < bounds[2] - 1 & u > bounds[1] - 1, exp(u), 0) * (1 / sigma2)
        },
        name="raking"
    )
class(cal.raking) <- "calfun"
