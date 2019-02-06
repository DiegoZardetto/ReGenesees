`cal.linear` <-
list(
        Fm1 = function(u, bounds, sigma2){
              u <- u / sigma2
              pmin(pmax(u + 1, bounds[1]), bounds[2]) - 1
            },
        dF = function(u, bounds, sigma2){
             u <- u / sigma2
             as.numeric(u < bounds[2] - 1 & u > bounds[1] - 1) * (1 / sigma2)
            },
        name = "linear calibration"
    )
class(cal.linear) <- "calfun"
