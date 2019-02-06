`FROM` <- function(fun, pkg)
#########################################################
# Retrieve function 'fun' from package 'pkg' namespace, #
# even if it is not exported.                           #
#########################################################
{
    get(as.character(fun), envir = asNamespace(as.character(pkg)), inherits = FALSE)
}




############################
# ReGenesees FROM register #
############################

#### pkg: stats
# in: gvf.R
`stats_plot.lm` <- FROM(fun = "plot.lm", pkg = "stats")

# in: gvf.R
`stats_rstandard.lm` <- FROM(fun = "rstandard.lm", pkg = "stats")

# in: gvf.R
`stats_rstudent.lm` <- FROM(fun = "rstudent.lm", pkg = "stats")

# in: gvf.R
`stats_anova.lm` <- FROM(fun = "anova.lm", pkg = "stats")

# in: gvf.R
`stats_coef.default` <- FROM(fun = "coef.default", pkg = "stats")

# in: gvf.R
`stats_effects.lm` <- FROM(fun = "effects.lm", pkg = "stats")

# in: gvf.R
`stats_fitted.default` <- FROM(fun = "fitted.default", pkg = "stats")

# in: gvf.R
`stats_vcov.lm` <- FROM(fun = "vcov.lm", pkg = "stats")

