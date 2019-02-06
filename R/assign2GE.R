##################################################
# Assign into the .GlobalEnv (in such a way that #
# R CMD check --as-cran is happy)                #
#                                                #
# NOTE: I don't like this, it makes me feel like #
#       I'm not actually programming in R but    #
#       rather in some CRAN dialect of the R     #
#       language...                              #
##################################################
`assign2GE` <- function(x, value){
     GE <- .GlobalEnv
     assign(x, value, envir = GE, inherits = FALSE, immediate = TRUE)
}
