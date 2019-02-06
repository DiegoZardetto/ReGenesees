#######################################################################
# Peter Dalgaard: "[...] the corner cases of model.matrix and friends #
# is some of the more impenetrable code in the R sources [...]"       #
#######################################################################

`contr.off` <-
 function(n, base = 1, contrasts = TRUE, sparse = FALSE){
####################################################################
# This function serves the purpose of switching off all contrasts  #
# treatments when using ReGenesees.                                #
# NOTE: Contrasts are generally useful when calibrating, as they   #
#       allow to drop linearly dependent auxiliary variables and   #
#       the related population totals. Anyway, there could be some #
#       user who may prefer to work with complete (albeit          #
#       redundant) population totals and auxiliary variables model #
#       matrices. Besides this, the standard R contrasts treatment #
#       can sometimes turn out to be not adequate under the        #
#       calibration perspective, especially for intricate          #
#       calibration models involving the interactions of a numeric #
#       variable (say X) and a factor variable (say A) with the    #
#       SAME factor variable (say D).                              #
####################################################################
  contr.treatment(n = n, base = base, contrasts = FALSE, sparse = sparse)
}


# ReGenesees user visible functions to switch off/on contrasts treatment
  # 1) ReGenesees STANDARD (i.e. no specific contrasts for ordered factors)
`contrasts.RG` <- function(){
  options("contrasts"=c(unordered="contr.treatment",ordered="contr.treatment"))
  cat("\n# Standard ReGenesees contrasts treatment has been set:\n\n")
  print(options("contrasts"))
}

  # 2) SWITCH OFF (!!!! Da stats non viene trovata la funzione: MUST FIX !!!!!!)
`contrasts.off` <- function(){
  options("contrasts"=c(unordered="contr.off",ordered="contr.off"))
  cat("\n# Contrasts treatment has been switched off:\n\n")
  print(options("contrasts"))
}

  # 3) RESTORE DEFAULTS  
`contrasts.reset` <- function(){
  options("contrasts"=c(unordered="contr.treatment",ordered="contr.poly"))
  cat("\n# Factory-fresh defaults for contrasts treatment have been restored:\n\n")
  print(options("contrasts"))
}
