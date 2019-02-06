.onLoad<-function(libname, pkgname){
  # Standard ReGenesees variance estimation options
  if (is.null(getOption("RG.lonely.psu")))
    options(RG.lonely.psu="fail")
  if (is.null(getOption("RG.ultimate.cluster")))
    options(RG.ultimate.cluster=FALSE)
  if (is.null(getOption("RG.adjust.domain.lonely")))
    options(RG.adjust.domain.lonely=FALSE)
  if (is.null(getOption("RG.warn.domain.lonely")))
    options(RG.warn.domain.lonely=FALSE)
  # Standard ReGenesees contrasts (i.e. no specific contrasts for ordered factors)
    options("contrasts"=c(unordered="contr.treatment",ordered="contr.treatment"))
}
