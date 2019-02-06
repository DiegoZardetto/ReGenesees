

ftable.svystat<-function(x, rownames=NULL, ...){

    m<-cbind(coef(x),SE(x))
    if (is.null(rownames))
        return(as.table(m))

    statname<-if (is.list(x)) attr(x[[1]],"statistic") else attr(x,"statistic")

    deff<-attr(x,"deff")
    has.deff<-!is.null(deff)
    if (has.deff)
      m<-cbind(m,diag(deff))
    
    rowdim<-sapply(rownames,length)

    if (has.deff){
      mm<-array(m,dim=c(rowdim,NCOL(m)),
                dimnames=c(as.list(rownames),
                  list(c(statname,"SE","Deff"))))
      
      ftable(mm,row.vars=length(rowdim)+0:1)
    } else {
      mm<-array(m,dim=c(rowdim,NCOL(m)),
                dimnames=c(as.list(rownames),
                  list(c(statname,"SE"))))
     
      ftable(mm,row.vars=length(rowdim)+0:1)
    }

}

ftable.svrepstat<-ftable.svystat


ftable.svyby <- function (x, ...)
###############################################################
# MODIFIED things:                                            #
# 1. Whenever the svyby object has been built with the option #
#    'drop.empty.groups = TRUE' .and. there are effectively   #
#    empty groups in data, the ftable.svyby output was        #
#    WRONGLY FORMATTED (due to the wrong way in wich          #
#    statistics and variances filled the array rval).         #
###############################################################
{
    info <- attr(x, "svyby")
    margins <- info$margins
    dimnames <- lapply(x[, margins, drop = FALSE], levels)
    dims <- sapply(dimnames, length)
    dims <- c(dims, variable = info$nstats)
    senames<-c(se="SE",cv="cv",cvpct="cv%",var="Var")[info$vartype]
    if (info$vars || info$deffs) {
        dims <- c(dims, 1 + info$vars + info$deffs)
        dimnames <- c(dimnames,
                      list(sub("^statistic\\.(.*)$", "\\1", info$variables)),
                      list(c(info$statistic,
                             if (info$vars) senames,
                             if (info$deffs) "DEff")))
      }
    else if (info$nstats == 1) {
        dimnames <- c(dimnames, list(info$statistic))
    }
    else {
        dimnames <- c(dimnames, list(info$variables))
    }
    ## fix by Sergio Calva for ordering bug.
    x <- x[do.call("order",x[,rev(margins),drop = FALSE]),]
#################
# MODIFIED part #
#################
    ## Build empty array with right dim and dimnames
    rval <- array(NA, dim = dims, dimnames = dimnames)
    ## if there were empty groups in svyby object, fill rval in the right way
    margins.tab <- table(x[, margins])
    rval[margins.tab!=0] <- as.matrix(x[, -margins, drop = FALSE])
#################
# END           #
#################
    ftable(rval, row.vars = c(1, length(dim(rval))))
}
