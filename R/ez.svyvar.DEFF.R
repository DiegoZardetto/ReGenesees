#######################################################
# - survey pkg problem:                               #
# WHEN COMPUTING THE DEFF, svymean gets called from   #
# svyvar (generic) and here THERE IS NO POINT IN      #
# ESTIMATING THE VARIANCE OF THE MEAN (i.e. svyrecvar #
# should not be called). Lumley's solution is highly  #
# memory and CPU time hungry WITHOUT REASON!          #
# - ReGenesees pkg solution:                          #
# A modified z.svyvar which calls a modified function #
# z.svymean (generic, since survey.design2 and        #
# cal.analytic objects differ ONLY in variance        #
# estimation) is a good alternative.                  #
# HUGE SAVINGS OF MEMORY AND CPU TIME ARE OBTAINED!!  #
#######################################################

z.svyvar <- function (x, design, na.rm = FALSE, ...) 
{
    if (inherits(x, "formula")) 
        x <- model.frame(x, design$variables, na.action = na.pass)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, design$variables)
    ## If z.svyvar has been called on a subset get n from the domain index
    ## else compute it for the whole sample:
    if (is.null(di <- attr(design, "domain.index"))) n <- NROW(x) else n <- length(di)
    xbar <- z.svymean(x, design, na.rm = na.rm)
    if (NCOL(x) == 1) {
        x <- x - xbar
        v <- z.svymean(x * x * n/(n - 1), design, na.rm = na.rm)
        attr(v, "statistic") <- "variance"
        return(v)
    }
    x <- t(t(x) - xbar)
    p <- NCOL(x)
    a <- matrix(rep(x, p), ncol = p * p)
    b <- x[, rep(1:p, each = p)]
    v <- z.svymean(a * b * n/(n - 1), design, na.rm = na.rm)
    v <- matrix(v, ncol = p)
    attr(v, "statistic") <- "variance"
    v
}

z.svymean <- function(x,design, na.rm=FALSE,deff=FALSE,...){
################################################################
# Similar to the original svymean, without Variance estimation #
# which is not necessary for DEFF estimation.                  #
################################################################

  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  }
  else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else if(is.data.frame(x) && any(sapply(x,is.factor))){
          xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
          cols<-sapply(xx,NCOL)
          scols<-c(0,cumsum(cols))
          cn<-character(sum(cols))
          for(i in 1:length(xx))
              cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
          x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
          for(i in 1:length(xx)){
              x[,scols[i]+1:cols[i]]<-xx[[i]]
          }
          colnames(x)<-cn
      }
    }
  x<-as.matrix(x)

    ## Missing values treatment
    nas<-rowSums(is.na(x))
    if (na.rm && sum(nas)>0){
        # If domain has some non-NA values, use them for estimation:
        if (length(design[nas==0,]$prob) > 0) {
            design<-design[nas==0,]
            if (length(nas)>length(design$prob))
                x<-x[nas==0,,drop=FALSE]
            else
                x[nas>0,]<-0
        }
        # If domain has only NAs, cannot do anything (i.e. behave as na.rm=FALSE)
        else {
            na.rm <- FALSE
        }
    }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  ######################################################################
  # What follows would give the "standard" formula for the HT estimate #
  # of the variance of a total under a SRSWOR of units (with fpc and   #
  # population size "rescaled" if the design is calibrated).           #
  # SHOULD DISCUSS THIS SOLUTION!!!                                    # 
  ######################################################################
  ## If z.svyvar has been called on a subset get n from the domain index
  ## else compute it for the whole sample:
  #  if (is.null(di <- attr(design, "domain.index"))) n <- NROW(x) else n <- length(di)
  #  average<-colSums(x/n)
  return(average)
  }
