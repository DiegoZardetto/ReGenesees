#######################################################
# - survey pkg problem:                               #
# WHEN COMPUTING ESTIMATES AND VARIANCES OF NONLINEAR #
# ESTIMATORS, svytotal is called to estimate values   #
# around which the Taylor linearization takes place:  #
# HERE THERE IS NO POINT IN ESTIMATING THE VARIANCE   #
# OF THE TOTAL (i.e. svyrecvar should not be called   #
# from svytotal). Lumley's solution is highly memory  #
# and CPU time hungry WITHOUT REASON!                 #
# - ReGenesees pkg solution:                          #
# A modified z.svytotal (generic, since cal.analytic  #
# and survey.design2 differ ONLY in variance          #
# estimation) is a good alternative.                  #
# HUGE SAVINGS OF MEMORY AND CPU TIME ARE OBTAINED!!  #
#######################################################

z.svytotal<-function(x,design,na.rm=FALSE,...){
  .svycheck(design)
  UseMethod("z.svytotal",design)
}

z.svytotal.survey.design2<-function(x,design, na.rm=FALSE, deff=FALSE,...){
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
    } else{
        if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        else {
            if(is.data.frame(x) && any(sapply(x,is.factor))){
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

    total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
    class(total)<-"svystat"
    attr(total,"statistic")<-"total"
    return(total)
}
