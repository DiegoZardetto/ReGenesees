\name{write.svystat}
\alias{write.svystat}
\title{
Export Survey Statistics
}
\description{
Prints survey statistics to a file or connection.
}
\usage{
write.svystat(x, ...)
}
\arguments{
  \item{x}{An object containing survey statistics.}
  \item{\dots}{Arguments to \code{write.table}}
}
\details{
This function is just a convenience wrapper to \code{\link{write.table}}, designed to export objects which have been returned by survey statistics functions (e.g. \code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}}, \code{\link{svySigma}}, \code{\link{svySigma2}}).
}
\author{
Diego Zardetto
}
\seealso{
\code{\link{write.table}} and the 'R Data Import/Export' manual. 
}
\examples{
# Creation of a design object:
data(sbs)
des<-e.svydesign(data=sbs,ids=~id,strata=~strata,weights=~weight,
     fpc=~fpc)

# Estimation of the average value added per employee
# for economic activity region and macro-sectors,
# with SE, CV% and standard confidence intervals:
stat <- svystatR(des,~va.imp2,~emp.num,by=~region:nace.macro,
        vartype=c("se","cvpct"),conf.int=TRUE)
stat

# In order to export the summary statistics above
# into a CSV file for input to Excel one can use:
\dontrun{
write.svystat(stat,file="stat.csv",sep=";")
}

# ...and to read this file back into R one needs
\dontrun{
stat.back <- read.table("stat.csv",header=TRUE,sep=";",
             check.names=FALSE)
stat.back
}

# Notice, however, that the latter object has
# lost a lot of meta-data as compared to the
# original one, so that e.g.:
\dontrun{
confint(stat.back)
}

# ...while, on the contrary:
confint(stat)
}