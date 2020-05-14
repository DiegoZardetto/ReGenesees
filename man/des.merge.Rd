\name{des.merge}
\alias{des.merge}

\title{Merge New Survey Data into Design Objects}

\description{
Modifies an analytic object by joining the original survey data with a new data frame via a common key.
}

\usage{
des.merge(design, data, key)
}

\arguments{
  \item{design}{Object of class \code{analytic} (or inheriting from it) containing survey data and sampling design metadata.}
  \item{data}{Data frame containing a key variable, plus new variables to be merged to \code{design} data.}
  \item{key}{Formula identifying the common key variable to be used for merging.}
}

\details{
This function updates the survey variables contained into \code{design} (i.e. \code{design$variables}), by merging the original data with those contained into the \code{data} data frame. The merge operation exploits a \emph{single} variable \code{key}, which must be common to both \code{design} and \code{data}.

The function preserves both the original \emph{ordering} of the survey data stored into \code{design}, as well as all the original sampling design \emph{metadata}.

The variable referenced by \code{key} must be a \emph{valid} unique key for both \code{design} and \code{data}: it must not contain duplicated values, nor \code{NA}s. Moreover, the values of \code{key} in \code{design} and \code{data} must be in \emph{1:1 correspondence}. These requirements are meant to ensure that the \emph{new} survey data (that is the merged ones) will have exactly the same number of rows as the \emph{old} survey data stored into \code{design}.

Should \code{design} and \code{data} contain further common variables besides the \code{key}, only their original \code{design} version will be retained. Thus, \code{des.merge} cannot modify any pre-existing \code{design} columns. This an intentional feature intended to safeguard the integrity of the relations between survey data and sampling design metadata stored in \code{design}.
}

\section{Practical Purpose}{
In the field of Official Statistics, it is not infrequent that calibration weights must be computed even several months before the target variables of the survey are made available for estimation. Such a time lag follows from the fact that target variables typically undergo much more thorough editing and imputation procedures than auxiliary variables.

In such production scenarios, function \code{des.merge} allows to tackle the task of computing estimates and errors for the fresh-released target variables \emph{without} any need of \emph{repeating} the calibration step. Indeed, by using the function, one can join the data contained into an already calibrated \code{design} object with new \code{data} made available only after the calibration step. The merge operation is made easy and safe, and preserves all the original calibration metadata (e.g. those needed for variance estimation).
}

\value{
An object of the same class of \code{design}, containing additional survey data but supplied with exactly the same metadata.
}

\author{Diego Zardetto}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata, \code{\link{e.calibrate}} for calibrating weights, \code{\link{des.addvars}} to add new variables to design objects.
}

\examples{
data(data.examples)

# Create a design object:
des<-e.svydesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
     weights=~weight)

# Create a calibrated design object as well (e.g. using population totals
# stored inside pop03p):
cal<-e.calibrate(design=des,df.population=pop03p,
                 calmodel=~marstat-1,partition=~sex,calfun="logit",
                 bounds=bounds)

# Lastly create a new data frame to be merged into des and cal:
set.seed(12345)    # RNG seed fixed for reproducibility
new.data<-example[,c("income","key")]
new.data$income <- 1000 + new.data$income    # altered income values
new.data$NEW.f<-factor(sample(c("A","B"),nrow(new.data),rep=TRUE))
new.data$NEW.n<-rnorm(nrow(new.data),10,2)
new.data <- new.data[sample(1:nrow(new.data)), ]    # rows ordering changed
head(new.data)

###########################################################
# Example 1: merge new data into a non calibrated design. #
###########################################################

# Merge new data inside des (note the warning on income):
des2<-des.merge(design=des,data=new.data,key=~key)

# Compare visually:
## before:
head(des$variables)
## after:
head(des2$variables)

# New data can be used as usual:
svystatTM(des2,~NEW.n,~NEW.f,vartype="cvpct")

# Old data are unaffected, as it must be:
svystatTM(des,~income,estimator="Mean",vartype="cvpct")
svystatTM(des2,~income,estimator="Mean",vartype="cvpct")

#######################################################
# Example 2: merge new data into a calibrated design. #
#######################################################

# Merge new data inside cal (note the warning on income):
cal2<-des.merge(design=cal,data=new.data,key=~key)

# Compare visually:
## before:
head(cal$variables)
## after:
head(cal2$variables)

# New data can be used as usual:
svystatTM(cal2,~NEW.n,~NEW.f,vartype="cvpct")

# Old data are unaffected, as it must be:
svystatTM(cal,~income,estimator="Mean",vartype="cvpct")
svystatTM(cal2,~income,estimator="Mean",vartype="cvpct")
}
\keyword{survey}