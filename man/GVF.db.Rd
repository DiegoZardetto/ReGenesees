\name{GVF.db}
\alias{GVF.db}

\title{
Archive of Registered GVF Models
}

\description{
\code{GVF.db} is the archive of \emph{registered} (i.e. built-in and/or user-defined) Generalized Variance Functions models supported by \pkg{ReGenesees}. Special accessor functions allow to customize, maintain, extend, update, save and reset such archive.
}

\format{
Each row of the \code{GVF.db} data frame represents a registered GVF model, with relevant information on the following 4 variables:
\describe{
    \item{\code{Model.id}}{A unique integer key identifying the GVF model, \code{integer}.}
    \item{\code{GVF.model}}{A character string specifying the GVF model formula, \code{character}. See also \sQuote{Details}.}
    \item{\code{Estimator.kind}}{A character string identifying the kind of estimators for which the GVF model is deemed to be appropriate, \code{character}. See also \sQuote{Details}.}
    \item{\code{Resp.to.CV}}{A character string which represents the function mapping the response of the GVF model (namely: variable 'resp') to the coefficient of variation (namely: variable 'CV'), \code{character}. See also \sQuote{Details}.}
  }
}

\usage{
GVF.db

\special{GVF.db$insert(GVF.model, Estimator.kind = NA, Resp.to.CV = NA, verbose = TRUE)}

\special{GVF.db$delete(Model.id, verbose = TRUE)}

\special{GVF.db$get(verbose = TRUE)}

\special{GVF.db$assign(value, verbose = TRUE)}

\special{GVF.db$reset(verbose = TRUE)}
}

\arguments{
  \item{GVF.model}{A GVF model, expressed as a formula object or as a character string (see \sQuote{Details}).}
  \item{Estimator.kind}{Character string identifying the kind of estimators for which the GVF model is deemed to be appropriate (see \sQuote{Details}).}
  \item{Resp.to.CV}{Character string representing the function which maps the response of the GVF model (namely: variable 'resp') to the coefficient of variation (namely: variable 'CV'), see \sQuote{Details}.}
  \item{Model.id}{Unique integer key identifying the GVF model.}
  \item{value}{An exported copy of \code{GVF.db}, as returned by \code{GVF.db$get()}.}
  \item{verbose}{Enables printing of a summary description of the result (the default is \code{TRUE}).}
}

\details{
\code{GVF.db} stores information about Generalized Variance Functions models supported by \pkg{ReGenesees}. When starting a new work session with \pkg{ReGenesees}, \code{GVF.db} contains few built-in GVF models (currently 5, see sections \sQuote{Source} and \sQuote{Examples}). The content of \code{GVF.db} can be customized by means of special accessor functions:
\preformatted{
  ACCESSOR FUNCTION       PURPOSE
  GVF.db$insert...........Register a new GVF model by adding a new row to the
                          GVF.db archive
  GVF.db$delete...........Unregister a GVF model by deleting the corresponding
                          row from GVF.db
  GVF.db$get..............Get the current version of GVF.db (e.g. to copy/save a
                          customized archive for later usage)
  GVF.db$assign...........Overwrite the current version of GVF.db (e.g. to use a
                          customized archive which was exported in a previous
                          ReGenesees session)
  GVF.db$reset............Reset GVF.db to its default version (i.e. the one with
                          built-in GVF models only)
}

Information about registered GVF models stored inside \code{GVF.db} will be accessed and used by \pkg{ReGenesees} Generalized Variance Functions facilities, e.g. functions \code{\link{fit.gvf}} or \code{\link{predictCV}}.\cr\cr

\strong{GVF.db$insert()}\cr
Function \code{GVF.db$insert} has just a single mandatory argument: \code{GVF.model}. This can be either a two-sided formula or a character string which would be transformed into a (well formed) two-sided formula by function \code{\link{as.formula}}.

The \code{GVF.model} formula to be inserted into \code{GVF.db} must be \emph{new} (i.e. not already present into the archive) and can involve only variables contained inside \code{\link{gvf.input}} objects, namely:
\preformatted{
 (1)  'Y'
 (2)  'SE'
 (3)  'CV'
 (4)  'VAR'
 (5)  'DEFF'
 }
Moreover, since GVF models are intended to model variances in terms of estimates, the response term of \code{GVF.model} must involve some of \code{'SE'}, \code{'CV'}, \code{'VAR'}, and the linear predictor must involve \code{'Y'}.

Optional argument \code{Estimator.kind} can be used to specify the kind of estimators for which the \code{GVF.model} is deemed to be appropriate. There are currently only 11 valid values for \code{Estimator.kind}, namely:
\preformatted{
 (1)  'Total'
 (2)  'Mean'
 (3)  'Frequency'
 (4)  'Absolute Frequency'
 (5)  'Relative Frequency'
 (6)  'Ratio'
 (7)  'Share'
 (8)  'Share Ratio'
 (9)  'Regression Coefficient'
(10)  'Quantile'
(11)  'Complex Estimator'
 }
Note that category \code{'Frequency'} has to be understood as an aggregation of categories \code{'Absolute Frequency'} and \code{'Relative Frequency'}, thus being appropriate for GVF models which are deemed to work well for estimators of \emph{both} kind of frequencies.

One of the primary motivations for building and fitting a GVF model is to exploit the fitted model to \emph{predict} the sampling error associated to a given estimate, instead of having to \emph{compute} directly an estimate of such sampling error. Optional argument \code{Resp.to.CV} is relevant to that scope.

Indeed, different GVF models can actually specify as response term (call it \code{'resp'} for definiteness) different functions of variables \code{'SE'}, \code{'CV'}, and \code{'VAR'}, but \pkg{ReGenesees} will always adopt variable \code{'CV'} as a \emph{pivot}. Thus, when registering a new GVF model, the user can provide via argument \code{Resp.to.CV} the function which transforms the response of the model, \code{'resp'}, into the pivot measure of variability, \code{'CV'}. A look to the default content of \code{GVF.db} should make the latter statement clear (see \sQuote{Examples}).

Note that while \code{Resp.to.CV} is passed as a character string, that string is expected to represent a well-formed mathematical expression (otherwise function \code{\link{predictCV}} would not work). Moreover, only variables \code{'resp'} and \code{'Y'} are allowed to appear inside \code{Resp.to.CV} (which is enough, since \code{'VAR'} and \code{'SE'} can be expressed in terms of \code{'CV'} and \code{'Y'}).

If the user does not specify \code{Resp.to.CV} when registering a new GVF model, he will be not able to use function \code{\link{predictCV}} for predicting CV values based on the fitted GVF model.

Lastly, note that the \code{Model.id} of a newly inserted GVF model will automatically be set, by adding 1 to the previous maximum of \code{Model.id}.\cr\cr


\strong{GVF.db$delete()}\cr
Function \code{GVF.db$delete} has just a single mandatory argument: \code{Model.id}. It must match the integer key of the (already existing) GVF model you want to drop from \code{GVF.db}.

Note that, after deleting a GVF model from \code{GVF.db}, values of column \code{Model.id} will be automatically renumbered, so as to range always from \code{1} to \code{nrow(GVF.db)}.\cr\cr


\strong{GVF.db$get()}\cr
Function \code{GVF.db$get} has no mandatory arguments. When invoked, the function returns the \emph{current} content of \code{GVF.db}, so that it can be assigned and saved/exported for later usage (see \sQuote{Examples}). Should the \emph{current} content of \code{GVF.db} happen to be empty, the function would inform the user and return \code{NULL}. The return value of \code{GVF.db$get} has class \code{"GVF.db_exported"}, and inherits from class \code{"data.frame"}.\cr\cr


\strong{GVF.db$assign()}\cr
Function \code{GVF.db$assign} has just a single mandatory argument: \code{value}. The object passed to argument \code{value} can only be a previously exported copy of \code{GVF.db}, i.e. an object of class \code{GVF.db_exported}. The function overwrites the \emph{current} version of \code{GVF.db} with \code{value}. As a result, after invoking \code{GVF.db$assign}, the content of \code{GVF.db} is \code{value}.\cr\cr


\strong{GVF.db$reset()}\cr
Function \code{GVF.db$reset} has no mandatory arguments and simply restores the default version of \code{GVF.db} (i.e. the one containing built-in GVF models only).\cr\cr
}

\source{
Built-in GVF models for frequencies (i.e. those with \code{Model.id} \code{1}, \code{2}, and \code{3}) are discussed in Chapter 7 of [Wolter 07], along with their theoretical justification. Built-in GVF models for totals (i.e. those with \code{Model.id} \code{4}, and \code{5}) lack a rigorous justification, but have sometimes been used successfully on a purely empirical basis. For instance, Istat surveys on structural business statistics adopted models of that kind to summarize standard errors in publications and to allow their approximate evaluation on a custom basis.
}

\references{
Wolter, K.M. (2007) \emph{\dQuote{Introduction to Variance Estimation}}, Second Edition, Springer-Verlag, New York.
}

\author{
Diego Zardetto
}

\seealso{
\code{\link{estimator.kind}} to assess what kind of estimates are stored inside a survey statistic object, \code{\link{gvf.input}} and \code{\link{svystat}} to prepare the input for GVF model fitting, \code{\link{fit.gvf}} to fit GVF models, \code{\link{plot.gvf.fit}} to get diagnostic plots for fitted GVF models, \code{\link{drop.gvf.points}} to drop alleged outliers from a fitted GVF model and simultaneously refit it, and \code{\link{predictCV}} to predict CV values via fitted GVF models.
}

\examples{
# Print the current content of GVF.db (invoking
# print(GVF.db) would do the same):
GVF.db

# Inspect the structure of the GVF.db data frame:
data.class(GVF.db) 
str(GVF.db)
dim(GVF.db)
nrow(GVF.db)


######################
# Accessor functions #
######################

# Delete the 3rd model:
GVF.db$delete(3)
# Print GVF.db (note that Model.id has been renumbered,
# so as to range always from 1 to nrow(GVF.db))
GVF.db

# Now delete the 1st model:
GVF.db$delete(1)
GVF.db

# Reset GVF.db to its default values:
GVF.db$reset()
GVF.db

# Insert a new tentative GVF model for Totals:
GVF.db$insert(CV ~ I(1/Y^2) + I(1/Y) + Y + I(Y^2), "Total", "resp")
GVF.db
# (notice that invoking GVF.db$insert() with first argument of type character,
# i.e. GVF.model="CV~I(1/Y^2)+I(1/Y)+Y+I(Y^2)", would have obtained exactly the
# same result)

# Now suppose you have somehow validated your newly added model,
# and you want to save your current, enhanced GVF.db in order to
# be able to use it later in a subsequent ReGenesees session.
### This can be achieved as follows:
### START
  # 1. You must first get a copy of it, by using accessor function
  #    GVF.db$get:
myGVF.db <- GVF.db$get()
myGVF.db
data.class(myGVF.db)

  # 2. Then, you must save the copy to a .RData workspace, in order
  #    to be able to load it later when needed, e.g.:
\dontrun{
save(myGVF.db, file="custom.GVF.Archive.RData")
}

  # 3. Starting a new ReGenesees session will set the default GVF.db,
  #    which we can simulate in this example as follows:
GVF.db$reset()
GVF.db

  # 4. Now you can load your previously saved customized GVF.db...
\dontrun{
load("custom.GVF.Archive.RData")
}
  # ...so that myGVF.db is back into your .GlobalEnv:
myGVF.db

  # 5. Lastly, you must overwrite GVF.db with your custom
  #    GVF archive myGVF.db via function GVF.db$assign:
GVF.db$assign(myGVF.db)
GVF.db

### Now your custom GVF archive is ready to be used by ReGenesees.
### STOP

# Illustrate some GVF.db$insert checks by trying crazy models
# or ill-specified attributes

# Examples start: reset GVF.db to its default values
GVF.db$reset()
GVF.db

# GVF model must be "syntactically new"...
\dontrun{
GVF.db$insert(log(CV^2) ~ log(Y))
}
# ...if this is the case, it can even be "equivalent" to old ones: e.g.
# the following is identical to model number 5 and will produces identical 
# estimates and predictions (as you may want to check):
GVF.db$insert(I(sqrt(VAR)/Y) ~ I(1/Y) + Y, "Total", Resp.to.CV = "resp")
GVF.db

# GVF model must have a response term
\dontrun{
GVF.db$insert(~ log(Y))
}

# GVF model response must involve some of 'SE', 'CV', 'VAR'
\dontrun{
GVF.db$insert(DEFF ~ log(Y))
}

# GVF model predictor must involve 'Y'
\dontrun{
GVF.db$insert(VAR ~ SE)
}

# If passed, Resp.to.CV can only involve 'resp' and 'Y'
\dontrun{
GVF.db$insert(I(sqrt(VAR)/Y) ~ I(1/Y) + Y + I(Y^2), Resp.to.CV = "sqrt(VAR)/Y")
}

# Examples end: reset GVF.db to its default values:
GVF.db$reset()
}

\keyword{survey}
