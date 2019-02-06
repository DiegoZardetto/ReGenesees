\name{sbs}
\alias{sbs}
\alias{sbs.frame}
\docType{data}

\title{
Artificial Structural Business Statistics Data
}

\description{
The \code{sbs} data frame stores artificial sbs-like sampling data, while \code{sbs.frame} is the artificial sampling frame from which the \code{sbs} units have been drawn. They allow to run \R code contained in the \sQuote{Examples} section of the \pkg{ReGenesees} package help pages.
}

\usage{data(sbs)}
\format{
The \code{sbs} data frame mimics data observed in a Structural Business Statistics survey, under a one-stage stratified unit sampling design. The sample is made up of 6909 units, for which the following 22 variables were observed:
\describe{
    \item{\code{id}}{Identifier of the sampling units (enterprises), \code{numeric}}
    \item{\code{public}}{Does the enterprise belong to the Public Sector? \code{factor} with levels \code{0} (No) and \code{1} (Yes)}
    \item{\code{emp.num}}{Number of employees, \code{numeric}}
    \item{\code{emp.cl}}{Number of employees classified into 5 categories, \code{factor} with levels \code{[6,9]} \code{(9,19]} \code{(19,49]} \code{(49,99]} \code{(99,Inf]} (notice that small enterprises with less than 6 employees fell outside the scope of the survey)}
    \item{\code{nace5}}{Economic Activity code with 5 digits, \code{factor} with \code{596} levels} 
    \item{\code{nace2}}{Economic Activity code with 2 digits, \code{factor} with \code{57} levels}
    \item{\code{area}}{Territorial Division, \code{factor} with \code{24} levels}
    \item{\code{cens}}{Flag identifying statistical units to be censused (hence defining take-all strata), \code{factor} with levels \code{0} (No) and \code{1} (Yes)}
    \item{\code{region}}{Macroregion, \code{factor} with levels \code{North} \code{Center} \code{South}}
    \item{\code{va.cl}}{Class of Value Added, \code{factor} with \code{27} levels}
    \item{\code{va}}{Value Added, \code{numeric} (contains \code{NA}s)}
    \item{\code{dom1}}{A planned estimation domain, \code{factor} with \code{261} levels (\code{dom1} crosses \code{nace2} and \code{emp.cl})}
    \item{\code{nace.macro}}{Economic Activity Macrosector, \code{factor} with levels \code{Agriculture} \code{Industry} \code{Commerce} \code{Services}}
    \item{\code{dom2}}{A planned estimation domain, \code{factor} with \code{12} levels (\code{dom2} crosses \code{nace.macro} and \code{region})}
    \item{\code{strata}}{Stratification Variable, a \code{factor} with \code{664} levels (obtained by crossing variables \code{region}, \code{nace2}, \code{emp.cl} and \code{cens})} 
    \item{\code{va.imp1}}{Value Added Imputed1, \code{numeric} (\code{NA}s were replaced with average values computed inside imputation strata obtained by crossing \code{region}, \code{nace.macro}, \code{emp.cl})}
    \item{\code{va.imp2}}{Value Added Imputed2, \code{numeric} (\code{NA}s were replaced with median values computed inside imputation strata obtained by crossing \code{region}, \code{nace.macro}, \code{emp.cl})}
    \item{\code{y}}{A \code{numeric} variable correlated with \code{va}}
    \item{\code{weight}}{Direct weights, \code{numeric}}
    \item{\code{fpc}}{Finite Population Corrections (given as sampling fractions inside strata), \code{numeric}}
    \item{\code{ent}}{Convenience \code{numeric} variable identically equal to \code{1} (sometimes useful, e.g. to estimate the total number of enterprises)}
    \item{\code{dom3}}{An unplanned estimation domain, \code{factor} with \code{4} levels}
  }

The \code{sbs.frame} sampling frame (from which \code{sbs} units have been drawn) contains 17318 units.
}

\examples{
data(sbs)
head(sbs)
str(sbs)
str(sbs.frame)
}
\keyword{datasets}