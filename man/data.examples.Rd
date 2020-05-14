\name{data.examples}
\alias{data.examples}
\alias{example}
\alias{pop01}
\alias{pop02}
\alias{pop03}
\alias{pop03p}
\alias{pop04}
\alias{pop04p}
\alias{pop05}
\alias{pop05p}
\alias{pop06p}
\alias{pop07}
\alias{pop07p}
\alias{pop07pp}
\alias{bounds}
\docType{data}
\title{Artificial Household Survey Data}
\description{
Example data frames. Allow to run \R code contained in the \sQuote{Examples} section of the \pkg{ReGenesees} package help pages.
}
\usage{data(data.examples)}
\format{
The main data frame, named \code{example}, contains (artificial) data from a two stage stratified cluster sampling design. The sample is made up of 3000 final units, for which the following 21 variables were observed:
  \describe{
    \item{\code{towcod}}{Code identifying "variance PSUs": towns (PSUs) in not-self-representing (NSR) strata, families (SSUs) in self-representing (SR) strata, \code{numeric}}
    \item{\code{famcod}}{Code identifying families (SSUs), \code{numeric}}
    \item{\code{key}}{Key identifying final units (individuals), \code{numeric}}
    \item{\code{weight}}{Initial weights, \code{numeric}}
    \item{\code{stratum}}{Stratification variable, \code{factor} with levels \code{801} \code{802} \code{803} \code{901} \code{902} \code{903} \code{904} \code{905} \code{906} \code{907} \code{908} \code{1001} \code{1002} \code{1003} \code{1004} \code{1005} \code{1006} \code{1007} \code{1008} \code{1009} \code{1101} \code{1102} \code{1103} \code{1104} \code{3001} \code{3002} \code{3003} \code{3004} \code{3005} \code{3006} \code{3007} \code{3008} \code{3009} \code{3010} \code{3011} \code{3012} \code{3101} \code{3102} \code{3103} \code{3104} \code{3105} \code{3106} \code{3107} \code{3108} \code{3201} \code{3202} \code{3203} \code{3204} \code{5401} \code{5402} \code{5403} \code{5404} \code{5405} \code{5406} \code{5407} \code{5408} \code{5409} \code{5410} \code{5411} \code{5412} \code{5413} \code{5414} \code{5415} \code{5416} \code{5501} \code{5502} \code{5503} \code{5504} \code{9301} \code{9302} \code{9303} \code{9304} \code{9305} \code{9306} \code{9307} \code{9308} \code{9309} \code{9310} \code{9311} \code{9312}}
    \item{\code{SUPERSTRATUM}}{Collapsed strata variable (eliminates lonely PSUs), \code{factor} with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22} \code{23} \code{24} \code{25} \code{26} \code{27} \code{28} \code{29} \code{30} \code{31} \code{32} \code{33} \code{34} \code{35} \code{36} \code{37} \code{38} \code{39} \code{40} \code{41} \code{42} \code{43} \code{44} \code{45} \code{46} \code{47} \code{48} \code{49} \code{50} \code{51} \code{52} \code{53} \code{54} \code{55}}
    \item{\code{sr}}{Strata type, \code{integer} with values \code{0} (NSR strata) and \code{1} (SR strata)}
    \item{\code{regcod}}{Code identifying regions, \code{factor} with levels \code{6} \code{7} \code{10}}
    \item{\code{procod}}{Code identifying provinces, \code{factor} with levels \code{8} \code{9} \code{10} \code{11} \code{30} \code{31} \code{32} \code{54} \code{55} \code{93}}
    \item{\code{x1}}{Indicator variable (integer), \code{numeric}}
    \item{\code{x2}}{Indicator variable (integer), \code{numeric}}
    \item{\code{x3}}{Indicator variable (integer), \code{numeric}}
    \item{\code{y1}}{Indicator variable (integer), \code{numeric}}
    \item{\code{y2}}{Indicator variable (integer), \code{numeric}}
    \item{\code{y3}}{Indicator variable (integer), \code{numeric}}
    \item{\code{age5c}}{Age variable with 5 classes, \code{factor} with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
    \item{\code{age10c}}{Age variable with 10 classes, \code{factor} with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10}}
    \item{\code{sex}}{Sex variable, \code{factor} with levels \code{f} \code{m}}
    \item{\code{marstat}}{Marital status variable, \code{factor} with levels \code{married} \code{unmarried} \code{widowed}}
    \item{\code{z}}{A continuous quantitative variable, \code{numeric}}
    \item{\code{income}}{Income variable, \code{numeric}}
  }
}
\details{
Objects \code{pop01}, \ldots, \code{pop07pp} contain known population totals for various calibration models. Object pairs with names differing in the '\code{p}' suffix (such as \code{pop03} and \code{pop03p}) refer to the \emph{same} calibration problem but pertain to \emph{different} solution methods (global and partitioned respectively, see \code{\link{e.calibrate}}). The two-component numeric vector \code{bounds} expresses a possible choice for the allowed range for the ratios between calibrated weights and direct weights in the aforementioned calibration problems.
}
\section{Warning}{
\strong{Data in the \code{example} data frame are artificial.} The \emph{structure} of \code{example} intentionally resembles the one of typical household survey data, but the \emph{values} it stores are unreliable. The only purpose of such data is that they can be fruitfully exploited to illustrate the syntax and the working mechanism of the functions provided by the \pkg{ReGenesees} package.
}
\examples{
data(data.examples)
head(example)
str(example)
}
\keyword{datasets}