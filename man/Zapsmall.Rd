\name{Zapsmall}
\alias{Zapsmall}
\alias{Zapsmall.default}
\alias{Zapsmall.data.frame}
\title{Zapsmall Data Frame Columns and Numeric Vectors}
\description{
Put to zero values "close" to zero.
}
\usage{
\method{Zapsmall}{default}(x, digits = getOption("digits"), ...)
\method{Zapsmall}{data.frame}(x, digits = getOption("digits"), except = NULL, ...)
}
\arguments{
  \item{x}{A \code{data.frame} with numeric columns or a numeric vector.}
  \item{digits}{Integer indicating the precision to be used.}
  \item{except}{Indices of columns not to be zapped (if any).}
  \item{\dots}{Arguments for future expansion.}
}
\details{
This function "extends" to \code{data.frame} objects function \code{\link[base]{zapsmall}} from the package \pkg{base}.\cr
The method for class data.frame 'zaps' values close to zero occurring in columns of \code{x}. Argument \code{except} can be used to prevent specific columns from being zapped.\cr
The default method is a bare copy of the original function from package \pkg{base}.
}
\value{
An object of the same class of \code{x}, with values "close" to zero zapped to zero.
}
\author{
Diego Zardetto
}
\seealso{
The original function \code{\link[base]{zapsmall}} from package \pkg{base}.
}
\examples{
# Create a test data frame with columns containing
# values of different orders of magnitude:
data <- data.frame(a = pi*10^(-8:1), b = c(rep(1000,8), c(1E-5, 1E-6)))

# Print on screen the test data frame:
data

# Compare with its zapped version:
Zapsmall(data)
}
\keyword{survey}