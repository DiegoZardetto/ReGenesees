\name{PlotCI}
\alias{PlotCI}
\alias{BarPlotCI}

\title{
Visualize Domain Estimates and Confidence Intervals by Simple Plots and Bar Charts
}

\description{
These functions plot point estimates and confidence intervals by domains (subpopulations).
}

\usage{
PlotCI(stat, level = 0.95, eb.len = 0.05, eb.col = "black", eb.lwd = 1,
       xlab = "Domain", ylab = "Estimate", labels = NULL, ...)

BarPlotCI(stat, level = 0.95, eb.len = 0.05, eb.col = "black", eb.lwd = 1,
          xlab = "Domain", ylab = "Estimate", names.arg = NULL, pch = "", ...)
}

\arguments{
  \item{stat}{An object containing survey statistics.}
  \item{level}{Confidence level for the confidence intervals to be plotted as error bars: the default value is \code{0.95}.}
  \item{eb.len}{Half length of the error bars cap (in inches). Defaults to \code{0.05}.}
  \item{eb.col}{Color of the error bars. Defaults to \code{'black'}.}
  \item{eb.lwd}{Line width of the error bars. Defaults to \code{1}.}
  \item{xlab}{A label for the x axis. Defaults to \code{'Domain'}.}
  \item{ylab}{A label for the y axis. Defaults to \code{'Estimate'}.}
  \item{labels}{Labels to be placed at the tick points (see \code{\link[graphics]{axis}} for more details). Only relevant for function \code{\link{PlotCI}}.}
  \item{names.arg}{Vector of names to be plotted below each bar of the bar chart (see \code{\link[graphics]{barplot}} for more details). Only relevant for function \code{\link{BarPlotCI}}.}
  \item{pch}{An optional plotting symbol to mark point estimates. By default no symbol is plotted. Only relevant for function \code{\link{BarPlotCI}}.}
  \item{\dots}{Other parameters to be passed through to basic plotting functions (\code{\link[graphics]{plot}} for function \code{\link{PlotCI}}) and \code{\link[graphics]{barplot}} for function \code{\link{BarPlotCI}}).}
}

\details{
Given a survey statistic object \code{stat} containing survey estimates and estimated sampling errors - possibly by domains - these functions produce a simple plot (\code{\link{PlotCI}}) or a bar chart (\code{\link{BarPlotCI}}) depicting both point estimates and confidence intervals. Confidence intervals are depicted as error bars (eb).

Argument \code{stat} can be any object which has been returned by calling a survey statistics function (e.g. \code{\link{svystatTM}}, \code{\link{svystatR}}, \code{\link{svystatS}}, \code{\link{svystatSR}}, \code{\link{svystatB}}, \code{\link{svystatQ}}, \code{\link{svystatL}}, \code{\link{svySigma}}, and \code{\link{svySigma2}}).
}

\author{
Diego Zardetto
}

\seealso{
The underlying workhorse functions \code{\link[graphics]{plot}} and \code{\link[graphics]{barplot}}.
}

\examples{
##########
# PlotCI #
##########

# Load household data:
data(data.examples)

# Create a design object:
des <- e.svydesign(data = example, ids = ~towcod + famcod,strata = ~SUPERSTRATUM,
                   weights = ~weight)

# Estimate average income by province:
M_income_procod <- svystatTM(des, ~income, ~procod, estimator = "Mean")

# Take a look:
M_income_procod

# Visualize estimates and confidence intervals:
PlotCI(M_income_procod, ylab = "Average Income", xlab = "Province Code",
       pch = 19, labels = levels(des$variables$procod),
       main = "Average Income by Province")


# Estimate the total of variable y1 by age classes:
T_y1_age5c <- svystatTM(des, ~y1, ~age5c)

# Take a look:
T_y1_age5c

# Visualize estimates and confidence intervals. Since age class (the domain variable) is
# an ordered categorical variable, one may add a dashed tendency line:
PlotCI(T_y1_age5c, ylab = "Total y1", xlab = "Age Class", pch = 19,
       labels = levels(des$variables$age5c), main = "Total y1 by Age Class",
       ty = "o", lty = 2)


# Load sbs data:
data(sbs)

# Create a design object:
sbsdes <- e.svydesign(data = sbs, ids = ~id, strata = ~strata,
                      weights = ~weight, fpc = ~fpc)

# Estimate average value added per employee by region:
VA_per_employee_reg <- svystatR(sbsdes, num = ~va.imp2, den = ~emp.num, by = ~region)

# Take a look:
VA_per_employee_reg

# Compute population parameter from the business register (sbs.frame)
pop.par <- with(sbs.frame, tapply(va.imp2, region, sum) / tapply(emp.num, region, sum))

# Take a look:
pop.par

# Compare estimates and confidence intervals with population parameters:
PlotCI(VA_per_employee_reg, ylab = "Average Value Added per Employee",
       xlab = "Region", pch = 19, labels = levels(sbsdes$variables$region),
       main = "Average Value Added per Employee by Region")
points(pop.par, pch = 19, col = "red")
legend("topleft", 
  legend = c("Sample Estimate", "Population Parameter"), 
  col = c("black", "red"), 
  pch = rep(19, 2),
  text.col = "black", 
  inset = c(0.05, 0.05)
)

# NOTE: This approach can be used to, e.g., compare uncalibrated estimates with the
#       corresponding population totals before the calibration step.


#############
# BarPlotCI #
#############

# Build a bar plot to visualize estimates and confidence intervals of average value added
# per employee by region (computed above):
BarPlotCI(VA_per_employee_reg, ylab = "Average Value Added per Employee",
          xlab = "Region", pch = 19, col = adjustcolor(2:4, alpha.f = 0.3),
          names.arg = levels(sbsdes$variables$region),
          main = "Average Value Added per Employee by Region")


# Estimate population count by 10 age classes:
Nind_age10c <- svystatTM(des, ~age10c)

# Take a look:
Nind_age10c

# Visualize estimates and confidence intervals with a bar plot:
BarPlotCI(Nind_age10c, ylab = "Number of Individuals", xlab = "Age Class",
          names.arg = levels(des$variables$age10c), pch = 19,
          col = "lightgoldenrod1", eb.col = "blue", border = "blue",
          main = "Population by Age Class")

}

\keyword{survey}
