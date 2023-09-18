\name{ReGenesees-package}
\alias{ReGenesees-package}
\alias{ReGenesees}
\docType{package}
\title{
ReGenesees: a Package for Design-Based and Model-Assisted Analysis of Complex
Sample Surveys
}

\description{
\pkg{ReGenesees} is an \R package for design-based and model-assisted analysis
of complex sample surveys. It handles multistage, stratified, clustered,
unequally weighted survey designs. Sampling variance estimation for nonlinear
(smooth) estimators is obtained by Taylor-series linearization. Sampling
variance estimation for multistage designs can be obtained both under the
Ultimate Cluster approximation or by means of an actual multistage computation.
\pkg{ReGenesees} offers comprehensive and advanced functionalities for
calibration of survey weights. Estimates, standard errors, confidence intervals
and design effects are provided for Horvitz-Thompson and Calibration estimators
of: Totals, Means, absolute and relative Frequency Distributions (marginal,
conditional and joint), Ratios, Shares and Ratios of Shares, Multiple Regression
Coefficients and Quantiles (variance via the Woodruff method). \pkg{ReGenesees}
also handles Complex Estimators, i.e. any user-defined estimator that can be
expressed as an analytic function of Horvitz-Thompson or Calibration estimators
of Totals or Means, by automatically linearizing them. The Design Covariance and
Correlation between Complex Estimators is also provided. Moreover, the package
can compute estimates and sampling errors of complex Measures of Change derived
from two not necessarily independent samples. All the analyses above can be
carried out for arbitrary subpopulations. In addition, \pkg{ReGenesees} can trim
calibration weights while preserving all the calibration constraints, and
perform \sQuote{special purpose calibration} tasks, i.e. can calibrate on
complex population parameters like Multiple Regression Coefficients.
\pkg{ReGenesees} also offers a Generalized Variance Functions (GVF)
infrastructure, i.e. facilities for defining, fitting, testing and plotting GVF
models, and to exploit them to predict variance estimates. Lastly, the package
offers simple survey planning tools to estimate sample size requirements and
perform power calculations.

The \pkg{ReGenesees} \emph{package} is the fundamental building block of a
full-fledged \R-based software \emph{system}: the \strong{ReGenesees System}.
The latter has a clear-cut two-layer architecture. The application layer of the
system is embedded into package \pkg{ReGenesees}. A second \R package, called
\pkg{ReGenesees.GUI}, implements the presentation layer of the system,
namely a user-friendly Tcl/Tk GUI.
}

\section{A Quick Reading Guide to the Reference Manual}{
This reference manual reports a documentation entry for each (user visible)
function of package \pkg{ReGenesees}. As you may have noticed by reading section
\sQuote{\R topics documented} (page 1 of the pdf manual), these documentation
entries are automatically sorted according to the alphabetic ordering of the
names of the functions. Such an ordering doesn't provide any clue about where
should a user start reading, nor on the best way to proceed further.

In section \sQuote{Table of Contents}, I tried to cluster the most important
topics documented in the reference manual into few broad groups, based on both
the statistical goals and on the software design of the underlying functions.

Moreover, I provided a \emph{relevance code} for each documented topic/function.
The meaning of such codes, along with the corresponding \emph{reading
suggestions}, are reported in the following table:\cr

\strong{Relevance Codes Legend}
\preformatted{
CODE    RELEVANCE           READING SUGGESTION
 ***    Very Important......Read these topics as soon as possible. A clear
                            understanding of these functions is mandatory
                            in order to start using profitably the package.

  **    Important...........Read these topics once you have been experiencing
                            for a while with (at least some of) the 'Very
                            Important' functions.

   *    Useful..............These functions are ancillary (albeit in
                            different ways) to the 'Very Important' and
                            'Important' ones (and their usage is generally
                            simpler).

   .    Advanced............These topics are very relevant but, unfortunately
                            quite difficult. As they involve technical
                            details, you should postpone their reading until
                            you become familiar with the package.
}

\strong{Important Notice}\cr
It goes without saying that the \sQuote{\strong{Examples}} sections at the end
of each documented topic \strong{represent a crucial part of this reference
manual}.\cr\cr
}

\section{TABLE OF CONTENTS}{

\subsection{Survey Design}{
\preformatted{
***  e.svydesign..........Specification of a Complex Survey Design
  *  weights..............Retrieve Sampling Units Weights
  *  find.lon.strata......Find Strata with Lonely PSUs
 **  collapse.strata......Collapse Strata Technique for Eliminating
                          Lonely PSUs
  *  des.addvars..........Add Variables to Design Objects
  *  des.merge............Merge New Survey Data into Design Objects
 **  smooth.strat.jump....Smooth Weights to Cope with Stratum Jumpers
}
}

\subsection{Calibration}{
\preformatted{
 **  pop.template.........Template Data Frame for Known Population Totals
  *  population.check.....Compliance Test for Known Totals Data Frames
  *  pop.desc.............Natural Language Description of Known Totals
                          Templates
 **  fill.template........Fill the Known Totals Template for a
                          Calibration Task
  *  pop.plot.............Plot Calibration Control Totals vs Current
                          Estimates
  *  bounds.hint..........A Hint for Range Restricted Calibration
***  e.calibrate..........Calibration of Survey Weights
  *  check.cal............Calibration Convergence Check
 **  trimcal..............Trim Calibration Weights while Preserving
                          Calibration Constraints
  *  g.range..............Range of g-Weights
  .  get.residuals........Calibration Residuals of Interest Variables
  .  get.linvar...........Linearized Variable(s) of Complex Estimators
                          by Domains
  *  ext.calibrated.......Make ReGenesees Digest Externally Calibrated
                          Weights
  .  contrasts.RG.........Set, Reset or Switch Off Contrasts for
                          Calibration Models
  .  \%into\%...............Compress Nested Factors
}
}

\subsection{Special Purpose Calibration}{
\preformatted{
  .  prep.calBeta.........Prepare a Survey Design to Calibration on
                          Multiple Regression Coefficients
  .  pop.calBeta..........Prepare Control Totals for Calibration on
                          Multiple Regression Coefficients
  .  pop.fuse.............Fuse Control Totals Data Frames for Special
                          Purpose and Ordinary Calibration Tasks
}
}

\subsection{Estimates and Sampling Errors}{
\preformatted{
***  svystatTM............Estimation of Totals and Means in
                          Subpopulations
***  svystatR.............Estimation of Ratios in Subpopulations
***  svystatS.............Estimation of Shares in Subpopulations
***  svystatSR............Estimation of Share Ratios in Subpopulations
***  svystatB.............Estimation of Population Regression Coefficients in
                          Subpopulations
***  svystatQ.............Estimation of Quantiles in Subpopulations
***  svystatL.............Estimation of Complex Estimators in
                          Subpopulations
***  svySigma.............Estimation of the Population Standard Deviation of
                          a Variable
***  svySigma2............Estimation of the Population Variance of a Variable
***  svyDelta.............Estimation of a Measure of Change from Two
                          Not Necessarily Independent Samples
  *  details..............Details on svyDelta results
 **  aux.estimates........Quick Estimates of Auxiliary Variables Totals
 **  CoV, Corr............Design Covariance and Correlation of Complex
                          Estimators in Subpopulations
  *  write.svystat........Export Survey Statistics
  *  extractors...........Extractor Functions for Variability Statistics
  .  ReGenesees.options...Variance Estimation Options for the ReGenesees
                          Package
}
}

\subsection{Generalized Variance Functions Method}{
\preformatted{
***  GVF.db...............Archive of Registered GVF Models
***  gvf.input............Prepare Input Data to Fit GVF Models
***  svystat..............Compute Many Estimates and Errors in Just a
                          Single Shot
***  fit.gvf..............Fit GVF Models
 **  plot.gvf.fit.........Diagnostic Plots for Fitted GVF Models
 **  drop.gvf.points......Drop Outliers and Refit a GVF Model
  *  getR2, AIC, BIC......Quality Measures on Fitted GVF Models
  *  getBest..............Identify the Best Fit GVF Model
***  predictCV............Predict CV Values via Fitted GVF Models
  *  gvf.misc.............Miscellanea: Methods for Fitted GVF Models 
  *  estimator.kind.......Which Estimator Did Generate these
                          Survey Statistics?
}
}

\subsection{Sample Size and Power}{
\preformatted{
 **  n.prop................Sample Size Requirements for the Estimation
                           of a Proportion
 **  prec.prop.............Expected Precision Level in the Estimation
                           of a Proportion
 **  n.comp2prop...........Power Calculations for a Test that Compares
                           Two Estimated Proportions: Sample Size
 **  pow.comp2prop.........Power Calculations for a Test that Compares
                           Two Estimated Proportions: Expected Power
 **  mde.comp2prop.........Power Calculations for a Test that Compares
                           Two Estimated Proportions: Expected Minimum
                           Detectable Effect
 **  n.mean................Sample Size Requirements for the Estimation
                           of a Mean
 **  prec.mean.............Expected Precision Level in the Estimation
                           of a Mean
 **  n.comp2mean...........Power Calculations for a Test that Compares
                           Two Estimated Means: Sample Size
 **  pow.comp2mean.........Power Calculations for a Test that Compares
                           Two Estimated Means: Expected Power
 **  mde.comp2mean.........Power Calculations for a Test that Compares
                           Two Estimated Means: Expected Minimum
                           Detectable Effect						   
}
}

\subsection{Diagnostics and Utilities}{
\preformatted{
  *  UWE...................Unequal Weighting Effect
  *  Zapsmall..............Zapsmall Data Frame Columns and Numeric Vectors
}
}

\subsection{Data Sets}{
\preformatted{
 **  data.examples........Artificial Household Survey Data
 **  fpcdat...............A Small But Not Trivial Artificial Sample
                          Data Set
 **  sbs..................Artificial Structural Business Statistics Data
 **  Delta.el.............Two Artificial Samples of Elementary Units for
                          Estimation of Change
 **  Delta.clus...........Two Artificial Cluster Samples for Estimation
                          of Change
 **  AF.gvf...............Example Data for GVF Model Fitting
}
}

The ordering of the above \sQuote{Table of Contents} reflects only loosely the
procedural sequence in which functions could be used. For instance, while you
cannot apply function \code{\link{e.calibrate}} unless you have previously built
a design object by using \code{\link{e.svydesign}}, you can exploit, e.g.,
function \code{\link{collapse.strata}} also after calibration. As a further
example, all functions in group \sQuote{Estimates and Sampling Errors} can be
used on objects created by \code{\link{e.svydesign}} (yielding estimates and
sampling errors for functions of Horvitz-Thompson estimators), as well as on
objects created by \code{\link{e.calibrate}} (yielding estimates and sampling
errors for functions of Calibration estimators).
}

\references{
Zardetto, D. (2015) \emph{\dQuote{ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys}}. Journal of Official Statistics, 31(2), 177-203. doi: https://doi.org/10.1515/jos-2015-0013.
}

\keyword{package}