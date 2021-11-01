\name{smooth.strat.jump}
\alias{smooth.strat.jump}

\title{Smooth Weights to Cope with Stratum Jumpers}

\description{
Given a stratified one-stage unit sampling design object, this function smooths survey weights to mitigate estimation issues that may arise from stratum jumpers.
}

\usage{
smooth.strat.jump(design, curr.strata, method = c("MinChange", "Beaumont"))
}

\arguments{
  \item{design}{Object of class \code{analytic} containing the weights to be smoothed. It must be a one-stage unit sampling design object. Moreover, it must be a non-calibrated object. See \sQuote{Details}.}
  \item{curr.strata}{Formula identifying the \emph{current} strata variable, as observed at survey-time (see \sQuote{Details}).}
  \item{method}{The smoothing method (see \sQuote{Details}). The default method \code{'MinChange'} smooths the weight of stratum jumpers almost without affecting the weights of other units. Method \code{'Beaumont'} adopts a much more aggressive smoothing strategy, which may result in significant modifications of the weights also for units that are \emph{not} stratum jumpers (see \sQuote{Examples}).}
}

\details{
In business surveys, stratum jumpers are sampling units (e.g. firms or establishment) whose stratum information observed at survey-time happen to differ from the stratum information that was available in the sampling frame at design-time.

If the current value of the strata variable (i.e. the one observed at survey-time) is reliable, stratum jumpers are evidence of frame imperfections (typically, the frame was not up-to-date). Empirically, stratum jumpers are often units that underwent a fast growth in size from sampling-time to survey-time.

In most enterprise surveys the sampling design is such that smaller firms receive a smaller inclusion probability (and hence a larger design weight). Therefore, stratum jumpers often have a \dQuote{too large} weight, in the sense that they would have received a smaller weight had their actual size been known at sampling-time. When these units also happen to have a large value of interest variable \code{y}, they may become influential in estimation. Even though stratum jumpers may exist in household surveys too, business surveys are much more exposed to the risk that stratum jumpers unduly influence estimation, as their target populations are typically highly skewed with respect to many interest variables. As a consequence, in business surveys, stratum jumpers can result in inefficient (and, under some circumstances, even biased) design-based estimators.

Despite stratum jumpers are actually a concern only when both their design weight \code{w} and their \code{y} value (as measured at survey-time) conspire to yield an influential value of \code{w * y}, function \code{smooth.strat.jump} tries to mitigate their potential adverse impact by smoothing the weights \emph{without} using any information on \code{y}. Such a choiche is driven by two major aims. First, the methodology must preserve the universality of the weights (i.e. the same weights must be used to compute estimates for whatever interest variables \code{y}). Second, the methodology should not require explicit modeling efforts and be easy enough to scale to production settings that need automated and replicable procedures. 

Argument \code{design} identifies the survey design object that is (possibly) affected by stratum jumpers and contains the weights that will be smoothed. In case no stratum jumpers are found, the function will raise an error. Object \code{design} can only be a one-stage unit sampling design. Moreover, it must be a non-calibrated object. Should any of these conditions be false, the function would raise an error.

Formula \code{curr.strata} defines the \emph{current} strata variable, as observed at survey-time. This is different from the \emph{design} strata variable used to build object \code{design} using function \code{\link{e.svydesign}}, as the latter was measured at sampling-time. Function \code{smooth.strat.jump} will flag as stratum jumpers all the units whose \emph{current} stratum differs from the \emph{design} one. Note that, in case \emph{current} strata become available only after object \code{design} was created, you may use function \code{\link{des.addvars}} to add this new column to the old object.

The weight smoothing process entails two steps:
\enumerate{
  \item {The weights are smoothed according to a given method (see below): \code{w -> w1}}
  \item {The weights \code{w1} of all units are scaled by a global factor so as to preserve the initial overall sum of weights: \code{w1 -> w2 = scale * w1} with  \code{scale = sum(w)/sum(w1)}}
}

Argument \code{method} controls the smoothing algorithm. Two methods can be selected: \code{'MinChange'} (the default) and \code{'Beaumont'} (which implements the proposal of [Beaumont, Rivest 09]). Note that the methods only differ with respect to step 1., as step 2. is identical for both of them.

The step 1 working mechanism of these methods can be summarized as follows:\cr
\itemize{
 \item{\strong{\sQuote{MinChange}}}
}
Only the weights of stratum jumpers are smoothed, by setting their value to the average weight of units belonging to the same \emph{current} stratum. Weights of all other units are left unchanged.\cr
\itemize{
 \item{\strong{\sQuote{Beaumont}}}
}
The weights of all units are smoothed, by setting their value to the average weight of units belonging to the same \emph{current} stratum, with the exception of \emph{minimum weight} units, whose weights are left unchanged. Therefore all weights, exluding only minimum weights within each current stratum, are smoothed. Note that this smoothing affects all current strata, even those that do \emph{not} include any stratum jumper.\cr\cr

In summary, both methods often lead to very similar smoothed weights for units that are stratum jumpers. However method \code{'Beaumont'} smooths the weights of all other (i.e. \emph{non} stratum jumpers) units much more aggressively than method \code{'MinChange'} (which only minimally alters them in step 2. to preserve the overall sum of weights). Moreover, method \code{'MinChange'} treats stratum jumpers that grew in size and those that decreased in size on the same footing, whereas method \code{'Beaumont'} typically smooths the weights of the former more than those of the latter (owing to its minimum weight preservation constraint).

Note that every call to \code{smooth.strat.jump} generates, by side effect, a diagnostics data structure named \code{strat.jump.status} into the \code{.GlobalEnv} (see \sQuote{Examples}). This is a data frame with one row for each stratum jumper unit, whith the following columns:
\preformatted{
  Column                 Meaning
  IDS....................Unit identifier
  W......................Initial weight
  DES_STR................Design stratum
  DES_STR_W_AVG..........Average of initial weights within the design stratum
  CURR_STR...............Current stratum
  CURR_STR_W_AVG.........Average of initial weights within the current stratum
  N_JUMP_DES_CURR_STR....Number of stratum jumpers that jumped between the
                         design stratum and the current stratum (NOTE: in any
                         direction)
  W_SMOOTH_UNSC..........Unscaled smoothed weight (as obtained after step 1)
  W_SMOOTH...............Smoothed weight (scaled, as obtained after step 2)
}
}

\value{
An object of the same class as \code{design}. The data frame it contains (stored in its \code{$variables} slot) includes the smoothed weights columns and a column that flags the stratum jumpers. The name of the smoothed weights column is obtained by pasting the name of the initial weights column with the string \code{".smooth"}. Stratum jumpers are identified by a new (\code{logical}) column named \code{is.jumper}.  
}

\section{Methodological Remark}{
Smoothing survey weights is a model-based approach, see [Beaumont 08] (e.g.\code{method = 'Beaumont'} basically models the smoothed weights as a function of the current strata using a one-way ANOVA model plus constraints). Therefore Horvitz-Thompson-like estimators that use smoothed weights - instead of design weights - cannot be guaranteed to be design-unbiased. Of course, the need to smooth the weights arises precisely because the existence of stratum jumpers already signals a departure from the ideal conditions of design-based inference.

As the design-unbiasedness of Horvitz-Thompson estimators in probability sampling rests on using design weights that are reciprocals of inclusion probabilities, smoothing methods that change the design weights the least appear preferable in a design-based perspective. For this reason, function \code{smooth.strat.jump} adopts the \code{'MinChange'} method by default. One could, nonetheless, argue that the \code{'Beaumont'} method would sometimes perform better (e.g. lead to more efficient estimates) from a model-based perspective.

Regardless the choice of argument \code{method}, in order to reduce any possible design-bias introduced by smoothing the weights, users are advised to \emph{calibrate the smoothed weights} using any auxiliary information available from external sources that are \emph{more up-to-date} than the sampling frame.
}

\author{Diego Zardetto}

\references{
Beaumont, J. F. (2008). A new approach to weighting and inference in sample surveys. Biometrika, 95(3), 539-553.

Beaumont, J. F., Rivest, L. P. (2009). Dealing with outliers in survey data. In Handbook of statistics (Vol. 29, pp. 247-279). Elsevier.
}

\seealso{
\code{\link{e.svydesign}} to bind survey data and sampling design metadata and \code{\link{e.calibrate}} for calibrating smoothed survey weights by leveraging auxiliary information that is more up-to-date than the sampling frame (a warmly suggested option).
}

\examples{
##################################################################################
# Build [Beaumont, Rivest 09] example dataset, containing:                       #
# - a first stratum 'A' that, at survey time, contains one large weight stratum  #
#   jumper received from design-stratum 'B'                                      #
# - a second stratum 'B' that, at survey time, does not contain stratum jumpers  #
#                                                                                #
# and enhance it with:                                                           #
# - a third stratum 'C' that, at survey time, does not contain stratum jumpers   #
# - a fourth stratum 'D' that, at survey time, does not contain stratum jumpers  #
# - a fifth stratum 'E' that, at survey time, contains one medium weight stratum #
#   jumper received from design-stratum 'D'                                      #
##################################################################################
BR <- data.frame(
       id = 1:90,
       des.strata =  factor(rep(c("A", "B", "C", "D", "E"), c(9, 41, 10, 20, 10))),
       curr.strata = factor(rep(c("A", "B", "C", "D", "E"), c(10, 40, 10, 19, 11))),
       w = c(rep(c(1, 31), c(9, 41)), c(27, 28, 22, 26, 11, 12, 13, 30, 17, 21), 1:20, 10:1)
    )

# Have a look at the data:
BR

# Have a look at the jumps:
with(BR, table(des.strata, curr.strata))

# Use the BR data frame to build a one stage stratified unit sampling design: 
BRdes <- e.svydesign(data=BR, ids=~id, strata=~des.strata, weights=~w)


## Now smooth the weights:
## Method: MinChange (the default)
M.smooth <- smooth.strat.jump(BRdes, ~curr.strata)
M.smooth

# Have a look at the new columns:
head(M.smooth$variables)

# Inspect the effects of smoothing on the stratum jumpers:
strat.jump.status

## Method: Beaumont
B.smooth <- smooth.strat.jump(BRdes, ~curr.strata, method = "Beaumont")
B.smooth

# Inspect the effects of smoothing on the stratum jumpers:
strat.jump.status

## As anticipated, smoothed weights of stratum jumpers are almost identical for
## both methods. However the methods differ significantly when it comes to non
## stratum jumpers. This is clearly shown in the following plots.

## Plot 1 - START
opar <- par("mfcol" = c(1, 2))

# M.smooth
with(M.smooth$variables, plot(w, w.smooth, pch = 19, col = 1 + is.jumper,
     xlab = "Original Weights", ylab = "Smoothed Weights", main = "method: MinChange"))
abline(0:1, col = "limegreen", lwd = 2, lty = 2)

legend("topleft", 
legend = c("Stratum Jumper", "Non Stratum Jumper"), 
col = c("red", "black"), 
pch = c(19, 19),
bty = "n", 
text.col = "black", 
inset = c(0.05, 0.05)
)

# B.smooth
with(B.smooth$variables, plot(w, w.smooth, pch = 19, col = 1 + is.jumper,
     xlab = "Original Weights", ylab = "Smoothed Weights", main = "method: Beaumont"))
abline(0:1, col = "limegreen", lwd = 2, lty = 2)

legend("topleft", 
legend = c("Stratum Jumper", "Non Stratum Jumper"), 
col = c("red", "black"), 
pch = c(19, 19),
bty = "n", 
text.col = "black", 
inset = c(0.05, 0.05)
)

par(opar)
## Plot 1 - END

## Plot 2 - START
# M.smooth
with(M.smooth$variables, plot(w, pch = 20, col = curr.strata,
     ylab = "Original and Smoothed Weights",
     main = "method: MinChange \n(colors identify different current strata)"))
with(M.smooth$variables, points(w.smooth, pch = 0, col = curr.strata))

legend("topright", 
legend = c("Original Weight", "Smoothed Weight"), 
col = rep("grey", 2), 
pch = c(20, 0),
text.col = "black", 
inset = c(0.01, 0.01)
)

# B.smooth
with(B.smooth$variables, plot(w, pch = 20, col = curr.strata,
     ylab = "Original and Smoothed Weights",
     main = "method: Beaumont \n(colors identify different current strata)"))
with(B.smooth$variables, points(w.smooth, pch = 0, col = curr.strata))

legend("topright", 
legend = c("Original Weight", "Smoothed Weight"), 
col = rep("grey", 2), 
pch = c(20, 0),
text.col = "black", 
inset = c(0.01, 0.01)
)
## Plot 2 - END

}
\keyword{survey}