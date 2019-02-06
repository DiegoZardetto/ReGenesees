get.residuals <- function(cal.design, y, scale = c("no", "w", "d", "g"))
############################################################################
# This function computes (scaled) residuals of a set of interest variables #
# w.r.t. the calibration model adopted to build object 'cal.design'.       #
# 'scale' argument acts as follows:                                        #
# scale = "no" -> no scaling                  -> e                         #
# scale =  "w" -> calibration weights scaling -> w*e                       #
# scale =  "d" -> direct weights scaling      -> d*e                       #
# scale =  "g" -> g-weights scaling (g=w/d)   -> g*e                       #
############################################################################
{
  # First verify if the function has been called inside another function:
  # this is needed to correctly manage metadata when e.g. the caller is a
  # GUI stratum
  directly <- !( length(sys.calls()) > 1 )
  design.expr <- if (directly) substitute(cal.design)

  if (!inherits(cal.design, "cal.analytic")) 
      stop("Object 'cal.design' must be of class cal.analytic")
  if (!inherits(y, "formula")) 
      stop("Interest variables must be supplied as a formula")
  y.vars <- all.vars(y)
  if (length(y.vars) < 1)
      stop("Interest variables must belong to the survey dataframe")

  # Check for missing values in interest variables: CANNOT handle them!
  NA.estvars(design = cal.design, estvars = y.vars, draconian = TRUE)

  e.df <- cal.design$variables
  weights <- attr(cal.design, "weights")    
  w.cal.char <- all.vars(weights)
  # Get "previous step" weights

  # NOTE: Must behave differently for *trimmed* objects, i.e. previous step means
  #       *2 steps before*! 27/09/2016
  if (!is.trimmed(cal.design)) {
      w.char <- substr(w.cal.char, 0, nchar(w.cal.char) - 4)
    } else {
      w.char <- substr(w.cal.char, 0, nchar(w.cal.char) - 8)
    }

  ww <- e.df[, w.char]
  ww.cal <- e.df[, w.cal.char]
  g <- ww.cal/ww

  # Build y matrix (do the right thing with factors)
  mf <- model.frame(y, cal.design$variables, na.action=na.pass)
  yy <- lapply(attr(terms(y),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
  cols <- sapply(yy, NCOL)
  y <- matrix(nrow = NROW(yy[[1]]), ncol=sum(cols))
  scols <- c(0,cumsum(cols))
  for(i in 1:length(yy)){
    y[, scols[i]+1:cols[i]] <- yy[[i]]
  }
  colnames(y) <- do.call("c", lapply(yy,colnames))
  y <- as.matrix(y)

  if (is.character(scale)) {
      scale <- match.arg(scale)
    }

  postStrata <- cal.design$postStrata

  # Code below computes w-scaled residuals of y under the calibration model:
  # y -> y=w*e
  y <- y*ww.cal
  if(!is.null(postStrata)){
    for (psvar in postStrata){
        if (inherits(psvar, "analytic_calibration")){
        ####################################################
        # This is the code added to cope with cal.analytic #
        ####################################################
        nobs <- dim(cal.design$variables)[1]
        ## If svyrecvar has been called on a subset get
        ## observations index:
        domain.index <- attr(cal.design, "domain.index")
        ## else build it for the whole sample:
        if (is.null(domain.index)) domain.index <- 1:nobs
        if (psvar$stage!=0)
            stop("Analytic calibration must be at population level")
        ## residuals must have the same structure as y (thus a matrix)
        residuals <- y
        ## ...but filled by zeros
        residuals[, ] <- 0
        ## Now get the list of qr decompositions over partitions (if any)
        qr.list <- postStrata[[1]]$qr.list
          for (part in qr.list) {
               part.index <- part$group
               part.gwhalf <- part$gwhalf
               # Only partitions with non-empty overlap with the estimation domain matter
               if (any(domain.index %in% part.index)) {
                   # Check if any y value is NaN or infinite: this signals that the
                   # estimator gradient (recall y can actually be a linearized variable)
                   # is singular at the Taylor series expansion point.
                   # If this is the case, warn and avoid calling qr.resid (which
                   # would give an error, being Inf NA and NaN unmanageable by Fortran):
                   # return NaN "residuals" instead:
                   z.woodruff <- y[part.index, ,drop=FALSE]
                   if (any( is.infinite(z.woodruff) | is.nan(z.woodruff) )){
                       warning("Estimator gradient is singular at Taylor series expansion point!")
                       residuals[part.index, ] <- NaN
                       }
                   else {
                       residuals[part.index, ] <- qr.resid(part$qr,
                                                           z.woodruff/part.gwhalf) * part.gwhalf
                      }
                }
            }
        y <- residuals
        }
    }

 # Recall that y is w*e here: now scale it according to 'scale'
 e.scaled <- switch(scale, no = y/ww.cal, w = y, d = y/g, g = y/ww)
 attr(e.scaled, "design") <- design.expr
 attr(e.scaled, "y.vars") <- y.vars
 attr(e.scaled, "scale") <- scale
 return(e.scaled)
  }
  else {
     stop("This should not happen!")
    }
}
