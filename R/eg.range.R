`g.range` <-
function (cal.design)
#############################################################################
#  Dato un oggetto di classe cal.analytic, calcola il range degli g-weights #
#############################################################################
{
    if (!inherits(cal.design, "cal.analytic")) 
        stop("Object 'cal.design' must be of class cal.analytic")

    e.df <- cal.design$variables
    weights <- attr(cal.design, "weights")    
    w.cal.char <- all.vars(weights)
    w.char <- substr(w.cal.char, 0, nchar(w.cal.char) - 4)
    ww <- e.df[, w.char]
    ww.cal <- e.df[, w.cal.char]
    g <- ww.cal/ww
    # check: are there NaN g values?
    g.nan <- g[is.nan(g)]
    if (length(g.nan)==0) {
         # No NaN g values found
         grange <- range(g)
        }
    else{
         # NaN g values found
         grange <- range(g[!is.nan(g)])
         # check: are there NaN not arising from 0/0?
         if ( !all( ww[is.nan(g)]==ww.cal[is.nan(g)] & ww[is.nan(g)]==0 ) ) {
            warning("NaN g-weights not arising from 0/0!")
            }
         else {
               # Treat 0/0 as 1
               grange <- range(grange, 1)
              }
        }

    names(grange) <- c("g.min", "g.max")
    grange
}
