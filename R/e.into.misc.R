`%into%` <- function(inner, outer){
#####################################################
# Some tricks to reduce the sparseness of the known #
# totals dataframe in calibration, when auxiliary   #
# information are available for nested factors      #
# levels (e.g. for provinces and for regions).      #
#####################################################

# Some checks on supplied arguments
if (!is.factor(inner))
    stop("'inner' variable is not a factor")
if (any(is.na(inner)))
    stop("Missing values in 'inner' factor")

if (!is.factor(outer))
    stop("'outer' variable is not a factor")
if (any(is.na(inner)))
    stop("Missing values in 'outer' factor")

if (!identical(len <- length(inner), length(outer)))
    stop("Length mismatch between 'inner' and 'outer' factors")

# Drop unused levels, if any:
inner <- factor(inner)
outer <- factor(outer)

# Check that levels of factor 'inner' are actually 
# strictly nested inside those of factor 'outer'.            
    # Following conversions for speeding up tapply
    ext.unit <- as.numeric(outer)
    tab.o.in.i <- tapply(ext.unit, inner, function(x) length(unique(x)) )
    int.unit <- as.numeric(inner)
    tab.i.in.o <- tapply(int.unit, outer, function(x) length(unique(x)) )
    
    if (any(tab.o.in.i > 1) || all(tab.i.in.o==1))
        stop("'inner' isn't actually nested inside 'outer'!")

# end of checks.

# Start actual processing
# Build an empty output vector to be filled
res <- rep(NA, len)

# Build a key
id <- 1:len

# Build an order index for inner levels inside outer levels
# (an integer vector, see unlist)

in.seq <- unlist( tapply(int.unit, outer, function(x) as.integer(factor(x)) ) )

# Build the id sequence inside outer levels
# (an integer vector again, see unlist)
id.seq <- unlist( tapply(id, outer, function(x) x) )

# Now copy the order index into the appropriate elements of res
# (recall that factors are not generally ordered by outer)
res[id.seq] <- in.seq

# Last turn the latter into a factor (whose levels will be
# as much as the maximum number of inner inside a outer)
res <- factor(res)
res
}
