`pop.desc` <- function(pop.totals){
###########################################################################
# This function provides a natural language description of a known totals #
# data frame to be used for a calibration task.                           #
###########################################################################

if (data.class(pop.totals) != "pop.totals"){ 
     # I don't use inherits above to avoid objects of class 'aux.estimates',
     # which store estimated totals, rather than known ones (though their
     # structure would be perfectly described)
     stop("Input object must be of class pop.totals")
    }
if (!(length(pop.totals) > 0)) return(pop.totals)

# Get needed information from pop.totals attributes
pop.vars.type <- attr(pop.totals, "pop.vars.type")     # from pop.template 4 'pop.desc' #
var.by.term <- attr(pop.totals, "var.by.term")         # from pop.template 4 'pop.desc' #
cols <- attr(pop.totals, "calmodel.names")             # from pop.template 4 'pop.desc' #
which.term.col <- attr(pop.totals, "which.term.col")   # from pop.template 4 'pop.desc' #
partition <- attr(pop.totals, "partition")

# Is it a global calibration task?
is.global <- identical(partition, FALSE)
n.part.vars <- 0
aux.totals <- pop.totals
# If calibration task is partitioned, get partition variables
if (!is.global) {
     partition.vars <- all.vars(partition)
     n.part.vars <- length(partition.vars)
     # Drop partition columns to restrict to auxiliary info columns
     aux.totals <- pop.totals[, -(1:n.part.vars), drop = FALSE]
    }

# Does the template include an intercept term?
has.intercept <- any(which.term.col == 0)

# Identify blocks of population totals (namely those related to a single calmodel term)
n.aux <- ncol(aux.totals)
aux.blocks <- split(1:n.aux, which.term.col)
n.aux.blocks <- length(aux.blocks)

if (has.intercept) {
     names(aux.blocks)[names(aux.blocks) == 0] <- "Intercept"
    }

# Name each component of aux.blocks with the corresponding generating term
# Exclude model ~1, which as the intercept ONLY
if (has.intercept && n.aux.blocks == 1) {
     NULL
     }
else {
     names(aux.blocks)[names(aux.blocks) != "Intercept"] <- colnames(var.by.term[, as.integer(names(aux.blocks)[names(aux.blocks) != "Intercept"]), drop = FALSE])
    }

#-- Print on screen START --#
cat("# Data frame of known totals for a ", if (is.global) "*global* " else "*partitioned* ", "calibration task", "\n", sep = "")
cat("- Number of rows:         ", nrow(pop.totals), "\n", sep = "")
cat("- Number of columns:      ", ncol(pop.totals), "\n", sep = "")
cat("- Number of known totals: ", prod(dim(aux.totals)), "\n", sep = "")
cat("\n")
cat("# The data frame structure is as follows\n")

if (!is.global) {
     cat("## ", if (n.part.vars == 1) "Column 1 identifies " else paste("Columns 1-", n.part.vars, " identify ", sep = ""),
         nrow(pop.totals), " *calibration domains* (one for each row)", "\n", sep = "")
     cat("\n")
     print(pop.totals[, 1:n.part.vars, drop = FALSE])
     cat("\n")
    }

blocks.string <- paste( " organized into ", if (n.aux.blocks == 1) "*1 BLOCK* " else paste("*", n.aux.blocks, " BLOCKS* ", sep = ""), sep = "")
cat("## ", if (n.aux == 1) paste("Column ", n.part.vars + 1, " identifies known totals", sep = "") else paste("Columns ", n.part.vars + 1, "-", n.part.vars + n.aux, " identify known totals", sep = ""), blocks.string, "\n", sep = "")
cat("\n")
#-- Print on screen END --#

# Process population totals blocks
for (i in 1:n.aux.blocks){
     # Go on term by term...
     term.i <- names(aux.blocks)[i]
     if (term.i != "Intercept"){
         # Identify variables in the term...
         vars.i <- rownames(var.by.term)[var.by.term[, term.i] > 0]
         # ...and their types...
         types.i <- pop.vars.type[vars.i]
         names(types.i) <- vars.i
         # Normalize template column names arising from term:
         # strip factor names, leaving in place categories only
         factors.i <- vars.i[types.i == "categorical"]
         if (length(factors.i)) {
             for (f in factors.i){
                 cols[aux.blocks[[i]]] <- gsub(f, "", cols[aux.blocks[[i]]])
                }
            }
         attr(aux.blocks[[i]], "vars") <- types.i

        # Textual description of columns in block i
         rc.i <- unique(range(n.part.vars + aux.blocks[[i]]))
         c.i <- paste(if (length(rc.i) == 1) " - Column:             " else " - Columns:            ", paste(rc.i, collapse = "-", sep = " "))

        # Textual description of auxiliary information in block i
         num.i <- names(types.i)[types.i != "categorical"]
         cat.i <- names(types.i)[types.i == "categorical"]
         if (length(num.i) > 0) {
             if (length(cat.i) > 0) {
                 d.i <- paste(" - Benchmark:           Totals of", paste(num.i, collapse = "*"), "by", paste(cat.i, collapse = ", "), sep = " ")
               }
             else {
                 d.i <- paste(" - Benchmark:           Total of", paste(num.i, collapse = "*"), sep = " ")
               }
            }
         else {
             d.i <- paste(" - Benchmark:           Counts of", paste(cat.i, collapse = ":"), sep = " ")
            }

        # How many auxiliary variables and totals in block i
         n.aux.i <- diff(range(rc.i)) + 1
         auxvar.i <- paste(" - Auxiliary variables: ", n.aux.i, sep = "")
         kntot.i  <- paste(" - Known totals:        ", prod(n.aux.i, nrow(pop.totals)), sep = "")

         # Data.frame with columns order and corresponding auxiliary information in block i
         df.i <- rbind(aux = cols[aux.blocks[[i]]])
         colnames(df.i) <- n.part.vars + aux.blocks[[i]]
        }
     else {
         # Textual description of columns in block i
         c.i <- paste(" - Column:              ", n.part.vars + aux.blocks[[i]], sep = "")
         # Textual description of auxiliary information in block i
         d.i <-       " - Benchmark:           Count of elementary units"
         # How many auxiliary variables and totals in block i
         auxvar.i <- paste(" - Auxiliary variables: ", 1, sep = "")
         kntot.i  <- paste(" - Known totals:        ", nrow(pop.totals), sep = "")

         # Data.frame with columns order and corresponding auxiliary information in block i
         df.i <- rbind(aux = cols[aux.blocks[[i]]])
         colnames(df.i) <- n.part.vars + aux.blocks[[i]]
        }

     # If calibration is partitioned, shift columns number accordingly
     aux.blocks[[i]] <- n.part.vars + aux.blocks[[i]]
     attr(aux.blocks[[i]], "c.i") <- c.i
     attr(aux.blocks[[i]], "d.i") <- d.i
     attr(aux.blocks[[i]], "df.i") <- df.i

#-- Print on screen START --#
     # Vertical separator
     vsep <- paste(rep("-", 0.8 * getOption("width")), collapse = "")
     cat(" - BLOCK", i, vsep, "\n", sep = " ")
     cat(d.i, "\n")
     cat(kntot.i, "\n")
     cat(auxvar.i, "\n")
     cat(c.i, "\n")
     cat("\n")
     print(df.i)
     cat("\n\n")
#-- Print on screen END --#
    }

# Return pop.totals invisibly (as print would do)
invisible(pop.totals)
}
