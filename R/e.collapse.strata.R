collapse.strata <- function(design, block.vars = NULL, sim.score = NULL){
##############################################################################
# This function modifies a design object (hineriting from class analytic) in #
# such a way that lonely PSU strata (lonely strata, for short) are collapsed #
# by pairing them with the most similar lonely stratum in each block.        #
# NOTE: Blocks are defined by blocking variables trough factor crossing:     #
#       obviously they cannot cut across strata.                             #
# NOTE: Not aggregable strata, i.e. lonely strata which are a singleton      #
#       inside a block, generate an error.                                   #
# NOTE: Strata similarity score should be positive (othervise its abs() will #
#       be used).                                                            #
##############################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum. (affects ONLY clps.strata.status)
directly <- !( length(sys.calls()) > 1 )
this.call <- sys.call()

if (!inherits(design, "analytic")) 
    stop("Design object must be of class analytic")

if (!is.null(block.vars)) {
     if (!inherits(block.vars, "formula")) 
         stop("Blocking variables must be supplied as a formula")
     block.vars <- all.vars(block.vars)
     na.Fail(design$variables, block.vars)
     typetest <- sapply(block.vars, function(v) is.factor(design$variables[, v]))
     if (!all(typetest)) 
         stop("Blocking variables must be factors")
    }

if (!is.null(sim.score)) {
     if (!inherits(sim.score,"formula"))
         stop("Strata similarity score must be passed as a formula")
     sim.score <- all.vars(sim.score)
     if (length(sim.score)>1) 
         stop("Similarity score formula must reference only one variable")
     na.Fail(design$variables, sim.score)
     if (!is.numeric(design$variables[, sim.score])) 
         stop("Similarity score variable must be numeric")
    }

# Check if design is actually stratified
if (!isTRUE(design$has.strata))
   stop("Design object is not stratified!")

old.strata <- design$strata[, 1]
strata.char <- all.vars(attr(design, "strata"))

## 1) Find names of strata with lonely PSUs (if any)

analyze.strata <- function(design){
########################################################
# Processes the fpc slot of a design (i.e. design$fpc) #
# in such a way that a two components list storing PSU #
# counts inside strata for both population and sample  #
# (popsize and sampsize slots) is returned.            #
########################################################
des.fpc <- design$fpc
des.strat <- design$strata[, 1]
strata.PSUcounts <- function(PSU, strata){
tapply(PSU, strata, unique)
}
lapply(des.fpc, function(el) { if (is.null(el))
                                  # i.e. for popsize component,
                                  # when no fpc have been specified
                                  return(el)
                               else strata.PSUcounts(el[, 1], des.strat) })
}

find.lPSU <- function(analysis){
####################################################
# Given the output of analyze.strata, finds strata #
# with lonely PSUs.                                #
# If no lonely PSUs are found gives an error.      #
####################################################
  if (is.null(pop <- analysis$popsize)) {
     # i.e. when no fpc have been specified
     samp <- analysis$sampsize
     if  (any(lonely <- (samp==1))) {
         str.lev <- names(samp)
         return(str.lev[lonely])
        }
     else {
         stop("No point in strata collapsing: no lonely PSUs found!")
        }
    }
  else {
         samp <- analysis$sampsize
         str.lev <- names(samp)
         # exlude certainty PSUs
         lonely <- ((samp==1) & (pop>1))
         if  (any(lonely)) {
             return(str.lev[lonely])
            }
         else{
              stop("No point in strata collapsing: no lonely PSUs found!")
            }
    }
stop("This should not happen! Please report this message to the author (zardetto@istat.it)")
}

des.lonelies <- find.lPSU(analyze.strata(design))

## 1) ENDS

## 2) Find and pre-process blocks containing strata with lonely PSUs

split.blocks <- function(design, block.vars = NULL, score = NULL, lonely){
##########################################################
# Given a stratified design, a set of blocking variables #
# names, a set of strata names with lonely PSUs, and a   #
# similarity score for strata, returns a list containing #
# strata names, lonely strata names and scores belonging #
# to the blocks.                                         #
# NOTE: If score varies inside strata gives an error.    #
# NOTE: If blocks cut across strata gives an error.      #
# NOTE: If any block exists with a single stratum (i.e.  #
#       a lonely stratum needing to be collapsed) gives  #
#       an error.                                        #
##########################################################
data <- design$variables
str  <- design$strata[, 1]
if (is.null(score)) {
     scr <- rep(1, nrow(data))
    }
else {
     # ensure that score is positive
     scr <- abs(data[, score])
     # Check that score does not vary inside strata
     tab <- tapply(scr, str, function(x) length(unique(x)))
     if (any(tab > 1)) {
         err.msg <- "Similarity score varies inside strata!"
         err.strata <- names(tab[tab > 1])
         if (directly) {
             # assign("clps.strata.status",
             #        list(error = err.msg, ko.strata = err.strata, call = this.call),
             #        envir = .GlobalEnv)
             assign2GE("clps.strata.status",
                       list(error = err.msg, ko.strata = err.strata, call = this.call))
            }
         else {
             # assign("clps.strata.status",
             #        list(error = err.msg, ko.strata = err.strata),
             #        envir = .GlobalEnv)
             assign2GE("clps.strata.status",
                       list(error = err.msg, ko.strata = err.strata))
            }
         stop(err.msg)
        }
    }

if (is.null(block.vars)) {
     interact <- factor(rep("G", nrow(data)))
    }
else {
     interact <- interaction(data[, block.vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
     # Check that strata are actually nested inside blocks
     # Following conversions for speeding up tapply
     blocker <- as.numeric(interact)
     tab <- tapply(blocker, str, function(x) length(unique(x)) )
     if (any(tab > 1)) {
         err.msg <- "Blocks cut across strata!"
         err.strata <- names(tab[tab > 1])
         if (directly) {
             # assign("clps.strata.status",
             #        list(error = err.msg, ko.strata = err.strata, call = this.call),
             #             envir = .GlobalEnv)
             assign2GE("clps.strata.status",
                       list(error = err.msg, ko.strata = err.strata, call = this.call))
            }
         else {
             # assign("clps.strata.status",
             #        list(error = err.msg, ko.strata = err.strata),
             #        envir = .GlobalEnv)
             assign2GE("clps.strata.status",
                       list(error = err.msg, ko.strata = err.strata))
            }
         stop(err.msg)
        }
    }

# blocks <- .Internal(split(1:nrow(data), interact))
blocks <- split(1:nrow(data), interact)
blocks.with.lonely <- blocks[sapply(blocks, function(bl) any(str[bl] %in% lonely))]
strata.blocks <- lapply(blocks.with.lonely, function(bl) as.character(unique(str[bl])))

# Check if some block has only one lonely PSU stratum (i.e. not aggregable)
check <- sapply(strata.blocks, function(bl) length(bl)==1)
if (any(check)) {
     err.msg <- "Some blocks contain just a single stratum: cannot collapse it!"
     err.blocks <- names(check[check])
     if (directly) {
         # assign("clps.strata.status",
         #        list(error = err.msg, ko.blocks = err.blocks, call = this.call),
         #        envir = .GlobalEnv)
         assign2GE("clps.strata.status",
                   list(error = err.msg, ko.blocks = err.blocks, call = this.call))
        }
     else {
         # assign("clps.strata.status",
         #        list(error = err.msg, ko.blocks = err.blocks),
         #        envir = .GlobalEnv)
         assign2GE("clps.strata.status",
                   list(error = err.msg, ko.blocks = err.blocks))
        }
     stop(err.msg)
    }

lonely.blocks <- lapply(strata.blocks, function(bl) lonely[lonely %in% bl])
score.blocks <- lapply(strata.blocks, function(bl) sapply(bl, function(stratum) unique(scr[str==stratum])))
list(strata=strata.blocks, lonely=lonely.blocks, score=score.blocks)
}

lon.blocks  <- split.blocks(design, block.vars, sim.score, des.lonelies)

## 2) ENDS

## 3) Collapse strata with lonely PSUs by pairing them with the most similar stratum in each block

collapse.block <- function(strata, lonely, score, block.flag){
###############################################################
# Given a vector of strata names, a vector of names of strata #
# containing a lonely PSU, and a vector of similarity scores  #
# for strata, returns a matrix mapping strata to SUPERSTRATA. #
# Strata are paired according to their similarity.            #
# NOTE: strata and score must be the same length.             #
# NOTE: lonely must have length >= 1 and <= of the length of  #
#       strata.                                               #
# NOTE: if length(strata)==1 (i.e. only a lonely PSU stratum  #
#       inside the block) gives an error.                     #
###############################################################
n <- length(strata)
if (n==1) stop("Just a single stratum in block: cannot collapse it!")
strata.ind <- 1:n
lonely.ind <- which(strata %in% lonely)
nlonely <- length(lonely.ind)
superstrata <- strata
busy <- NULL
assigned <- rep(NA, nlonely)
i.lev <- 0
j.tag <- 0
for (lon.i in lonely.ind) {
     if (lon.i %in% busy) next
     i.lev <- i.lev + 1 
     candidates <- setdiff(strata.ind, union(lon.i, busy))
     if (length(candidates) < 1) {
         all.others <- strata.ind[-lon.i]
         assigned[i.lev] <- cur <- all.others[which.min(abs(score[all.others] - score[lon.i]))]
         superstrata[lon.i] <- superstrata[cur]
        }
     else {
         j.tag <- j.tag + 1
         assigned[i.lev] <- cur <- candidates[which.min(abs(score[candidates] - score[lon.i]))]
         superstrata[lon.i] <- superstrata[cur] <- paste(block.flag, j.tag, sep=".clps.")
         busy <- union(busy, union(lon.i, cur))
        }
    }

# Check that each lonely stratum has been actually aggregated 
if (any(sapply(lonely.ind, function(i) !(superstrata[i] %in% superstrata[-i]))))
    stop("Strata aggregation failed!")

out <- matrix(c(strata, superstrata), nrow=n)
colnames(out) <- c("strata", "superstrata")
out
}

nblocks <- length(lon.blocks[[1]])
str.map <- lapply(1:nblocks, function(bl.i){
                                            collapse.block(lon.blocks$strata[[bl.i]],
                                                           lon.blocks$lonely[[bl.i]],
                                                           lon.blocks$score[[bl.i]],
                                                           names(lon.blocks$strata[bl.i])
                                                           )
                                           }
                 )
names(str.map) <- names(lon.blocks[[1]])

## 3) ENDS

## 4) Use the mapping above (str.map) to build superstrata

superstrata <- as.character(old.strata)
for (block in str.map) {
     mapply(block[, "strata"], block[, "superstrata"],
             FUN = function(stratum, superstratum) superstrata[superstrata == stratum] <<- superstratum
            )
    }
superstrata <- data.frame(factor(superstrata))
superstrata.name <- paste(strata.char,".collapsed",sep="")
names(superstrata) <- superstrata.name

## 4) ENDS

## 5) Build new strata and fpc slot for design by using superstrata,
##    add the superstrata column to internal data

    # interaction function
    interaction<-function (..., drop = TRUE) {
        args <- list(...)
        narg <- length(args)
        if (narg == 1 && is.list(args[[1]])) {
            args <- args[[1]]
            narg <- length(args)
        }
        
        ls<-sapply(args,function(a) length(levels(a)))
        ans<-do.call("paste",c(lapply(args,as.character),sep="."))
        ans<-factor(ans)
        return(ans)
        
    }

# Process and update strata slot
ids <- design$cluster
N <- ncol(ids)
new.strata <- superstrata
NS <- ncol(new.strata)
  if (N > 1) {
     for (i in 2:N) new.strata[, i] <- interaction(new.strata[, min(i, NS)], ids[, i - 1])
    }

design$strata <- new.strata

# Process and update fpc slot (only first stage is involved, i.e. PSUs)
new.fpc <- old.fpc <- design$fpc
for (i in seq_along(new.fpc)){
     if (!is.null(new.fpc[[i]])) {
         # i.e. don't process popsize component,
         # when no fpc have been specified
         new.fpc[[i]][, 1] <- ave( ave(new.fpc[[i]][, 1], old.strata, FUN = function(x) x/length(x)),
                                   superstrata, FUN = sum)

         # This way the updated column of the i-th slot of fpc would numeric,
         # whereas ordinary designs MAY have it integer (e.g. for popsize it
         # depends on the way fpc was specified also for ordinary designs),
         # hence make new.fpc type identical to the old:
         if (is.integer(old.fpc[[i]][, 1])) {
             new.fpc[[i]][, 1] <- round(new.fpc[[i]][, 1])
            }
        }
    }

design$fpc <- new.fpc

# Add superstrata column to inner data
design$variables[[superstrata.name]] <- design$strata[, 1]

## 5) ENDS

## 6) Build a table mapping collapsed strata to superstrata

################################################################################
#clps.table2 <- unique(design$variables[, c(strata.char, superstrata.name)])
#changed <- as.character(clps.table2[, 1]) != as.character(clps.table2[, 2])
#clps.table2 <- clps.table2[changed, ]
#lon <- as.integer(clps.table2[, 1] %in% des.lonelies)
#clps.table2 <- data.frame(str = clps.table2[[1]], lon = lon, sstr = clps.table2[[2]])
#names(clps.table2) <- c(strata.char, "lonely", superstrata.name)
#clps.table2 <- clps.table2[order(clps.table2[, superstrata.name]), ]
#rownames(clps.table2) <- NULL
#assign("clps.table2", clps.table2, .GlobalEnv)
################################################################################
blocks.v <- rep(names(str.map), sapply(str.map, nrow))
strata.v <- as.character(unlist(lapply(str.map, function(bl) bl[, "strata"])))
superstrata.v <- as.character(unlist(lapply(str.map, function(bl) bl[, "superstrata"])))
lonely.v <- as.integer(strata.v %in% des.lonelies)
clps.table <- data.frame(b = blocks.v,
                         s = strata.v,
                         l = lonely.v,
                         ss = superstrata.v,
                         stringsAsFactors = FALSE
                        )
# exclude unchanged strata (if any)
changed <- ( strata.v != superstrata.v )
clps.table <- clps.table[changed, ]
names(clps.table) <- c("block", strata.char, "lonely", superstrata.name)
# order by superstratum
clps.table <- clps.table[order(clps.table[, superstrata.name]), ]
rownames(clps.table) <- NULL

msg <- paste("All lonely strata (", length(des.lonelies), ") successfully collapsed!", sep="")
cat(paste("\n# ", msg, "\n\n", sep=""))
if (directly) {
     # assign("clps.strata.status",
     #        list(message = msg, clps.table = clps.table, call = this.call),
     #        envir = .GlobalEnv)
     assign2GE("clps.strata.status",
               list(message = msg, clps.table = clps.table, call = this.call))
    }
else  {
     # assign("clps.strata.status",
     #        list(message = msg, clps.table = clps.table),
     #        envir = .GlobalEnv)
     assign2GE("clps.strata.status",
               list(message = msg, clps.table = clps.table))
    }

## 6) ENDS

# Warn if sim.score has not been passed (and no other exceptions have been raised)
if (is.null(sim.score))
    warning("No similarity score specified: achieved strata aggregation depends on the ordering of sample data")

attr(design,"collapse.strata") <- this.call
design
}
