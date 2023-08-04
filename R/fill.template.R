`fill.template` <-
function(universe, template, mem.frac = 10)
################################################################################
# Given a complete list of units from a population ('universe') and a template #
# for the totals of the auxiliary variables in a calibration task ('template') #
# FILLS the template by ACTUAL population totals computed from 'universe'.     #
# NOTE: If the data to be processed (whose size depend on both universe and    #
#       template) become too big as compared to the maximum allocable memory,  #
#       universe gets split into smaller chunks (see below).                   #
#       Parameter 'mem.frac' fixes the memory threshold to 1/mem.frac of the   #
#       maximum allocable memory. Actual choice mem.frac=10 seems very good.   #
################################################################################
{
    if (!inherits(universe, "data.frame"))
        stop("Universe must be a data frame")
    # Prevent havoc caused by tibbles:
    if (inherits(universe, c("tbl_df", "tbl")))
        universe <- as.data.frame(universe)
    # Drop empty levels from factor variables (if any)
    universe <- emptylev.check(universe)

    if (!inherits(template, "pop.totals"))
        stop("Template must be of class pop.totals")

    if (!is.numeric(mem.frac) || (mem.frac < 0))
        stop("Parameter mem.frac must be numeric and non-negative")

    # Coherence check between template and universe
    universe.check(universe, template)

    # Missing data check on universe
    calmodel <- attr(template, "calmodel")
    calmodel.vars <- all.vars(calmodel)
    na.Fail(universe, calmodel.vars)

    partition <- attr(template, "partition")
    if (!identical(partition, FALSE)) {
        partition.vars <- all.vars(partition)
        na.Fail(universe, partition.vars)
        # Prepare partition names
        partition.names <- apply(template[,rev(partition.vars),drop=FALSE],1,paste,collapse=".")
    }
    
    #################################################################
    #  Il parametro need.fetch determina se i record di 'universe'  #
    #  debbano, o non debbano, essere processati in blocchetti.     #
    #  Il valore di soglia e' determinato stimando la dimensione    #
    #  massima delle model-matrix da calcolare: se la memoria       #
    #  necessaria eccede il 1/mem.frac della memoria massima        #
    #  allocabile il sistema splitta. Il numero di chunks generati  #
    #  e' tale che la matrice di cui sopra generata sui singoli     #
    #  chunks non ecceda 1/(2*mem.frac) della memoria massima.      #
    #  NOTA: Non dividere per il numero di domini di calibrazione   #
    #        nella stima della dimensione massima della             #
    #        model-matrix serve a contrastare potenziali casi       #
    #        "sfortunati" in cui le dimensioni dei domini di        #
    #        calibrazione siano molto squilibrate (cioe' quando     #
    #        esistano domini con un numero di unita' molto maggiore #
    #        che nell'equidistribuzione).                           #
    #  NOTA: Il default di mem.frac e' 10.                          #
    #################################################################
    need.fetch <- FALSE
    if (Sys.info()["sysname"] == "Windows"){
        naux  <- if (identical(partition, FALSE)) ncol(template) else (ncol(template) - length(partition.vars))
        nrec <- nrow(universe)
        # MEM.mega <- memory.limit()
        # See the NOTE on memory.limit() after 4.2.0 in e.calibrate
        MEM.mega <- 4096
    #   mem.frac <- 10 !Taken as default value: user must be able to change it if necessary!
        need.fetch <- ( ((8 * nrec * naux )/(1024^2))  > (MEM.mega / mem.frac) )
    }

    if (need.fetch) {
        nchunks <- ceiling( 2 * ((8 * nrec * naux )/(1024^2)) * (mem.frac / MEM.mega) )
        warning("Processing 'universe' as a whole would require more than 1/",
                mem.frac," of 0.4 GB:\n it will be split into ",
                nchunks, " chunks instead.", immediate. = TRUE)
        if (!is.numeric(nchunks) || (nchunks < 1) || (nchunks > nrec)){
            stop("Chunks number doesn't satisfy 1 <= nchunks <= ", nrec,"!")
        }

        # Prepare tools needed to split universe into chunks
        # 'splitter': factor i cui livelli identificano i chunks con un sequenziale
        splitter <- factor(rep(1:nchunks, each=ceiling(nrec/nchunks), length.out=nrec))
        # 'chunks': lista che contiene gli indici di riga delle osservazioni nei diversi chunks
        # chunks <- .Internal(split(1:nrec, splitter))
        chunks <- split(1:nrec, splitter)
    }


    if (identical(partition,FALSE)) {
    ########################
    #  Global calibration  #
    ########################
        if (!need.fetch){
            mm <- model.matrix(calmodel, model.frame(calmodel, data = universe))
            template[,] <- colSums(mm)
        }
        else {
            ##################################
            # Process chunk-by-chunk and add #
            # partial sums to template.      #
            ##################################
            # Initialize template total to zero
            template[,] <- 0
            for (ch in chunks) {
                mm.ch <- model.matrix(calmodel, model.frame(calmodel, data = universe[ch, , drop = FALSE]))
                template[,] <- template[,] + colSums(mm.ch)
            }
        }
    }
    else {
    #############################
    #  Partitioned Calibration  #
    #############################
        if (!need.fetch){
            # 'interact': factor i cui livelli identificano le partizioni
            interact <- interaction(universe[, rev(partition.vars), drop = FALSE], drop=TRUE)
            # 'groups': lista che contiene gli indici di riga delle osservazioni nelle diverse partizioni
            # groups <- .Internal(split(1:nrow(universe), interact))
            groups <- split(1:nrow(universe), interact)
            group.names <- names(groups)
            ######################################
            # Calcola i totali per le partizioni #
            ######################################
            i.g <- 0
            for (g in groups) {
                i.g <- i.g+1
                mm.g <- model.matrix(calmodel, model.frame(calmodel, data = universe[g, , drop = FALSE]))
                template[partition.names == group.names[i.g],
                         which(!(names(template) %in% partition.vars))] <- colSums(mm.g)
            }
        }
        else {
            ##################################
            # Process chunk-by-chunk and add #
            # partial sums to template.      #
            ##################################
            # Initialize template total to zero for auxiliary variables columns
            aux.ind <-which(!(names(template) %in% partition.vars))
            template[, aux.ind] <- 0
            for (ch in chunks) {
                universe.ch <- universe[ch, , drop = FALSE]
                # 'interact': factor i cui livelli identificano le partizioni
                interact.ch <- interaction(universe.ch[, rev(partition.vars), drop = FALSE], drop=TRUE)
                # 'groups': lista che contiene gli indici di riga delle osservazioni nelle diverse partizioni
                # groups.ch <- .Internal(split(1:nrow(universe.ch), interact.ch))
                groups.ch <- split(1:nrow(universe.ch), interact.ch)
                group.names.ch <- names(groups.ch)
                # Process group-by-group inside current chunk
                i.g <- 0
                for (g in groups.ch) {
                    i.g <- i.g+1
                    mm.g <- model.matrix(calmodel, model.frame(calmodel, data = universe.ch[g, , drop = FALSE]))
                    template[partition.names == group.names.ch[i.g], aux.ind] <-
                    template[partition.names == group.names.ch[i.g], aux.ind] + colSums(mm.g)
                }
            }
        }
    }
 class(template) <- unique(c("pop.totals", class(template)))
 return(template)
}


`universe.check` <-
function (universe, pop.totals)
{
    calmodel <- attr(pop.totals, "calmodel")
    partition <- attr(pop.totals, "partition")

    template <- pop.template(universe, calmodel, partition)

    if (!identical(dim(pop.totals), dim(template))){
        stop.dim <- paste("Dimension mismatch: auxiliary variables in universe do not agree with template\n(to solve the problem use pop.template)")
        stop(stop.dim)
    }

    if (!identical(names(pop.totals), names(template))){
        stop.names <- paste("Names of auxiliary variables in universe do not agree with template\n(to solve the problem use pop.template)")
        stop(stop.names)
    }

    if (!identical(partition, FALSE)) {
        test.var.class <- function(df, class) sapply(names(df), 
            function(v) inherits(df[, v], class))
        template.factor <- data.frame(template[, test.var.class(template, 
            "factor"), drop = FALSE])
        pop.totals.factor <- data.frame(pop.totals[, test.var.class(pop.totals, 
            "factor"), drop = FALSE])
        if (!identical(as.matrix(pop.totals.factor), as.matrix(template.factor))){
            stop.fact <- paste("Template columns defining calibration domains do not agree with universe\n(to solve the problem use pop.template)")
            stop(stop.fact)
        }
    }

    cat("\n# Coherence check between 'universe' and 'template': OK\n\n")
}
