`trimcal` <-
function (cal.design, w.range = c(-Inf, Inf),
          maxit = 50, epsilon = 1e-07, force = TRUE)
################################################################################
# This function trims calibration weights while preserving all the calibration #
# constraints.                                                                 #
#                                                                              #
# NOTE: The constrained trimming algorithm is not unlike the GEM (Generalized  #
#       Exponential Method) of Folsom & Singh (2000), but adopts the range     #
#       restricted euclidean distance - instead of the logit - for numerical   #
#       stability considerations.                                              #
# NOTE: As for variance estimation, the output object is treated as a *true*   #
#       calibrated object, and the trimmed calibrated weights are regarded as  #
#       if they were obtained in just a single calibration step:               #
#                                                                              #
#                     CAL           TRIM                                       #
#                  w ------> w.cal ------> w.cal.cal <--                       #
#                  |   g1            g2                 |                      #
#                  |   s2            s2                 |                      #
#                   ------------------------------------                       #
#                                CAL(+)TRIM                                    #
#                                   g.tld                                      #
#                                  s2.tld                                      #
#                                                                              #
################################################################################
{
    # First verify if the function has been called inside another function:
    # this is needed to correctly manage metadata when e.g. the caller is a
    # GUI stratum. (affects ONLY ecal.status)
    directly <- !( length(sys.calls()) > 1 )

    if (!inherits(cal.design, "cal.analytic")) 
        stop("Object 'cal.design' must be of class cal.analytic")

    # Trimming an already trimmed object must be avoided, in order to prevent
    # 'trimming cycles'. The user should be advised to use tighter constraints
    # on the calibration weights since the beginning:
    if (is.trimmed(cal.design)) 
        stop("Calibration weights of object 'cal.design' have already been trimmed!\n  If needed, you can trim the *original* object using a tighter 'w.range' from the beginning!")

    # Some checks on w.range
    if ( !is.numeric(w.range) || anyNA(w.range) )
        stop("Argument 'w.range' must be numeric and cannot contain NAs!")
    # If the specified range restrictions do not impose any actual constraints,
    # print a warning and return the input object
    curr.w.range <- range(weights(cal.design))
    if ( (w.range[1] <= curr.w.range[1]) && (curr.w.range[2] <= w.range[2]) ) {
        warning("No point in trimming:\ncalibration weights already fall within the limits specified by 'w.range'!\n", immediate. = TRUE)
        return(cal.design)
    }
    # Empty interval
    if ( w.range[1] >= w.range[2] )
        stop("Argument 'w.range' must satisfy w.range[1] < w.range[2]")
    # Empty overlap with current weights
    if ( (w.range[1] >= curr.w.range[2]) || (w.range[2] <= curr.w.range[1]) )
        stop("No overlap between the current calibration weights and the specified trimming interval 'w.range': cannot succeed!")
    # Trimmed weights must be positive
    if (any(w.range < 0)) {
         stop("If specified, both the trimming limits defined by 'w.range' must be positive!")
        }

    if (!is.logical(force)) 
        stop("Parameter 'force' must be logical")

    # NOTE: ONLY linear trimming metric enabled, TO DATE.
    #       This preference is essentially dictated by numerical stability
    #       considerations.
    # NOTE: Here trimfun is a character: the actual binding with function
    #       `trim.linear` (implementing the Generalized Euclidean Trimming) is
    #       delegated to internal function `ez.trim`
    trimfun <- "linear"

    # Get calibration metadata from cal.design attributes:
    calmodel  <- attr(cal.design, "calmodel")
    partition <- attr(cal.design, "partition")
    sigma2 <- attr(cal.design, "sigma2")
    aggregate.stage <- attr(cal.design, "aggregate.stage")

    # Desume known population totals from available calibration metadata:
    # REMINDER: SHOULD consider swithcing off printing on screen by
    #           aux.estimates, as the message is USELESS here! DONE 22/07/2016
    df.population <- aux.estimates(cal.design, calmodel = calmodel, partition = partition)

    e.df <- cal.design$variables
    calmodel.vars <- all.vars(calmodel)
    na.Fail(e.df, calmodel.vars)
    ids <- attr(cal.design, "ids")
    ids.char <- all.vars(ids)
    stages <- length(ids.char)
    weights <- attr(cal.design, "weights")
    weights.char <- all.vars(weights)
      # Initial weights of previous calibration step 
      weights.back.char <- substr(weights.char, 0, nchar(weights.char) - 4)
      weights.back <- as.formula(paste("~", weights.back.char, sep = ""), env = .GlobalEnv)
    cal.weights.char <- paste(weights.char, ".cal",sep = "")

    if (!is.null(sigma2)){
        sigma2.char <- all.vars(sigma2)
        variance <- e.df[, sigma2.char]
    } else sigma2.char <- NULL

    mk.ecal.status <- function(df.population,partition,partition.names=NULL){
    ############################################################
    #  Diagnostica sul processo di calibrazione.               #
    #  Crea nel .GlobalEnv la lista a 2/3 componenti           #
    #  'ecal.status':                                          #
    #  - 'call':                                               #
    #    identifica la chiamata di e.calibrate che ha          #
    #    generato la lista 'ecal.status';                      #
    #  - 'return.code':                                        #
    #    matrice che contiene i codici di ritorno dei singoli  #
    #    sub-task del processo di calibrazione complessivo:    #
    #     -1 -> task non ancora affrontato;                    #
    #     0  -> convergenza ottenuta;                          #
    #     1  -> convergenza NON ottenuta ma force=TRUE.        #
    # - 'fail.diagnostics':                                    #
    #    presente solo se almeno un return code e' 1, fornisce #
    #    utili informazioni di diagnostica per ogni partizione #
    #    in cui la convergenza sia stata forzata.              #
    ############################################################
        tasks <- 1
        parts <- nrow(df.population)
        rownam <- "code"
        colnam <- if (identical(partition,FALSE)) "global" else partition.names
        ret.cod <- matrix(-1, nrow = tasks, ncol = parts, dimnames = list(rownam, colnam))
        if   (directly) {
               # assign("ecal.status", list(call = sys.call(-1), return.code = ret.cod), envir = .GlobalEnv)
               assign2GE("ecal.status", list(call = sys.call(-1), return.code = ret.cod))
              }
        else  {
               # assign("ecal.status", list(return.code = ret.cod), envir = .GlobalEnv)
               assign2GE("ecal.status", list(return.code = ret.cod))
              }
    }
    upd.ecal.status <- function(n.sub.task,code){
    ###############################################
    #  Aggiorna la lista 'ecal.status' nel        #
    #  .GlobalEnv con il codice 'code' ritornato  #
    #  dal sub-task di ordine 'n.subtask' e, se   #
    #  code=1, con i dati diagnostici.            #
    ###############################################
        last.task <- col(t(ecal.status[["return.code"]]))[n.sub.task]
        last.part <- row(t(ecal.status[["return.code"]]))[n.sub.task]
        ecal.status[["return.code"]][last.task,last.part] <<- code
        # Add fail diagnostics data:
        if (code==1) {
            fd <- list(attr(code, "fail.diagnostics"))
            names(fd) <- colnames(ecal.status[["return.code"]])[last.part]
            ecal.status[["fail.diagnostics"]] <<- c(ecal.status[["fail.diagnostics"]], fd)
        }
    }
    #################################################################
    #  Il parametro need.gc determina se la garbage collection      #
    #  debba, o non debba, essere gestita dal programma.            #
    #  Il valore di soglia per la dimensione della model-matrix     #
    #  completa del modello di calibrazione e' fissato ad 1/10      #
    #  della memoria massima allocabile.                            #
    #  NOTA: Non dividere per il numero di domini di calibrazione   #
    #        nella stima della dimensione massima della             #
    #        model-matrix serve a contrastare potenziali casi       #
    #        "sfortunati" in cui le dimensioni dei domini di        #
    #        calibrazione siano molto squilibrate (cioe' quando     #
    #        esistano domini con un numero di unita' molto maggiore #
    #        che nell'equidistribuzione).                           #
    #################################################################
    need.gc <- FALSE
    if (Sys.info()["sysname"] == "Windows"){
        naux  <- if (identical(partition, FALSE)) {
                     ncol(df.population)
                    }
                 else {
                     (ncol(df.population) - length(all.vars(partition)))
                    }
        nrec <- nrow(e.df)
        # MEM.mega <- memory.limit()
        # See the NOTE on memory.limit() after 4.2.0 in e.calibrate
        MEM.mega <- 4096
        mem.frac <- 10                 #Default value
        need.gc <- ( ((8 * nrec * naux )/(1024^2))  > (MEM.mega / mem.frac) )
    }
    if (need.gc) 
        warning("Complete calibration model-matrix ", if (identical(partition, FALSE)) "takes" else "would take",
                " up more than 0.4 GB of allocable memory", immediate. = TRUE)

    gc.here <- function(doit) {
    #################################################
    #  Se doit=TRUE effettua la garbage collection  #
    #  quando viene invocata.                       #
    #################################################
        if (doit) 
            gc()
    }
    # Require MASS package to get ginv at hand ## No longer needed: ReGenesees now IMPORTS MASS
    # require(MASS)
    cal.sub.task <- 0

    check.totals <- function(tot)
    #######################################
    # "Type" check for population totals. #
    #######################################
    {
     num.col <- function(df) sapply(names(df), function(v) is.numeric(df[, v]))
     if ( any(is.na(as.matrix(tot))) || any(is.infinite(as.matrix(tot))) || !(all(num.col(tot))) ) 
          stop("Population totals must be numeric and finite")
    }

    if (identical(partition,FALSE)) {
    #####################################################
    #  Calibrazione globale (senza domini)              #
    #####################################################
        check.totals(df.population)
        N.cal.constr <- ncol(df.population)
        mk.ecal.status(df.population,partition)
        qr.list <- vector("list", length = 1)
        names(qr.list) <- "population"
        gc.here(need.gc)
        des <- ez.design(data = e.df, weights = weights, ids = ids, 
            variables = c(ids.char, weights.char, weights.back.char, calmodel.vars, sigma2.char))
        gc.here(need.gc)
        w.cal <- ez.trim(design = des, population = as.numeric(df.population), 
            formula = calmodel, trimfun = trimfun, w.range = w.range, 
            aggregate.stage = aggregate.stage, sigma2 = sigma2, weights.back = weights.back,
            maxit = maxit, epsilon = epsilon, force = force)
        gc.here(need.gc)
        cal.sub.task <- (cal.sub.task+1)
        upd.ecal.status(cal.sub.task,attr(w.cal,"ret.code"))
        qr.list[[1]] <- list(group = 1:nrow(e.df), qr = attr(w.cal,"qr"), gwhalf = attr(w.cal,"gwhalf"))
        gc.here(need.gc)
        w.cal <- as.numeric(w.cal)
        gc.here(need.gc)
    }
    else {
    #####################################################
    #  Calibrazione sui singoli domini di calibrazione  #
    #####################################################
        partition.vars <- all.vars(partition)
        check.totals(df.population[, -(1:length(partition.vars)), drop = FALSE])
        na.Fail(e.df, partition.vars)
        N.cal.constr <- prod(nrow(df.population), ncol(df.population)- length(partition.vars))
        partition.names <- apply(df.population[,partition.vars,drop=FALSE],1,paste,collapse=".")
        mk.ecal.status(df.population,partition,partition.names)
        #  'interact': factor i cui livelli identificano le partizioni
        interact <- interaction(e.df[, rev(partition.vars), drop = FALSE], drop=TRUE)
        #  'groups': lista che contiene gli indici di riga delle osservazioni nelle diverse partizioni
        # groups <- .Internal(split(1:nrow(e.df), interact))
        groups <- split(1:nrow(e.df), interact)
        qr.list <- vector("list", length = length(groups))
        names(qr.list) <- partition.names
        #  'w.cal' conterra' i pesi calibrati calcolati nel loop partizione per partizione
        w.cal <- rep(NA,nrow(e.df))
        i.g <- 0
        for (g in groups) {
            i.g <- i.g+1
            des.g <- ez.design(data = e.df[g,], weights = weights, ids = ids, 
              variables = c(ids.char, weights.char, weights.back.char, calmodel.vars, partition.vars, sigma2.char))
            pop.g <- df.population[i.g, which(!(names(df.population) %in% 
              partition.vars))]
            w.cal.g <- ez.trim(design = des.g, population = as.numeric(pop.g), 
              formula = calmodel, trimfun = trimfun, w.range = w.range,
              aggregate.stage = aggregate.stage, sigma2 = sigma2, weights.back = weights.back,
              maxit = maxit, epsilon = epsilon, force = force)
            gc.here(need.gc)
            cal.sub.task <- (cal.sub.task+1)
            upd.ecal.status(cal.sub.task,attr(w.cal.g,"ret.code"))
            qr.list[[i.g]] <- list(group = g, qr = attr(w.cal.g,"qr"), gwhalf = attr(w.cal.g,"gwhalf"))
            w.cal[g] <- as.numeric(w.cal.g)
            }
        gc.here(need.gc)
    }
# Aggiunge al dataframe interno a cal.design i pesi calibrati calcolati w.cal
cal.design$variables[[cal.weights.char]] <- w.cal
cal.design$prob <- 1/w.cal
# NOTE: All ReGenesees function that *change* the weights of design objects (e.g.
#       e.calibrate, trimcal, ...) do *not* update the *design$allprob* slot.
#       Therefore, *design$allprob* acts as a persistent memory of the *initial*
#       weights in arbitrary multi-step weights adjustment procedures!
attr(cal.design, "weights") <- as.formula(paste("~", cal.weights.char, sep = ""), env = .GlobalEnv)
if (!is.null(sigma2)){
    attr(cal.design, "sigma2")  <- sigma2
    }
# Add some component to qr.list (so it becomes conform to caldata in calibrate.survey.design2)
psvar <- list(qr.list = qr.list, stage = 0, index = NULL)
class(psvar) <- "analytic_calibration"
cal.design$postStrata <- list(psvar)
# ADD calibration model metadata (2016/07/21 ADDED TO BENEFIT trimcal)
attr(cal.design, "calmodel") <- calmodel
attr(cal.design, "partition") <- partition
if (!is.null(aggregate.stage)){
    attr(cal.design, "aggregate.stage")  <- aggregate.stage
    }
# Add a token to testify trimming
# NOTE: THIS TOKEN COULD (AND MUST) BE REMOVED BY ANY SUBSEQUENT CALL OF
#       e.calibrate
attr(cal.design, "trimmed") <- TRUE

# Add calibration diagnostics
attr(cal.design, "ecal.status") <- get("ecal.status", envir = .GlobalEnv)
attr(cal.design, "epsilon") <- epsilon
attr(cal.design, "N.cal.constr") <- N.cal.constr
# Define cal.analytic class
if (inherits(cal.design, "cal.analytic")) {
    cal.design$call <- c(sys.call(), cal.design$call)
    gc.here(need.gc)
    }
else {
    cal.design$call <- sys.call()
    class(cal.design) <- c("cal.analytic", class(cal.design))
    gc.here(need.gc)
    }
cal.design
}


is.trimmed <- function(design){
  isTRUE(attr(design, "trimmed"))
}
