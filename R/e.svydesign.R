`e.svydesign` <-
function(data, ids, strata = NULL, weights,
         fpc = NULL, self.rep.str = NULL, check.data = TRUE){
########################################################################
# Defines an object of class 'analytic' which is an extension of the   #
# 'survey.design2' class from package survey.                          #
# NOTE: If passed, all arguments referencing variables must be         #
#       formulas (unlike in the original svydesign function)           #
# NOTE: Some original arguments (i.e. nest, variables, check.strata,   #
#       probs) have been intentionally dropped for the purpose of      #
#       developing a GUI application. Conversely, new arguments        #
#       self.rep.str and check.data have been introduced.              #
########################################################################

# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )

# Just some checks on supplied arguments
    if (!inherits(data, "data.frame")) 
        stop("Survey data must be supplied as a data frame")

    data.expr <- if (directly) substitute(data)

    # Prevent havoc caused by tibbles:
    if (inherits(data, c("tbl_df", "tbl")))
        data <- as.data.frame(data)

    # Drop empty levels from factor variables (if any)
    data <- emptylev.check(data)

    if (!inherits(ids,"formula"))
        stop("Cluster identifiers must be supplied as a formula")
    ids.charvect <- all.vars(ids)
    if (length(ids.charvect)<1)
        stop("Cluster identifiers must reference survey data variables")
    na.Fail(data, ids.charvect)
    
    if (!inherits(weights,"formula"))
        stop("Weights must be passed as a formula")
    weights.char <- all.vars(weights)
    if (length(weights.char)<1) 
        stop("Weights formula must reference a survey data variable")
    if (length(weights.char)>1) 
        stop("Weights formula must reference only one variable")
    na.Fail(data, weights.char)
    if (!is.numeric(data[, weights.char])) 
        stop("Weights variable ", weights.char, " is not numeric")

    ### BLOCK BELOW IS COMMENTED BECAUSE IT DOES NOT WORK YET ###
    # Prevent negative direct weights, with the possible EXCEPTION of those
    # coming from an invokation of ext.calibrated(), as they are actually
    # CALIBRATED external weights (see ext.calibrated)
    # if (any(data[, weights.char] <= 0) && !isTRUE(attr(data, "negw.pass"))) 
    ### BLOCK ABOVE IS COMMENTED BECAUSE IT DOES NOT WORK YET ###

    if (any(data[, weights.char] <= 0)) 
        stop("Direct weights must be positive!")


    strata.expr <- FALSE
    if (!is.null(strata)){
        if (!inherits(strata,"formula"))
            stop("If supplied, strata must be passed as a formula")
        strata.expr <- strata
        }
    
    fpc.expr <- FALSE
    fpc.charvect <- NULL
    few.fpc <- FALSE
    if (!is.null(fpc)){
        if (!inherits(fpc,"formula"))
            stop("If supplied, fpc must be passed as a formula")
        fpc.charvect <- all.vars(fpc)
        if (length(fpc.charvect)<1) 
            stop("fpc formula must reference survey data variables")
        na.Fail(data, fpc.charvect)
        typetest <- sapply(fpc.charvect, function(y) is.numeric(data[, 
            y]))
        if (!all(typetest)) 
            stop("'fpc' variables must be numeric")
        # Check for fpc stages
        if (length(fpc.charvect) > length(ids.charvect)) stop("Too many stages in fpc formula")
        if (length(fpc.charvect) < length(ids.charvect)) few.fpc <- TRUE
        fpc.expr <- fpc
    }

    self.rep.str.expr <- FALSE

    check.nest <- function(data, ids, strata) {
    #################################################
    #  Controlla che le unita' di campionamento di  #
    #  stadio k+1 siano correttamente innestate     #
    #  all'interno di quelle di stadio k.           #
    #  Nota: Il controllo si estende agli strati,   #
    #        che sono trattati come unita' di       #
    #        di stadio k=0.                         #
    #################################################
        nested.unit <- ids
        unit <- c(strata, nested.unit)
        for (i in 1:length(nested.unit)) {
            # Following conversion for speeding up tapply
            ext.unit <- if ( is.factor(data[, unit[i]]) ) as.numeric(data[, unit[i]]) else data[, unit[i]]
            if (any(tapply(ext.unit, data[, 
                nested.unit[i]], function(x) length(unique(x)) > 
                1))) {
                if (i == 1) {
                  stop("There is at least one PSU belonging to different strata")
                }
                else {
                  stop("There is at least one stage ", 
                    i, " cluster belonging to different stage ", 
                    i - 1," clusters")
                }
            }
        }
    }

    selfrep <- function(data, ids.charvect, fpc.charvect = NULL, self.rep.str, strata.char) {
    ####################################################
    # Crea la variabile var.PSU (identificativo delle  #
    # unita' che forniscono il contributo leading alla #
    # varianza) per disegni di campionamento che       #
    # includono strati auto-rappresentativi (uno per   #
    # ognuna PSU campionata con certezza).             #
    # Se specificate, aggiorna le fpc: negli strati    #
    # autorappresentativi il contributo di primo       #
    # stadio e' dato dalle fpc del secondo stadio.     #
    ####################################################
    # Check that self.rep.str has only values 0 and 1 when it is numeric or factor
    srs <- data[, self.rep.str]
    if (is.numeric(srs)) {
         if (!all(srs %in% 0:1))
             stop("self.rep.str values can be only 0 and 1 or TRUE and FALSE")
        }
    if (is.factor(srs)) {
         if (!all(as.numeric(levels(srs)) %in% 0:1))
             stop("self.rep.str values can be only 0 and 1 or TRUE and FALSE")
         srs <- as.numeric(levels(srs))[srs]
        }
    # Check that each stratum is either SR or not-SR
    if (any(err <- tapply(srs, data[, strata.char],
                          function(x) length(unique(x)) > 1))) {
            first.err <- names(err[err])[1]
            stop("There is at least one stratum with varying self.rep.str! (e.g. stratum ", first.err,")")            
            }
        var.PSU <- as.character(data[, ids.charvect[1]])
        self <- as.logical(srs)
        self.ri <- which(self)
        var.PSU[self.ri] <- paste(data[self.ri, ids.charvect[1]],
                                  data[self.ri, ids.charvect[2]], sep=".")
        var.FPC <- NULL
        if (!is.null(fpc.charvect)) {
            var.FPC <- data[, fpc.charvect[1]]
            var.FPC[self.ri] <- data[self.ri, fpc.charvect[2]]
        }
        list(var.PSU = var.PSU, var.FPC = var.FPC)
    }
    
# Check for correct nesting of PSUs in strata (if any) and of clusters at different stages
    if (!is.null(strata)){
        strata.char  <- all.vars(strata)
        if (length(strata.char)<1) 
            stop("Strata formula must reference a survey data variable")
        if (length(strata.char)>1) 
            stop("Strata formula must reference only one variable")
        na.Fail(data, strata.char)
        if (!is.factor(data[, strata.char])) 
            stop("Strata variable ", strata.char, " is not a factor")
        if (check.data) 
            check.nest(data, ids.charvect, strata.char)
        # If there are self representing strata and the ultimate cluster approximation applies,
        # build variance PSUs and modify first stage fpc data (if any)
        # NOTE: Even though (1) next to first stage fpc data (if any) are kept
        #       as they were in the user's call, AND (2) even if the ultimate
        #       cluster approximation is *switched off*, design with variance PSUs
        #       will in any case undergo a *single stage variance estimation* (see
        #       svyrecvar).
        if (!is.null(self.rep.str)) {
            if (!inherits(self.rep.str, "formula")) 
                stop("If supplied, self.rep.str must be passed as a formula")
            srs.char <- all.vars(self.rep.str)
            if (length(srs.char) < 1) 
                stop("self.rep.str formula must reference a survey data variable")
            if (length(srs.char) > 1) 
                stop("self.rep.str formula must reference only one variable")
            na.Fail(data, srs.char)
            if (!is.logical(data[, srs.char]) && !is.numeric(data[, srs.char]) && !is.factor(data[, srs.char]))
                stop("Variable self.rep.str must be logical or numeric or factor")
            if (length(ids.charvect)==1) {
                self.rep.str <- NULL
            }
            else {
                if (few.fpc) 
                    stop("If both fpc and self.rep.str are passed, fpcs must be specified for each sampling stage!")
                this.selfrep <- selfrep(data, ids.charvect, fpc.charvect, srs.char, strata.char)
                data[["var.PSU"]] <- this.selfrep$var.PSU
                if (!is.null(fpc)) {
                    data[[fpc.charvect[1]]] <- this.selfrep$var.FPC
                    }
            }
        }
    }
    else {
        data[["strata.default"]] <- as.factor(1)
        if (check.data && (length(ids.charvect) > 1)) 
            check.nest(data, ids.charvect, "strata.default")
        self.rep.str <- NULL
    }

# If variance PSUs have been built, warn user and update ids formula
    if (!is.null(self.rep.str)) {       
        id.PSU <- "var.PSU"
        ids <- as.formula(paste("~", paste(c(id.PSU, ids.charvect[-1]), collapse = "+")), env = .GlobalEnv)
        self.rep.str.expr <- self.rep.str
    }

    edes <- svydesign(ids = ids, probs = NULL, strata = strata, variables = NULL,
                      fpc = fpc, data = data, nest = FALSE, check.strata = FALSE,
                      weights = weights)

# Just some little modification of the survey.design2 object

    # Catch the actual call
    edes$call<-sys.call()

    # Attach to the would-be survey.design2 object some metadata
    # that migth be used by e.calibrate (if needed)
    attr(edes,"data")         <- data.expr
    attr(edes,"ids")          <- ids 
    attr(edes,"strata")       <- strata.expr
    attr(edes,"weights")      <- weights
    attr(edes,"fpc")          <- fpc.expr
    attr(edes,"self.rep.str") <- self.rep.str.expr

    # Eventually extend survey.design2 class
    class(edes) <- c("analytic", class(edes))
    
    # If SR strata specified (and no errors have been generated till now)
    # put a reminder for the user:
    if (!is.null(self.rep.str)){
       warning("Sampling variance estimation for this design will take into account only leading contributions, i.e. PSUs in not-SR strata and SSUs in SR strata (see ?e.svydesign and ?ReGenesees.options for details)")
       }
    edes
}
