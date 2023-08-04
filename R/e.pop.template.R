`pop.template` <-
function (data, calmodel, partition = FALSE)
###################################################################
#  Costruisce il template del dataframe dei totali noti per un    #
#  determinato modello di calibrazione.                           #
#  La funzione e.calibrate() esige dataframe dei totali noti      #
#  conformi a tale template.                                      #
#  NOTA: 'calmodel' specifica la formula che definisce il         #
#         modello di calibrazione (in particolare le variabili    #
#         ausiliarie).                                            #
#  NOTA: 'partition' specifica la formula che definisce i factor  #
#        la cui interazione individua le sottopopolazioni per le  #
#        quali sono noti i totali di TUTTE le variabili           #
#        ausiliarie (i "domini di calibrazione").                 #
#        La formula deve essere del tipo:                         #
#          partition=~var1:...:varn                               #
#        (ogni operatore nella formula verra' comunque            #
#        interpretato come ":")                                   #
#        Il valore FALSE per il parametro 'partition' (opzione di #
#        default) indica un modello di calibrazione globale.      #
#  NOTA: Se l'argomento 'partition' viene fornito in modo         #
#        esplicito e' ERRATO riportarne le variabili nella        #
#        formula 'calmodel'. Il parametro 'calmodel' deve         #
#        indicare la formula che definisce il modello di          #
#        calibrazione comune a tutti i domini di calibrazione.    #
#  NOTE: Some instructions were added to enable new function      #
#        'pop.desc'. Such instructions are preceded by string     #
#        # 4 'pop.desc' #.                                        #
###################################################################
{
    # First verify if the function has been called inside another function:
    # this is needed to correctly manage metadata when e.g. the caller is a
    # GUI stratum
    directly <- !( length(sys.calls()) > 1 )

    if (!inherits(data, "data.frame") && !inherits(data, "analytic")) 
        stop("Survey data must be supplied as a data frame or as an object of class 'analytic'")
    data.expr <- if (directly) substitute(data)
    if (inherits(data, "analytic")){
        data <- data$variables
        }
    else {
          # Prevent havoc caused by tibbles:
          if (inherits(data, c("tbl_df", "tbl")))
              data <- as.data.frame(data)
          # Drop empty levels from factor variables (if any)
          # (recall this is unnecessary for analytic objects thanks to e.svydesign)
          data <- emptylev.check(data)
        }

    if (!inherits(calmodel, "formula")) 
        stop("Parameter 'calmodel' must be supplied as a formula")
    calmodel.vars <- all.vars(calmodel)
    na.Fail(data, calmodel.vars)
    cal.mm <- model.matrix(calmodel, model.frame(calmodel, data[1, ])) # 4 'pop.desc' #
    calmodel.names <- colnames(cal.mm)

    partition.expr <- FALSE
    if (!identical(partition, FALSE)) {
        if (!inherits(partition, "formula")) 
            stop("Parameter 'partition' must be supplied as a formula")
        partition.expr <- partition
        partition.vars <- all.vars(partition)
        if (length(partition.vars)<1)
            stop("Parameter 'partition' must reference survey data variables")
        na.Fail(data, partition.vars)
        typetest <- sapply(partition.vars, function(v) is.factor(data[, 
            v]))
        if (!all(typetest)) 
            stop("Partition variables must be factors")
        if (any(partition.vars %in% calmodel.vars)) 
            stop("Calibration model formula cannot reference partition variables")
        u <- unique(data[, partition.vars, drop = FALSE])
        us <- u[do.call(order, u), , drop = FALSE]
        calmodel.tot <- as.data.frame(matrix(data = as.numeric(NA), 
            nrow = nrow(us), ncol = length(calmodel.names), dimnames = list(NULL, 
                calmodel.names)))
        template <- cbind(us, calmodel.tot)
    }
    else {
        template <- as.data.frame(matrix(data = as.numeric(NA), 
            nrow = 1, ncol = length(calmodel.names), dimnames = list(NULL, 
                calmodel.names)))
    }
    row.names(template) <- NULL
    ############################################################
    # Prepare some data structures for new function 'pop.desc' #
    ############################################################
    # 4 'pop.desc' START #
      # Get survey variables appearing in calmodel formula with their "type"
      calmodel.vars.type <- sapply(calmodel.vars, function(var) ifelse(is.numeric(data[[var]]), "numeric", "categorical"))
      # Get the "variable" by "term" matrix
      # NOTE: here "variable" can be a non-survey, computed variable
      var.by.term <- attr(terms(calmodel), "factors")
      # Get computed variables, if any
      comp.vars <- rownames(var.by.term)[!(rownames(var.by.term) %in% calmodel.vars)]
      # Get all variables contributing to template columns related to totals
      # (i.e. exclude partition variables, if any), and determine their "type"
      # NOTE: type could be 'numeric', 'categorical' or 'computed'
      pop.vars <- rownames(var.by.term)
      pop.vars.type <- ifelse(pop.vars %in% comp.vars, "computed", calmodel.vars.type[pop.vars])
      names(pop.vars.type) <- pop.vars
      # For each auxiliaries column of the template, what term did generate it?
      which.term.col <- attr(cal.mm, "assign")
    # 4 'pop.desc' END #

    attr(template, "calmodel")  <- calmodel

    attr(template, "pop.vars.type")  <- pop.vars.type    # 4 'pop.desc' #
    attr(template, "var.by.term")  <- var.by.term        # 4 'pop.desc' #
    attr(template, "calmodel.names")  <- calmodel.names  # 4 'pop.desc' #
    attr(template, "which.term.col")  <- which.term.col  # 4 'pop.desc' #

    attr(template, "partition") <- partition.expr
    attr(template, "data")      <- data.expr
    class(template) <- c("pop.totals", class(template))
    template
}
