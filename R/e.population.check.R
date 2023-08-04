`population.check` <-
function (df.population, data, calmodel, partition = FALSE)
####################################################################
#  Verifica se un dataframe di totali noti di popolazione sia,     #
#  o non sia, conforme allo standard richiesto da e.calibrate      #
#  (la struttura standard e' quella generata dalla funzione        #
#  pop.template).                                                  #
#  Nel primo caso ritorna (senza stamparlo a video) il dataframe   #
#  dei totali noti dopo averlo convertito in un oggetto di classe  #
#  pop.totals, nel secondo caso interrompe l'elaborazione e stampa #
#  un messaggio di errore.                                         #
####################################################################
{
    if (!inherits(df.population, "data.frame")) 
        stop("Known population totals must be supplied as a data frame")
    # Prevent havoc caused by tibbles:
    if (inherits(df.population, c("tbl_df", "tbl")))
        df.population <- as.data.frame(df.population)

    if (!inherits(data, "data.frame") && !inherits(data, "analytic")) 
        stop("Survey data must be supplied as a data frame or as an object of class 'analytic'")
    # Prevent havoc caused by tibbles:
    if (inherits(data, c("tbl_df", "tbl")))
        data <- as.data.frame(data)
    if (inherits(data, "analytic"))
        data <- data$variables

    if (!inherits(calmodel, "formula")) 
        stop("Calibration model must be supplied as a formula")
    if (!identical(partition, FALSE)) {
        if (!inherits(partition, "formula")) 
            stop("Partition variables must be supplied as a formula")
    }
    template <- pop.template(data, calmodel, partition)

    if (!identical(dim(df.population), dim(template))){
        stop.dim <- "Dimension of dataframe 'df.population' does not agree with 'calmodel' and 'partition' formulas\n
                    (to solve the problem use pop.template)"
        stop(stop.dim)
    }
    if (!identical(names(df.population),names(template))){
        stop.names <- "Columns names of data frame 'df.population' do not agree with 'calmodel' and 'partition' formulas\n
                      (to solve the problem use pop.template)"
        stop(stop.names)
    }
    if (!identical(partition, FALSE)) {
        test.var.class <- function(df, class) sapply(names(df), 
            function(v) inherits(df[, v], class))
        template.factor <- data.frame(template[, test.var.class(template, 
            "factor"), drop = FALSE])
        df.population.factor <- data.frame(df.population[, test.var.class(df.population, 
            "factor"), drop = FALSE])
        df.fm <- as.matrix(df.population.factor)
        # Strip rownames from df.population (those of template are NULL, and
        # a rownames mismatch - if any - has to be tolerated)
        row.names(df.fm) <- NULL
        if (!identical(df.fm, as.matrix(template.factor))){
            stop.fact <- "Columns of data frame 'df.population' defining calibration domains\n
                         do not agree with 'calmodel' and 'partition' formulas\n
                         (to solve the problem use pop.template)"
            stop(stop.fact)
        }
    }
    template[, ] <- df.population
    cat("\n# Checking Known Totals dataframe: OK\n\n")
    invisible(template)
}
