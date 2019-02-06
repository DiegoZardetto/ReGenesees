`na.Fail` <-
function (data, varnames)
###########################################################
# Controlla che le variabili i cui nomi sono specificati  #
# da 'varnames' siano 1) effettivamente PRESENTI nel      #
# dataframe 'data' e 2) siano prive di valori mancanti.   #
###########################################################
{
    absent.vars <- varnames[!(varnames %in% names(data))]
    if (length(absent.vars) > 0)
        stop("Variables not found: ", paste(absent.vars, collapse = ", "))
    has.na <- sapply(varnames, function(var) any(is.na(data[, 
        var])))
    if (any(has.na)) 
        stop("Missing values in: ", paste(varnames[has.na], 
              collapse = ", "), "\n")
}
