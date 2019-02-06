`des.addvars` <-
function(design, ...)
##########################################################
#  Aggiorna i dati contenuti in un oggetto analytic.     #
#  Gli argomenti ... devono essere del tipo tag = expr   #
#  (un tag puo' essere o un identificatore o una stringa #
#  di caratteri).                                        #
#  La nuova variabile creata avra' nome 'tag' e valori   #
#  ottenuti valutando 'expr' su design.                  #
#  NOTA: Eventuali espressioni sprovviste di tag in ...  #
#        vengono ignorate e non hanno, dunque, nessun    #
#        effetto sul valore di ritorno di des.addvars.   #
##########################################################
{
    if (!inherits(design, "analytic")) 
        stop("Object 'design' must inherit from class analytic")
    dots <- substitute(list(...))[-1]
    newnames <- names(dots)
    if (length(dots) < 1) 
        return(design)
    if (is.null(newnames)) {
        # Solo espressioni prive di tag in (...): esco.
        warning("Untagged input expressions have been dropped")
        return(design)
    }
    if (any(newnames %in% names(design$variables)))
        stop("Cannot modify pre-existing ", 
            substitute(design)," variables")
    if (any(newnames == "")) {
        # Anche espressioni prive di tag in (...): le rimuovo.
        warning("Untagged input expressions have been dropped")
    }
    full.newnames <- newnames[newnames != ""]
    full.dots <- dots[newnames != ""]
    for (j in seq(along = full.dots)) {
        design$variables[, full.newnames[j]] <- eval(full.dots[[j]], design$variables,
            parent.frame())
    }
    design
}
