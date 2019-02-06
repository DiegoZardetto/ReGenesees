`NA.estvars` <- 
function(design, estvars, na.rm = FALSE, draconian = FALSE)
##################################################
# Check of missing values in interest variables, #
# plus reaction policy.                          #
################################################## 
{
    data <- design$variables
    estvars <- unique(estvars)
    absent.vars <- estvars[!(estvars %in% names(data))]
    if (length(absent.vars) > 0) 
        stop("Variables not found: ", paste(absent.vars, collapse = ", "))
    has.na <- sapply(estvars, function(var) any(is.na(data[, 
        var])))
    # Se sono presenti NA...
    if (any(has.na)){
        # Draconiano: i.e. in questa funzione NON posso gestire gli NA:
        if (isTRUE(draconian)) {
            stop("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                 "\n", "(this function can't handle NAs)\n")
            }
        # Una variabile di interesse: ha senso na.rm=TRUE
        if (length(estvars)==1){
            if (!na.rm){
                # Disegno non calibrato: gli NA non causano errori anche se na.rm=FALSE:
                if (!is.calibrated(design)) {
                    warning("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                            "\n(you may want to specify na.rm=TRUE)")
                }
                # Disegno calibrato: gli NA se na.rm=FALSE causano errori:
                else{
                    stop("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                         "\nyou need to specify na.rm=TRUE")
                }
            }
        }
        # Variabili di interesse multiple: non ha senso na.rm=TRUE
        else{
            # Disegno calibrato: gli NA causano errori
            if (is.calibrated(design)) {
                # Non posso gestire variabili multiple (neanche con na.rm=TRUE)
                stop("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                     "\n", if (na.rm) "(na.rm = TRUE can't handle multiple interest variables)\n")
            }
            # Disegno non calibrato: posso gestire solo na.rm=FALSE 
            else {
                if (na.rm) {
                    stop("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                     "\n(na.rm = TRUE can't handle multiple interest variables)\n")
                }
                else {
                    warning("Missing values in: ", paste(estvars[has.na], collapse = ", "), 
                            "\n")
                }
            }
        }
    }
}
