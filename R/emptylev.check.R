emptylev.check <- function (data){
#######################################################
# Checks wheter any factor variable in dataframe data #
# has empty levels, and (in case) drops them.         #
# NOTE: ReGenesees function e.svydesign will ALWAYS   #
#       apply emptylev.check to its data argument:    #
#       thus data bound in design objects will be     #
#       empty-levels free.                            #
# NOTE: The function heavily outperforms droplevels() #
#       especially when data is huge (avoids to call  #
#       factor() unless it is actually needed, and    #
#       handles gc() directly).                       #
#######################################################

    #############################################################
    #  Il parametro need.gc determina se la garbage collection  #
    #  debba, o non debba, essere gestita dal programma.        #
    #  Il valore di soglia per la dimensione di data e'         #
    #  fissata ad 1/10 della memoria massima allocabile.        #
    #############################################################
    need.gc <- FALSE
    if (Sys.info()["sysname"] == "Windows"){
        need.gc <- ((object.size(data)/(1024^2)) > 4096/10)
        # See the NOTE on memory.limit() after 4.2.0 in e.calibrate
    }
    if (need.gc) 
        warning("Input analytic object takes up more than 0.4 GB of allocable memory", 
                immediate. = TRUE)
    gc.here <- function(doit) {
    #################################################
    #  Se doit=TRUE effettua la garbage collection  #
    #  quando viene invocata.                       #
    #################################################
        if (doit)
            gc()
    }

  # On exit, collect garbage always! (good when data is very close to need.gc
  # but doesn't actually exceed the threshold)
  on.exit(gc())

  # Check for factors with empty levels...
  has.emptylev <- function(f) {
                     ( !anyNA(f) & is.factor(f) & ( length(levels(f)) > length(unique(f)) ) ) |
                     (  anyNA(f) & is.factor(f) & ( length(levels(f)) > (length(unique(f)) - 1) ) ) # DEBUG 26/06/2023
                    }
  with.empty <- vapply(data, has.emptylev, NA)
  gc.here(need.gc)
  emptylev.vars <- names(data)[with.empty]
  #... if any, drop empty levels,
  if ( length(emptylev.vars)>0 ) {
     cat("\n# Empty levels found in factors:", paste(emptylev.vars, collapse = ", "))
     data[with.empty] <- lapply(data[with.empty], factor)
     cat("\n# Empty levels have been dropped!\n\n")
     gc.here(need.gc)
     data
    }
   # else, leave data unchanged.
   else data
}
