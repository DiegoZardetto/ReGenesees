FUN.ACCUMULATOR <- function(FUN) {
  acc.name <- paste("ACCUMULATOR.", FUN, sep="")
  if (!exists(acc.name, envir = .GlobalEnv)){
     eval(parse(text=paste("`", acc.name, "` <<- 0", sep="")))
    }
  else {
     eval(parse(text=paste("`", acc.name, "` <<- `", acc.name, "` + 1", sep="")))
    }
}

# TO BE PASTED IN THE FUN BODY
# FUN.ACCUMULATOR( FUN = deparse(match.call()[[1]]) )
