`ez.ranger` <-
function (design, formula, population)
###############################################################
# La funzione calcola, per ogni singoli sub-task del processo #
# di calibrazione complessivo, il range dei rapporti fra      #
# totali noti e corrispondenti stime dirette.                 #
###############################################################
{
    mm <- model.matrix(formula, model.frame(formula, data = design$variables))
    ww <- design$dir.weights
    sample.total <- colSums(mm * ww)
    if (length(sample.total) != length(population)) 
        stop("Population and sample totals are not the same length.")
    if (any(sample.total == 0)) {
#       zz <- (population == 0) & (apply(mm, 2, function(x) all(x == 0))) # OLD: It's much slower than NEW
        zz <- (population == 0) & (colSums(abs(mm)) == 0)
        mm <- mm[, !zz, drop = FALSE]
        population <- population[!zz]
        sample.total <- sample.total[!zz]
    }
    smallest.int <- range(population/sample.total)
    if ( any(population < 0) || any(sample.total < 0) ) {
         # If HT estimates or population benchmarks are *negative*, then
         # suggested interval is not expected to be reliable.
         # Thus must track those cases...
         attr(smallest.int, "negatives") <- TRUE
     }
    smallest.int
}
