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
    smallest.int
}
