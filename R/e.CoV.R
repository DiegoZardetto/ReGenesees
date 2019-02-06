#################################################
# Design Covariance and Correlation Estimation. #
#################################################

CoV <- function(design, expr1, expr2,
                by = NULL, na.rm = FALSE){
# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Compute basic Estimators (this will also perform the fundamental checks,
# e.g. length 1 expressions...)
E1 <- svystatL(design, expr1, by = by, na.rm = na.rm)
E2 <- svystatL(design, expr2, by = by, na.rm = na.rm)

# Get the strings representing the Complex Estimators (could also take
# those from SEs' names...)
s.e1 <- deparse(expr1[[1]])
s.e2 <- deparse(expr2[[1]])

# Build the Estimator of the Sum
expr12 <- parse( text = paste(s.e1, "+", s.e2) )
E12 <-svystatL(design, expr12, by = by, na.rm = na.rm)

# Compute basic SEs
SE1 <- SE(E1)
SE2 <- SE(E2)

# Compute SE of the Estimator of the Sum
SE12 <- SE(E12)

# Compute Covariance
CoV <- ( SE12^2 - SE1^2 - SE2^2 ) / 2

# Give Covariance Estimator a name
CoV.name <- paste("CoV(", s.e1, ", ", s.e2,")", sep="")

# Check if by != NULL and (in case) get domain columns from E12,
# in order to build output data frame
if (!is.null(by)) {
     by.cols <- attr(E12, "svyby")$margins
     CoV.df <- E12[, by.cols, drop=FALSE]
     CoV <- as.data.frame(CoV)
     colnames(CoV) <- CoV.name
     .CoV <- cbind(CoV.df, CoV)
     class(.CoV) <- c("CoV.by", "CoV", class(.CoV))
    }
else {
     colnames(CoV) <- CoV.name
     .CoV <- as.data.frame(CoV)
     class(.CoV) <- c("CoV", class(.CoV))
    }

attr(.CoV,"design") <- design.expr
.CoV
}


Corr <- function(design, expr1, expr2,
                 by = NULL, na.rm = FALSE){
# First verify if the function has been called inside another function:
# this is needed to correctly manage metadata when e.g. the caller is a
# GUI stratum
directly <- !( length(sys.calls()) > 1 )
design.expr <- if (directly) substitute(design)

# Compute basic Estimators (this will also perform the fundamental checks,
# e.g. length 1 expressions...)
E1 <- svystatL(design, expr1, by = by, na.rm = na.rm)
E2 <- svystatL(design, expr2, by = by, na.rm = na.rm)

# Get the strings representing the Complex Estimators (could also take
# those from SEs' names...)
s.e1 <- deparse(expr1[[1]])
s.e2 <- deparse(expr2[[1]])

# Build the Estimator of the Sum
expr12 <- parse( text = paste(s.e1, "+", s.e2) )
E12 <-svystatL(design, expr12, by = by, na.rm = na.rm)

# Compute basic SEs
SE1 <- SE(E1)
SE2 <- SE(E2)

# Compute SE of the Estimator of the Sum
SE12 <- SE(E12)

# Compute Correlation
CoV <- ( SE12^2 - SE1^2 - SE2^2 ) / 2
Corr <- CoV / (SE1*SE2)

# Give Covariance Estimator a name
Corr.name <- paste("Corr(", s.e1, ", ", s.e2,")", sep="")

# Check if by != NULL and (in case) get domain columns from E12,
# in order to build output data frame
if (!is.null(by)) {
     by.cols <- attr(E12, "svyby")$margins
     Corr.df <- E12[, by.cols, drop=FALSE]
     Corr <- as.data.frame(Corr)
     colnames(Corr) <- Corr.name
     .Corr <- cbind(Corr.df, Corr)
     class(.Corr) <- c("Corr.by", "Corr", class(.Corr))
    }
else {
     colnames(Corr) <- Corr.name
     .Corr <- as.data.frame(Corr)
     class(.Corr) <- c("Corr", class(.Corr))
    }

attr(.Corr,"design") <- design.expr
.Corr
}
