#' BIC for sparse weight based PCA methods 
#'
#' A function that returns the baysian information criterion (BIC) given a set of tuning parameters and a sparse weight based PCA method 
#'
#' @param X A data matrix of class 'matrix'
#' @param ncomp An integer specifying the number of components
#' @param FUN A pointer to a function (i.e. the function name with no brackets) that performs sparse weight based PCA. it should return a list containing a matrix object called "W" that contains the component weights 
#' @param ... specify all the arguments the function in FUN needs 
#' @return The following items in a list \cr
#' \code{BIC} The BIC given the set of tuning parameters 
#' \code{nNonZeroCoef} The number of non-zero coefficients in the model 
#' @examples
#'  
#' @export
BICforPCAwithSparseWeights <- function(X, ncomp, FUN, ...)  {

    I <- nrow(X)
    p <- ncol(X)

    W <- svd(X)$v[, 1:ncomp]

    T_hat0 <- X %*% W
    P_hat0 <- W
    residual_variance0 <- sum((X - T_hat0 %*% t(P_hat0))^2)

    VarSelect <- FUN(X, ...) 

    P_hat <- VarSelect$P
    W_hat <- VarSelect$W
    T_hat <- X %*% W_hat 
    nonZeroCoefInW_hat <- sum(W_hat != 0)

    residual_variance <- sum((X - T_hat %*% t(P_hat))^2)
    nNonZeroCoef <- nonZeroCoefInW_hat

    BIC  <- (residual_variance / residual_variance0) +
        (nonZeroCoefInW_hat * (log(I) / I))

    output <- list(BIC=BIC, nNonZeroCoef) 
    return(output)

}



