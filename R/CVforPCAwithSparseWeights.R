#' Cross-validation with the EigenVector method (CV) for sparse weight based PCA methods 
#'
#' A function that returns the mean squared prediction error (MSPE) given a set of tuning parameters and a sparse weight based PCA method 
#'
#' @param X A data matrix of class 'matrix'
#' @param ncomp An integer specifying the number of components
#' @param FUN A pointer to a function (i.e. the function name with no brackets) that performs sparse weight based PCA. it should return a list containing a matrix object called "W" that contains the component weights 
#' @param ... specify all the arguments the function in FUN needs 
#' @return The following items in a list \cr
#' \code{MSPE} The MSPE given the tuning parameters 
#' \code{MSPEstdError} The standard error of the MSPE 
#' \code{nNonZeroCoef} The number of non-zero coefficients in the model 
#' @examples
#'  
#' @export
CVforPCAwithSparseWeights <- function(X, ncomp, nrFolds, FUN, ...)  {

    folds <- rep_len(1:nrFolds, nrow(X))
    cvError  <- matrix(NA, nrow(X), ncol(X))
    MSEkthFold <- rep(NA, nrFolds)

    for (a in 1:nrFolds){

        fold <- which(folds == a)
        trainDat <- X[-fold, ]
        testDat <- X[fold, , drop = FALSE]

        res <- FUN(trainDat, ...) 
        pred <- matrix(NA, nrow(testDat), ncol(testDat))

        for (b in 1:ncol(X)) {
            TMinusVariableJ <- testDat[, -b] %*% res$W[-b, ]
            pred[, b] <- TMinusVariableJ %*% res$P[b, ] 
        }

        cvError[fold, ] <- (testDat - pred)^2
        MSEkthFold[a]  <-  mean(cvError[fold, ]) 
    }

    MSPE <- mean(MSEkthFold)
    MSPEstdError <- sd(MSEkthFold) / sqrt(nrFolds)

    # Determine the number of non-zero weights 
    res <- FUN(X, ...) 
    nNonZeroCoef <- sum(res$W != 0)

    output <- list(MSPE=MSPE, MSPEstdError=MSPEstdError, 
                   nNonZeroCoef=nNonZeroCoef)
    return(output)
}







