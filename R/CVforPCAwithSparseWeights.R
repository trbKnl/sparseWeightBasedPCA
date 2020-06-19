#' Cross-validation with the EigenVector method (CV) for sparse weight based PCA methods 
#'
#' A function that returns the mean squared prediction error (MSPE) given a set of tuning parameters and a sparse weight based PCA method 
#'
#' @param X A data matrix of class 'matrix'
#' @param nrFolds The number of folds that should be used. This should be less than \code{nrow(X)}. 
#' @param FUN A pointer to a function (i.e. the function name with no brackets) that performs sparse weight based PCA. It should return a list containing a matrix object called \code{W} that contains the component weights 
#' @param ... specify all the arguments the function in \code{FUN} needs 
#' @return The following items in a list \cr
#' \code{MSPE} The \code{MSPE} given the tuning parameters \cr
#' \code{MSPEstdError} The standard error of the \code{MSPE} \cr
#' \code{nNonZeroCoef} The number of non-zero coefficients in the model 
#' @examples
#'  
#' J  <- 10
#' ncomp <- 3
#' X <- matrix(rnorm(100 * J), 100, J)
#' 
#' CVforPCAwithSparseWeights(X = X, nrFolds = 10, FUN = scads, ncomp, ridge = 0, lasso = rep(0.1, ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#' 
#' 
#' #Extended example: select LASSO parameter using CVforPCAwithSparseWeights and the 1 standard error rule with scads()
#' #create sammple data
#' ncomp <- 3 
#' J <- 30
#' comdis <- matrix(1, J, ncomp)
#' comdis <- sparsify(comdis, 0.7) #set 70 percent of the 1's to zero
#' variances <- makeVariance(varianceOfComps = c(100, 80, 90), J = J, error = 0.05) #create realistic eigenvalues
#' dat <- makeDat(n = 100, comdis = comdis, variances = variances)
#' X <- dat$X
#' 
#' 
#' #Use cross-validation to look for the data generating structure 
#' lasso <- seq(0, 1, length.out = 500)
#' MSPE <- rep(NA, length(lasso))
#' MSPEstdError <- rep(NA, length(lasso))
#' nNonZeroCoef <- rep(NA, length(lasso))
#' for (i in 1:length(lasso)) {
#'     print(i)
#'     res <- CVforPCAwithSparseWeights(X = X, nrFolds = 10, FUN = scads, ncomp, ridge = 0, lasso = rep(lasso[i], ncomp),
#'                     constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                     itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#'     MSPE[i]  <- res$MSPE #store MSPE for each lasso
#'     MSPEstdError[i]  <- res$MSPEstdError #store standard error for each lasso
#'     nNonZeroCoef[i] <- res$nNonZeroCoef #store the number of non-zero coefficients for each model
#' }
#' 
#' #make a plot of all the MSPE's
#' x <- 1:length(MSPE)
#' plot(x, MSPE)
#' #add errorbars
#' arrows(x, MSPE - MSPEstdError, x , MSPE + MSPEstdError, length = 0.05, angle = 90, code = 3)
#' 
#' #Select all models within one standard error of the best model
#' eligibleModels <- MSPE < MSPE[which.min(MSPE)] + MSPEstdError[which.min(MSPE)]
#' 
#' #Selected from those models the models with the lowest amount of non-zero coeffients
#' best <- which.min(nNonZeroCoef[eligibleModels])
#' 
#' #Do the analysis with the "winning" lasso
#' results <- scads(X = X, ncomp = ncomp, ridge = 0, lasso = rep(lasso[best], ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = TRUE , tol = 10^-5)
#' 
#' results$W #inspect results of the estimation
#' dat$P[, 1:ncomp] #inspect data generating model
#' 
#' @export
CVforPCAwithSparseWeights <- function(X, nrFolds, FUN, ...)  {

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




