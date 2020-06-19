#' BIC for sparse weight based PCA methods 
#'
#' A function that returns the Bayesian information criterion (BIC) given a set of tuning parameters and a sparse weight based PCA method 
#'
#' @param X A data matrix of class 'matrix'
#' @param ncomp An integer specifying the number of components
#' @param FUN A pointer to a function (i.e. the function name with no brackets) that performs sparse weight based PCA. The function should return a list containing a matrix object called "W" that contains the component weights 
#' @param ... specify all the arguments the function in FUN needs 
#' @return The following items in a list \cr
#' \code{BIC} The BIC given the set of tuning parameters \cr
#' \code{nNonZeroCoef} The number of non-zero coefficients in the model \cr
#' @examples
#'  
#' J  <- 10
#' ncomp <- 3
#' X <- matrix(rnorm(100 * J), 100, J)
#' 
#' BICforPCAwithSparseWeights(X = X, ncom = ncomp, FUN = scads, ncomp, ridge = 0, lasso = rep(0.1, ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#' 
#' 
#' #Extended example: select LASSO parameter using BICforPCAwithSparseWeights with scads()
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
#' BIC <- rep(NA, length(lasso))
#' nNonZeroCoef <- rep(NA, length(lasso))
#' for (i in 1:length(lasso)) {
#'     print(i)
#'     res <- BICforPCAwithSparseWeights(X = X, ncomp = ncomp, FUN = scads, ncomp, ridge = 0, lasso = rep(lasso[i], ncomp),
#'                     constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                     itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#'     BIC[i]  <- res$BIC #store BIC for each lasso
#' }
#' 
#' #make a plot of all the BIC's choose the minimal BIC or a model with a small BIC
#' #BIC tends to select very sparse models
#' x <- 1:length(BIC)
#' plot(x, BIC)
#' 
#' #Select the model with the lowest BIC (or alternatively select a model around the model with the lowest BIC)
#' best <- which.min(BIC)
#' 
#' #Do the analysis with the "winning" structure
#' results <- scads(X = X, ncomp = ncomp, ridge = 0, lasso = rep(lasso[best], ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = TRUE , tol = 10^-5)
#' 
#' results$W #inspect results of the estimation
#' dat$P[, 1:ncomp] #inspect data generating model
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

    output <- list(BIC=BIC, nNonZeroCoef = nNonZeroCoef) 
    return(output)

}



