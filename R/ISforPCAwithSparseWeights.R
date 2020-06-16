#' Index of sparseness (IS) for sparse weight based PCA methods 
#'
#' A function that returns the IS given a set of tuning parameters and a sparse weight based PCA method 
#'
#' @param X A data matrix of class 'matrix'
#' @param ncomp An integer specifying the number of components
#' @param FUN A pointer to a function (i.e. the function name with no brackets) that performs sparse weight based PCA. it should return a list containing a matrix object called "W" that contains the component weights 
#' @param ... specify all the arguments the function in FUN needs 
#' @return The following items in a list \cr
#' \code{IS} The IS given the tuning parameters \cr
#' \code{nNonZeroCoef} The number of non-zero coefficients in the model \cr
#' @examples
#'  
#' 
#' J  <- 10
#' ncomp <- 3
#' X <- matrix(rnorm(100 * J), 100, J)
#' 
#' ISforPCAwithSparseWeights(X = X, ncom = ncomp, FUN = scads, ncomp, ridge = 0, lasso = rep(0.1, ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#' 
#' 
#' #Extended example: select LASSO parameter using ISforPCAwithSparseWeights with scads()
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
#' IS <- rep(NA, length(lasso))
#' nNonZeroCoef <- rep(NA, length(lasso))
#' for (i in 1:length(lasso)) {
#'     print(i)
#'     res <- ISforPCAwithSparseWeights(X = X, ncomp = ncomp, FUN = scads, ncomp, ridge = 0, lasso = rep(lasso[i], ncomp),
#'                     constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                     itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)
#'     IS[i]  <- res$IS #store IS for each lasso
#' }
#' 
#' #make a plot of all the IS's choose the maximal IS or a model that is less sparse around the maximal IS
#' #IS tends to select very sparse models
#' x <- 1:length(IS)
#' plot(x, IS)
#' 
#' #select the maximal IS (alternative select a model that is around the maximum, but results in less non-zero weights)
#' best <- which.max(IS)
#' best
#' 
#' #Do the analysis with the "winning" structure
#' results <- scads(X = X, ncomp = ncomp, ridge = 0, lasso = rep(lasso[best], ncomp),
#'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
#'                 itr = 100000, nStarts = 1, printLoss = TRUE , tol = 10^-8)
#' 
#' results$W #inspect results of the estimation
#' dat$P[, 1:ncomp] #inspect data generating model
#' 
#' @export
ISforPCAwithSparseWeights <- function(X, ncomp, FUN, ...)  {

    I <- nrow(X)
    p <- ncol(X)
    totalCoef <- p*ncomp
    W <- svd(X)$v[, 1:ncomp]

    T_hat0 <- X %*% W
    P_hat0 <- W

    V_s <- sum((T_hat0%*%t(P_hat0))^2) 
    V_oo <- sum(X^2)  

    VarSelect <- FUN(X, ...) 

    P_hat <- VarSelect$P
    W_hat <- VarSelect$W
    T_hat <- X %*% W_hat 
    nNonZeroCoef <- sum(W_hat != 0)

    V_a <- sum((T_hat %*% t(P_hat))^2)  
    IS <- V_a * V_s / V_oo^2 * (sum(W_hat == 0) / totalCoef)

    output <- list(IS = IS, nNonZeroCoef = nNonZeroCoef)
    return(output)
}




