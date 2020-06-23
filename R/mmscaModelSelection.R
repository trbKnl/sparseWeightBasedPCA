#Function to get all combinations of tuning parameters given sequences of ridge, lasso etc, and the number of components ncomp 
combOfTuningParams <- function(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq, ncompSeq = NULL, printProgress = TRUE) {

    if (is.null(ridgeSeq) || is.null(lassoSeq) || 
       is.null(grouplassoSeq) || is.null(elitistlassoSeq)) {stop("A sequence for either the ridge, lasso, grouplasso,
   or elitistlasso is missing. If a penalty is not wanted specifiy the sequence to be a numeric value of zero")}


    if (is.null(ncompSeq)) {stop("Give the number of components: an integer or an integer sequence")}

    outlist  <- vector("list", length = length(ncompSeq))
    lambdaList <- list(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq)
    totalNumberOfCombs <- prod(unlist(lapply(lambdaList, length)))

    if (printProgress) {
        cat(paste("Total number of models: ", totalNumberOfCombs * length(ncompSeq), "\n", sep=""))
    }

    combmat <- as.list(expand.grid(lambdaList))
    names(combmat) <- NULL

    for (i in 1:length(ncompSeq)) {
        outlist[[i]] <- lapply(combmat, function(x){out <- matrix(NA, length(x), ncompSeq[i]); out[] <- x; return(out)})
    }
    return(outlist)
    
}

#' Model selection for MMSCA  
#'
#' A function that performs model selection, for the regularizers and the number of components for \code{mmsca()} 
#'
#' @param X A data matrix of class \code{matrix}
#' @param ridgeSeq A range of values for the ridge penalty that need to be examined. Specify a zero if the tuning parameter is not wanted.
#' @param lassoSeq A range of values for the lasso penalty that need to be examined. Specify a zero if the tuning parameter is not wanted.
#' @param grouplassoSeq A range of values for the group lasso penalty that need to be examined. Specify a zero if the tuning parameter is not wanted.
#' @param elitistlassoSeq A range of values for the elitist lasso penalty that need to be examined. Specify a zero if the tuning parameter is not wanted.
#' @param ncompSeq A range of integers for the number of components that need to be examined. 
#' @param tuningMethod A string indicating which model selection method should be used. "BIC" enables the Bayesian information criterion, "IS" enables the index of sparseness. "CV" enables cross-validation (CV) with the EigenVector method, if CV is used, the number of folds nrFolds needs to be chosen. The number of folds should be an integer less than \code{nrow(X)}. The data are then split in equal sized chunks if order of appearance. 
#' @param groups A vector specifying which columns of X belong to what block. Example: \code{c(10, 100, 1000)}. The first 10 variables belong to the first block, the 100 variables after that belong to the second block etc. 
#' @param nrFold An integer that specify the number of folds that Cross-validation should use if tuningmethod == "CV", the number of folds needs to be lower then \code{nrow(X)}. 
#' @param itr The maximum number of iterations (a positive integer) 
#' @param nStarts The number of random starts the analysis should perform. The first start will be a warm start. You can not give custom starting values. 
#' @param tol The convergence is determined by comparing the loss function value after each iteration, if the difference is smaller than \code{tol}, the analysis is converged. Default value is \code{10e-8}
#' @param coorDes A boolean with the default \code{FALSE}. If coorDes is \code{FALSE} the estimation of the majorizing function to estimate the component weights W conditional on the loadings P will be found using matrix inverses which can be slow. If set to \code{TRUE} the marjozing function will be optimized (or partially optimized) using coordinate descent, in some cases coordinate descent will be faster
#' @param coorDesItr An integer specifying the maximum number of iterations for the coordinate descent algorithm, the default is set to 1. You do not have to run this algorithm until convergence before alternating back to the estimation of the loadings. The tolerance for this algorithm is hardcoded and set to \code{10^-8}.
#' @param printProgress A boolean: \code{TRUE} will print the progress of the model selection
#' @return A list containing: \cr
#' \code{results} A list with \code{ncomp} elements each containing the following items in a list \cr
#' \itemize{
#'  \item{"BIC, IS or MSPE"}{ The index chosen in tuningMethod for all combinations of ridge, lasso, grouplasso and elististlasso}
#'  \item{"bestBIC, bestIS, bestMSPE or bestMSPE1stdErrorRule"}{ The best index according to the chosen tuning method}
#'  \item{"nNonZeroCoef"}{ The number of non zero weights in the best model}
#'  \item{"ridge"}{ The value of the ridge penalty corresponding to the best model}
#'  \item{"lasso"}{ The value of the lasso penalty corresponding to the best model}
#'  \item{"grouplasso"}{ The value of the group lasso penalty corresponding to the best model}
#'  \item{"elististlasso"}{ The value of the elitist lasso penalty corresponding to the best model}
#'  \item{"ncomp"}{ The number of component that was used for these items}
#'  \item{"ridge1stdErrorRule"}{ In case tuningMethod == "CV", the value of the ridge penalty according to the 1 standard error rule: the most sparse model within one standard error of the model with the lowest MSPE}
#'  \item{"lasso1stdErrorRule"}{ In case tuningMethod == "CV", the value of the lasso penalty according to the 1 standard error rule: the most sparse model within one standard error of the model with the lowest MSPE}
#'  \item{"grouplasso1stdErrorRule"}{ In case tuningMethod == "CV", the value of the group lasso penalty according to the 1 standard error rule: the most sparse model within one standard error of the model with the lowest MSPE}
#'  \item{"elitistlasso1stdErrorRule"}{ In case tuningMethod == "CV", the value of the elitist lasso penalty according to the 1 standard error rule: the most sparse model within one standard error of the model with the lowest MSPE}
#'  \item{"ridge1stdErrorRule"}{ In case tuningMethod == "CV", the value of the ridge according to the 1 standard error rule: the most sparse model within one standard error of the model with the lowest MSPE}
#'  }
#' \code{bestNcomp} The number of component with the best value for the chosen tuning index \cr
#' @examples
#'  
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' 
#' out <- mmscaModelSelection(X, 
#'             ridgeSeq = seq(0, 1, by = 0.1), 
#'             lassoSeq = 0:100, 
#'             grouplassoSeq = 0,
#'             elitistlassoSeq = 0, 
#'             ncompSeq = 1:3, 
#'             tuningMethod = "CV", 
#'             groups = ncol(X), 
#'             nrFolds = 10, 
#'             itr = 100000, 
#'             nStart = 1, 
#'             coorDes = FALSE, 
#'             coorDesItr = 100, 
#'             printProgress = TRUE)
#'
#' #Inspect the results of the model selection for the optimal number of components according to the tuning method
#' out$results[[out$bestNcomp]]
#' 
#' @export
mmscaModelSelection <- function(X, ridgeSeq, lassoSeq, grouplassoSeq, 
                                elitistlassoSeq, ncompSeq,
                                tuningMethod = "BIC",
                                groups, 
                                nrFolds = NULL, 
                                itr = 10e5, nStart = 1,
                                tol = 10e-8,
                                coorDes = TRUE, coorDesItr = 100,
                                printProgress = TRUE) {

    combs <- combOfTuningParams(ridgeSeq, lassoSeq, grouplassoSeq,
                       elitistlassoSeq, ncompSeq, printProgress)


    #Function is written so that new model selection procedures can be added/changed in the if-else blocks
    #Not pretty, but easy to understand and maintain
    ################## BIC ################## 
    if (tuningMethod == "BIC") {

        I <- nrow(X)
        p <- ncol(X)
        oldBest <- Inf

        tuningIndexGivenComp <- vector("list", length(ncompSeq))

        for (j in 1:length(ncompSeq)) {

            ncomp <- ncompSeq[j]
            total_coef <- sum(groups) * ncomp

            names(tuningIndexGivenComp)[j] <- paste("ncomp_", ncomp, sep ="")

            BIC <- rep(NA, nrow(combs[[1]][[1]]))
            nNonZeroCoef <- rep(NA, nrow(combs[[1]][[1]]))

            VarSelect0 <- mmsca(X, 
                               ridge = rep(0, ncomp),
                               lasso = rep(0, ncomp),
                               constraints = matrix(1, p, ncomp),
                               grouplasso = rep(0, ncomp),
                               elitistlasso = rep(0, ncomp),
                               groups = groups, 
                               ncomp = ncomp, 
                               nStart = nStart,
                               tol = tol,
                               itr = 100000, 
                               printLoss = FALSE,
                               Wstart = matrix(0, p, ncomp), 
                               coorDes = coorDes, 
                               coorDesItr = coorDesItr)

            T_hat0 <- X %*% VarSelect0$W
            P_hat0 <- VarSelect0$P
            residual_variance0 <- sum((X - T_hat0 %*% t(P_hat0))^2)

            for (i in 1:nrow(combs[[1]][[1]])) {

                if (printProgress == TRUE) {
                    cat(paste("Number of components: ", ncomp,
                                ", model number: ", i, "\n", sep = ""))
                }

                VarSelect <- mmsca(X, 
                                   ridge = combs[[j]][[1]][i, ],
                                   lasso = combs[[j]][[2]][i, ],
                                   constraints = matrix(1, p, ncomp),
                                   grouplasso = combs[[j]][[3]][i, ],
                                   elitistlasso = combs[[j]][[4]][i, ],
                                   groups = groups, 
                                   ncomp = ncomp, 
                                   nStart = nStart,
                                   tol = tol,
                                   itr = 100000, 
                                   printLoss = FALSE,
                                   Wstart = matrix(0, p, ncomp), 
                                   coorDes = coorDes, 
                                   coorDesItr = coorDesItr)

                P_hat <- VarSelect$P
                W_hat <- VarSelect$W
                T_hat <- X %*% W_hat 
                nonZeroCoefInW_hat <- sum(W_hat != 0)

                residual_variance <- sum((X - T_hat %*% t(P_hat))^2)
                nNonZeroCoef[i] <- nonZeroCoefInW_hat

                BIC[i]  <- (residual_variance / residual_variance0) +
                    (nonZeroCoefInW_hat * (log(I) / I))
            }
            bestTuningSeqsBIC <- lapply(combs[[1]],
                                        function(x){return(x[which.min(BIC), ])}) 
            bestBIC <- BIC[which.min(BIC)]

            if (bestBIC < oldBest) {
                bestNcomp <- ncomp
                oldBest <- bestBIC
            }

            ridge <- bestTuningSeqsBIC[[1]]
            lasso <- bestTuningSeqsBIC[[2]]
            grouplasso <- bestTuningSeqsBIC[[3]]
            elitistlasso <- bestTuningSeqsBIC[[4]]

            output <- list(BIC = BIC, bestBIC = bestBIC, 
                           nNonZeroCoef = nNonZeroCoef[which.min(BIC)],
                           ridge = ridge, lasso = lasso, grouplasso = grouplasso,
                           elitistlasso = elitistlasso, ncomp = ncomp) 
            tuningIndexGivenComp[[j]] <- output 
        }
        return(list(results = tuningIndexGivenComp, bestNcomp = bestNcomp))

    ################## Index of Sparseness ################## 
    } else if (tuningMethod == "IS") {


        I <- nrow(X)
        p <- ncol(X)
        oldBest <- -Inf

        tuningIndexGivenComp <- vector("list", length(ncompSeq))

        for (j in 1:length(ncompSeq)) {

            ncomp <- ncompSeq[j]
            total_coef <- sum(groups) * ncomp

            names(tuningIndexGivenComp)[j] <- paste("ncomp_", ncomp, sep ="")

            IS <- rep(NA, nrow(combs[[1]][[1]]))
            nNonZeroCoef <- rep(NA, nrow(combs[[1]][[1]]))

            VarSelect0 <- mmsca(X, 
                               ridge = rep(0, ncomp),
                               lasso = rep(0, ncomp),
                               constraints = matrix(1, p, ncomp),
                               grouplasso = rep(0, ncomp),
                               elitistlasso = rep(0, ncomp),
                               groups = groups, 
                               ncomp = ncomp, 
                               tol = tol,
                               nStart = nStart,
                               itr = 100000, 
                               printLoss = FALSE,
                               Wstart = matrix(0, p, ncomp), 
                               coorDes = coorDes, 
                               coorDesItr = coorDesItr)

            T_hat0 <- X %*% VarSelect0$W
            P_hat0 <- VarSelect0$P

            V_s <- sum((T_hat0%*%t(P_hat0))^2) 
            V_oo <- sum(X^2)  

            for (i in 1:nrow(combs[[1]][[1]])) {

                if (printProgress == TRUE) {
                    cat(paste("Number of components: ", ncomp,
                                ", model number: ", i, "\n", sep = ""))
                }

                VarSelect <- mmsca(X, 
                                   ridge = combs[[j]][[1]][i, ],
                                   lasso = combs[[j]][[2]][i, ],
                                   constraints = matrix(1, p, ncomp),
                                   grouplasso = combs[[j]][[3]][i, ],
                                   elitistlasso = combs[[j]][[4]][i, ],
                                   groups = groups, 
                                   tol = tol,
                                   ncomp = ncomp, 
                                   nStart = nStart,
                                   itr = 100000, 
                                   printLoss = FALSE,
                                   Wstart = matrix(0, p, ncomp), 
                                   coorDes = coorDes, 
                                   coorDesItr = coorDesItr)

                P_hat <- VarSelect$P
                W_hat <- VarSelect$W
                T_hat <- X %*% W_hat 
                nonZeroCoefInW_hat <- sum(W_hat != 0)

                nNonZeroCoef[i] <- nonZeroCoefInW_hat

                V_a <- sum((T_hat %*% t(P_hat))^2)  
                IS[i] <- V_a * V_s / V_oo^2 * (sum(W_hat == 0) / total_coef)
            }
            bestTuningSeqsIS <- lapply(combs[[1]],
                                        function(x){return(x[which.max(IS), ])}) 
            bestIS <- IS[which.max(IS)]
            if (bestIS > oldBest) {
                bestNcomp <- ncomp
                oldBest <- bestIS
            }

            ridge <- bestTuningSeqsIS[[1]]
            lasso <- bestTuningSeqsIS[[2]]
            grouplasso <- bestTuningSeqsIS[[3]]
            elitistlasso <- bestTuningSeqsIS[[4]]

            output <- list(IS = IS, bestIS = bestIS, nNonZeroCoef = nNonZeroCoef[which.max(IS)], 
                           ridge = ridge, lasso = lasso, grouplasso = grouplasso,
                           elitistlasso = elitistlasso, ncomp = ncomp) 
            tuningIndexGivenComp[[j]] <- output 
        }
        return(list(results=tuningIndexGivenComp, bestNcomp = bestNcomp))

    ################## CV Eigenvector method ################## 
    } else if (tuningMethod == "CV") {

        I <- nrow(X)
        p <- ncol(X)
        oldBest <- Inf

        if (is.null(nrFolds)) { stop("Specify the number of folds: nrFolds") }
        if (nrFolds > I) { stop("nrFolds must be less than nrow(X)") }

        tuningIndexGivenComp <- vector("list", length(ncompSeq))

        for (j in 1:length(ncompSeq)) {

            ncomp <- ncompSeq[j]
            total_coef <- sum(groups) * ncomp

            names(tuningIndexGivenComp)[j] <- paste("ncomp_",
                                                    ncomp, sep ="")

            MSPE <- rep(NA, nrow(combs[[1]][[1]]))
            MSPEstdError <- rep(NA, nrow(combs[[1]][[1]]))
            nNonZeroCoef <- rep(NA, nrow(combs[[1]][[1]]))

            folds <- rep_len(1:nrFolds, nrow(X))

            for (i in 1:nrow(combs[[1]][[1]])) {

                if (printProgress == TRUE) {
                    cat(paste("Number of components: ", ncomp,
                                ", model number: ", i, "\n", sep = ""))
                }

                cvError  <- matrix(NA, nrow(X), ncol(X))
                MSEkthFold <- rep(NA, nrFolds)

                for(a in 1:nrFolds){

                    fold <- which(folds == a)
                    trainDat <- X[-fold, ]
                    testDat <- X[fold, , drop = FALSE]

                    res <- mmsca(trainDat, 
                                       ridge = combs[[j]][[1]][i, ],
                                       lasso = combs[[j]][[2]][i, ],
                                       constraints = matrix(1, p, ncomp),
                                       grouplasso = combs[[j]][[3]][i, ],
                                       elitistlasso = combs[[j]][[4]][i, ],
                                       groups = groups, 
                                       tol = tol,
                                       ncomp = ncomp, 
                                       nStart = nStart,
                                       itr = 100000, 
                                       printLoss = FALSE,
                                       Wstart = matrix(0, p, ncomp), 
                                       coorDes = coorDes, 
                                       coorDesItr = coorDesItr)

                    pred <- matrix(NA, nrow(testDat), ncol(testDat))
                    for (b in 1:ncol(X)) {
                        TMinusVariableJ <- testDat[, -b] %*% res$W[-b, ]
                        pred[, b] <- TMinusVariableJ %*% res$P[b, ] 
                    }
                    cvError[fold, ] <- (testDat - pred)^2
                    MSEkthFold[a]  <-  mean(cvError[fold, ]) 
                }

                MSPE[i] <- mean(MSEkthFold)
                MSPEstdError[i] <- sd(MSEkthFold) / sqrt(nrFolds)

                # Determine the number of non-zero weights 
                res <- mmsca(X, 
                                   ridge = combs[[j]][[1]][i, ],
                                   lasso = combs[[j]][[2]][i, ],
                                   constraints = matrix(1, p, ncomp),
                                   grouplasso = combs[[j]][[3]][i, ],
                                   elitistlasso = combs[[j]][[4]][i, ],
                                   groups = groups, 
                                   tol = tol,
                                   ncomp = ncomp, 
                                   nStart = nStart,
                                   itr = 100000, 
                                   printLoss = FALSE,
                                   Wstart = matrix(0, p, ncomp), 
                                   coorDes = coorDes, 
                                   coorDesItr = coorDesItr)

                nNonZeroCoef[i] <- sum(res$W != 0)
            }

            bestTuningSeqsCV <- lapply(combs[[1]],
                                        function(x){return(x[which.min(MSPE), ])}) 
            bestModel <- which.min(MSPE)
            bestMSPE <- MSPE[bestModel]

            if (bestMSPE < oldBest) {
                bestNcomp <- ncomp
                oldBest <- bestMSPE
            }

            MSPE1stdErrorRule <- MSPE[MSPE < (MSPE[bestModel] +
                                           MSPEstdError[bestModel])]
            bestModel1stdErrorRule <- which(MSPE1stdErrorRule[length(MSPE1stdErrorRule)]
                                            == MSPE)
            bestMSPE1stdErrorRule <- MSPE[bestModel1stdErrorRule]

            bestTuningSeqsCV1stdErrorRule <- lapply(combs[[1]],
                                                    function(x){return(x[bestModel1stdErrorRule, ])})

            ridge <- bestTuningSeqsCV[[1]]
            lasso <- bestTuningSeqsCV[[2]]
            grouplasso <- bestTuningSeqsCV[[3]]
            elitistlasso <- bestTuningSeqsCV[[4]]

            ridge1stdErrorRule <- bestTuningSeqsCV1stdErrorRule[[1]]
            lasso1stdErrorRule <- bestTuningSeqsCV1stdErrorRule[[2]]
            grouplasso1stdErrorRule <- bestTuningSeqsCV1stdErrorRule[[3]]
            elitistlasso1stdErrorRule <- bestTuningSeqsCV1stdErrorRule[[4]]

            output <- list(MSPE  =  MSPE, MSPEstdError = MSPEstdError, 
                           nNonZeroCoef = nNonZeroCoef[which.min(MSPE)],
                           bestMSPE = bestMSPE, ridge = ridge, lasso = lasso,
                           grouplasso = grouplasso, elitistlasso = elitistlasso,
                           ncomp = ncomp, bestMSPE1stdErrorRule = bestMSPE1stdErrorRule,
                           ridge1stdErrorRule = ridge1stdErrorRule,
                           lasso1stdErrorRule = lasso1stdErrorRule, 
                           grouplasso1stdErrorRule = grouplasso1stdErrorRule, 
                           elitistlasso1stdErrorRule = elitistlasso1stdErrorRule) 
            tuningIndexGivenComp[[j]] <- output 
        }
        return(list(results = tuningIndexGivenComp, bestNcomp = bestNcomp))
    }
}

