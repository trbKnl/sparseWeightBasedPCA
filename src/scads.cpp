#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


//' This function performs sparse PCA with constraints on the component weights and/or ridge and lasso regularization.
//' 
//' @param X A data matrix of class 'matrix'
//' @param ncomp The number of components to estimate (an integer)
//' @param ridge A numeric value containing the ridge parameter for the component weight matrix W
//' @param lasso A vector containing a ridge parameter for each column of W seperately, to set the same lasso penalty for the component weights W, specify: lasso = rep(value, ncomp)
//' @param constraints A matrix of the same dimensions as the component weights matrix W (ncol(X) x ncomp). A zero entry corresponds in constraints corresponds to an element in the same location in W that needs to be constraint to zero. A non-zero entry corresponds to an element in the same location in W that needs to be estimated.
//' @param itr The maximum number of iterations (an integer)
//' @param Wstart A matrix of ncomp columns and nrow(X) rows with starting values for the component weight matrix W, if Wstart only contains zeros, a warm start is used: the first ncomp right singular vectors of X
//' @param tol The convergence is determined by comparing the loss function value after each iteration, if the difference is smaller than tol, the analysis is converged. The default value is 10e-8.
//' @param nStarts The number of random starts the analysis should perform. The first start will be performed with the values given by Wstart. The consecutive starts will be Wstart plus a matrix with random uniform values times the current start number (the first start has index zero).
//' @param printLoss A boolean: TRUE will print the lossfunction value each 10th iteration.
//' @return A list containing: \cr
//' \code{W} A matrix containing the component weights \cr
//' \code{P} A matrix containing the loadings \cr
//' \code{loss} A numeric variable ccontaining the minumn loss function value of all the nStarts starts \cr
//' \code{converged} A boolean containing \code{TRUE} if converged \code{FALSE} if not converged. 
//' @export
//' @examples
//'
//' J <- 30
//' X <- matrix(rnorm(100*J), 100, J)
//' ncomp <- 3 
//' constraints <- matrix(1, J, ncomp) #no constraints 
//' 
//' scads(X, ncomp = ncomp, ridge = 10e-8, lasso = rep(1, ncomp), 
//'         constraints = constraints, Wstart = matrix(0, ncomp, J))
//'         
//' # Extended examples:
//' # Example 1: Perform PCA with elistastic net regularization no constraints 
//' #create sample dataset
//' ncomp <- 3 
//' J <- 30
//' comdis <- matrix(1, J, ncomp)
//' comdis <- sparsify(comdis, 0.7) #set 70% of the 1's to zero
//' variances <- makeVariance(varianceOfComps = c(100, 80, 70), J = J, error = 0.05) #create realistic eigenvalues
//' dat <- makeDat(n = 100, comdis = comdis, variances = variances)
//' X <- dat$X
//' 
//' results <- scads(X = X, ncomp = ncomp, ridge = 0.1, lasso = rep(0.1, ncomp),
//'                 constraints = matrix(1, J, ncomp), Wstart = matrix(0, J, ncomp),
//'                 itr = 100000, nStarts = 1, printLoss = TRUE , tol = 10^-8)
//' 
//' head(results$W) #inspect results of the estimation
//' head(dat$P[, 1:ncomp]) #inspect data generating model
//' 
//' 
//' # Example 2: Perform SCA with lasso regularization try out all common dinstinctive structures
//' # create sample data, with common and distinctive structure
//' ncomp <- 3 
//' J <- 30
//' comdis <- matrix(1, J, ncomp)
//' comdis[1:15, 1] <- 0 
//' comdis[15:30, 2] <- 0 
//' 
//' comdis <- sparsify(comdis, 0.2) #set 20 percent of the 1's to zero
//' variances <- makeVariance(varianceOfComps = c(100, 80, 90), J = J, error = 0.05) #create realistic eigenvalues
//' dat <- makeDat(n = 100, comdis = comdis, variances = variances)
//' X <- dat$X
//' 
//' #generate all possible common and distinctive structures
//' allstructures <- allCommonDistinctive(vars = c(15, 15), ncomp = 3, allPermutations = TRUE, filterZeroSegments = TRUE)
//' 
//' #Use cross-validation to look for the data generating structure 
//' index <- rep(NA, length(allstructures))
//' for (i in 1:length(allstructures)) {
//'     print(i)
//'     index[i] <- CVforPCAwithSparseWeights(X = X, nrFolds = 10, FUN = scads, ncomp, ridge = 0, lasso = rep(0.01, ncomp),
//'                 constraints = allstructures[[i]], Wstart = matrix(0, J, ncomp),
//'                 itr = 100000, nStarts = 1, printLoss = FALSE, tol = 10^-5)$MSPE
//' }
//' 
//' #Do the analysis with the "winning" structure
//' results <- scads(X = X, ncomp = ncomp, ridge = 0.1, lasso = rep(0.1, ncomp),
//'                 constraints = allstructures[[which.min(index)]], Wstart = matrix(0, J, ncomp),
//'                 itr = 100000, nStarts = 1, printLoss = TRUE , tol = 10^-5)
//' 
//' head(results$W) #inspect results of the estimation
//' head(dat$P[, 1:ncomp]) #inspect data generating model
//' @references
//' De Schipper, N. C., & Van Deun, K. (2018). Revealing the Joint Mechanisms in Traditional Data Linked With Big 					Data. Zeitschrift Für Psychologie, 226(4), 212–231. doi:10.1027/2151-2604/a000341
// [[Rcpp::export]]
Rcpp::List scads(arma::mat X, int ncomp, double ridge, arma::vec lasso, 
                       arma::mat constraints, int itr, arma::mat Wstart, 
                       double tol = 10e-8, int nStarts = 1,
                       bool printLoss = true){ 
                       
    /* some input checking */
    if ((int) lasso.n_elem != ncomp ) {
        Rcpp::Rcout << "Length lasso argument: " << lasso.n_elem << "\n";

        Rcpp::stop("Lasso argument is not equal to the number of components ncomp. The penalities are specified per component individually. \n");
    }
    if ((int) Wstart.n_cols != ncomp || Wstart.n_rows != X.n_cols) {
        Rcpp::stop("Wstart should have dimensions ncol(X) x ncomp");
    }
    if ((int) constraints.n_cols != ncomp || constraints.n_rows != X.n_cols) {
        Rcpp::stop("constraints should have dimensions ncol(X) x ncomp") ;
    }

    // Initialize objects
    bool converged = false;
    Rcpp::List ret;
    arma::mat W, P, U, V, XtX = X.t() * X;
    arma::vec XtXdiag = XtX.diag(), lossFunctionValue, loss, sumPRESS, D;
    int J = X.n_cols, I = X.n_rows, sign_CP;
    double wold, CP, wols, wnew, minLoss, minLossGlobal = arma::datum::inf; 
 
    // Multi-Start 
    for (int z = 0; z < nStarts; z++) {
        if (arma::accu(Wstart) == 0) {
            /*the first start of the algorithm will be "warm", the consecutive 
              starts will be random starts */
            arma::svd(U, D, V, X);
            W = V.cols(0, ncomp - 1);
            W += arma::randu<arma::mat>(size(W)) * z;
            W.elem(find(constraints == 0)).zeros();
        } else {
            /* If the starting value matrix does not sum to zero,
             * meaning the user wants custom starting values, do:
             */
            W = Wstart;
            W += arma::randu<arma::mat>(size(W)) * z;
            W.elem(find(constraints == 0)).zeros();
        }

        lossFunctionValue = arma::vec(itr + 1);
        lossFunctionValue.fill(arma::datum::inf);

        for (int n = 0; n < itr; n++) {

             // Update P: Procruse rotation
             arma::svd(U, D, V, (XtX * W));
             P = U.cols(0, (ncomp-1)) * V.t();

             // Update W: Coordinate Descent
             for (int q = 0; q < ncomp; q++) {
                
                sumPRESS =  X * (P.col(q) - W.col(q));
                
                for (int j = 0; j < J; j++) {
                    if (constraints(j, q) != 0) {
                        wold = W(j, q);
                        CP  = (1 / double(I)) * (dot((X.col(j)),
                                    sumPRESS.t()) + (wold * XtXdiag(j)));  
                        sign_CP = (CP > 0) ? 1 : ((CP < 0) ? -1 : 0);
                        wols = double(sign_CP) * (std::abs(CP) - lasso(q));
                        wnew = wols / (ridge + XtXdiag(j) * (1/double(I)));

                        if (std::abs(CP) < lasso(q) ) {
                            wnew = 0;
                            sumPRESS += X.col(j) * wold;
                        } else {
                            sumPRESS += (wold - wnew) * X.col(j);
                        }

                        W(j, q) = wnew;
                    }

                }

            }
           
            // Update loss function 
            lossFunctionValue(n) = (1 / double(I * 2)) * accu(pow(X - (X * W * P.t()) , 2)) +
                ((1 / double(2)) * ridge * accu(pow( W, 2))) + dot(lasso, sum(abs(W), 0 )); 

            // Print loss function 
            if (printLoss && n % 10 == 0) {
                Rcpp::Rcout << "Start: " << z+1 << " at itr: " << n+1 << " at loss val: "
                    << lossFunctionValue(n)  << "\n";
            } 
            
            // Evaluate condition to stop inner loop 
            if (n > 0 && lossFunctionValue(n - 1) - lossFunctionValue(n) < tol) {
                if (printLoss) {
                    Rcpp::Rcout << "converged" << "\n";
                }
                converged = true;
                break; 
            }
        
        } 

        loss = lossFunctionValue(find(lossFunctionValue != arma::datum::inf));
        minLoss = loss(loss.size() - 1); 
 
        if (minLoss < minLossGlobal) {
            minLossGlobal = minLoss; 
            ret["W"] = W; 
            ret["P"] = P; 
            ret["loss"] = minLoss; 
            ret["converged"] = converged;
        }   
    }
    return ret;
}




