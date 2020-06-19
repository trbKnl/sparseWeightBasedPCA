#include <RcppArmadillo.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]


void updatewCpp(arma::mat& W, const arma::mat& P, const int& Q, const arma::mat& XtX, 
    const arma::vec& nzeros, const double& a) {

    arma::mat Sxy = XtX * P;

    for (int i = 0; i < Q; i++) {
        W.col(i) = W.col(i) - (1/a) * (XtX * W.col(i) - Sxy.col(i)); 
        //add cardinality constraints
        arma::uvec indices = arma::sort_index(arma::abs(W.col(i)));
        for (int j = 0; j < nzeros[i]; j++) {
            W(indices[j], i) = 0;
        }
    }

}

//' ccpca: Sparse pca with cardinality constraints on the component weights
//'
//' This function performs PCA with cardinality constraints on the component weights.
//' 
//' @param X A data matrix of class 'matrix'
//' @param ncomp The number of components to estimate (an integer)
//' @param nzeros A vector of length ncomp containing the number of desired zeros in the columns of the component weight matrix \code{W}
//' @param itr The maximum number of iterations (an integer)
//' @param Wstart A matrix of \code{ncomp} columns and \code{nrow(X)} rows with starting values for the component weight matrix \code{W}, if \code{Wstart} only contains zeros, a warm start is used: the first \code{ncomp} right singular vectors of \code{X}
//' @param nStarts The number of random starts the analysis should perform. The first start will be performed with the values given by \code{Wstart}. The consecutive starts will be \code{Wstart} plus a matrix with random uniform values times the current start number (the first start has index zero). The default value is 1.
//' @param tol The convergence is determined by comparing the loss function value after each iteration, if the difference is smaller than \code{tol} the analysis is converged. The default value is \code{10e-8}
//' @param printLoss A boolean: \code{TRUE} will print the loss function value each 10th iteration.
//' @return A list containing: \cr
//' \code{W} A matrix containing the component weights \cr
//' \code{P} A matrix containing the loadings \cr
//' \code{loss} A numeric variable containing the minimum loss function value of all the \code{nStarts} starts \cr
//' \code{converged} A boolean containing \code{TRUE} if converged \code{FALSE} if not converged.
//' @export
//' @examples
//'
//' I <- 100
//' J <- 50 
//' ncomp <- 3
//' X <- matrix(rnorm(I*J), I, J)
//' 
//' ccpca(X = X, ncomp = ncomp,  nzeros = c(10, 20, 30), itr = 100000, 
//'      Wstart = matrix(0, J, ncomp), nStarts = 1, tol = 10^-8, printLoss = TRUE)
//' 
//' # Extended example: Perform CCPCA, with oracle information
//' # create sample data
//' ncomp <- 3 
//' J <- 30
//' comdis <- matrix(1, J, ncomp)
//' 
//' comdis <- sparsify(comdis, 0.5) #set 10 percent of the 1's to zero
//' variances <- makeVariance(varianceOfComps = c(100, 80, 90), J = J, error = 0.05) #create realistic eigenvalues
//' dat <- makeDat(n = 100, comdis = comdis, variances = variances)
//' X <- dat$X
//' 
//' #check how many zero's are in the data generating model
//' nzeros <- apply(dat$P[, 1:ncomp], 2, function(x) {return(sum(x == 0))} )
//' nzeros
//' 
//' #run the analysis with oracle information of the exact number of zero's in the component weights 
//' results <- ccpca(X = X, ncomp = ncomp,  nzeros = nzeros, itr = 10000000, 
//'       Wstart = matrix(0, J, ncomp), nStarts = 1, tol = 10^-8, printLoss = TRUE)
//' 
//' #inspect the results
//' head(results$W) 
//' head(dat$P[, 1:ncomp])
//' 
// [[Rcpp::export]]
Rcpp::List ccpca(const arma::mat& X, const int& ncomp, const arma::vec& nzeros, 
        const int &itr, arma::mat Wstart, int nStarts = 1, double tol = 10e-8, bool printLoss = true) {


    const arma::mat XtX = X.t() * X;
    arma::mat U, V, U2, V2, P, W;
    arma::vec D, D2;
    Rcpp::List ret;

    double prevLoss, curLoss, minLossGlobal = arma::datum::inf;
    bool converged;

    arma::svd(U, D, V, X);
    double a = std::pow(arma::max(D), 2);

    for (int k = 0; k < nStarts; k++) {
        prevLoss = arma::datum::inf;
        curLoss = std::pow(10, 100);
        converged = false;

        /* First run will be warm, or with the custom starting values */
        if(arma::accu(Wstart) == 0){
            arma::svd(U, D, V, X);
            W = V.cols(0, ncomp - 1);
            W += arma::randu<arma::mat>(size(W)) * k;
        } else {
            /* If the starting value matrix does not sum to zero,
             * meaning the user wants custom starting values, do:
             */
            W = Wstart;
            W += arma::randu<arma::mat>(size(W)) * k;
        }

        for (int i = 0; i < itr; i++) {
            Rcpp::checkUserInterrupt();
            //check condition if true break else continue
            if (prevLoss - curLoss < tol) {
                if (printLoss) {
                    Rcpp::Rcout << "converged" << "\n";
                }
                converged = true;
                break; 
            }
            //print loss 
            if (printLoss && i % 10 == 0) {
                Rcpp::Rcout << "Start: " << k+1 << " at itr: " << i << " at loss val: " << curLoss << "\n";
            }

            //procruste rotation least squares P given W
            arma::svd(U2, D2, V2, XtX * W);
            P = U2.cols(0, ncomp-1) * V2.t(); 
            
            //updateW
            updatewCpp(W, P, ncomp, XtX, nzeros, a);
            
            prevLoss = curLoss;
            curLoss =  arma::accu(arma::pow((X - X * W * P.t()), 2));
        }

        if(curLoss < minLossGlobal){
            minLossGlobal = curLoss;

            //set small elements to zero
            ret["W"] = W;
            ret["P"] = P; 
            ret["loss"] = minLossGlobal;
            ret["converged"] = converged;
        } 
    }

    return ret;
}


