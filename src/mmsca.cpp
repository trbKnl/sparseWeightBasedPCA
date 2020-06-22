#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]


/*Calculate the grouplasso portion of the loss function*/
double grouplassoLossCpp(const arma::uvec& groups, const arma::vec& grouplasso, const arma::mat& W, const int& Q) {
    if(arma::accu(grouplasso) == 0) {
        return 0;
    } else {
            double out = 0;
        for (int j = 0; j < Q; j++) {

            unsigned int from = 0;
            int to = -1;

            for (arma::uword i = 0; i < groups.n_elem; i++) { 
                to += groups[i];
                out += std::sqrt(groups[i]) * arma::norm(W.col(j).subvec(from, to)) * grouplasso[j];
                from += groups[i];
            } 
        }
        return out; 
    }
}



/*Calculate the grouplasso portion of the loss function
  note this is a void function*/
void grouplassoPenaltyCpp(const arma::uvec& groups, const arma::mat& W, arma::mat& DG, const int& Q){
    for (int j = 0; j < Q; j++) {

        unsigned int from = 0;
        int to = -1;
 
        for (arma::uword i = 0; i < groups.n_elem; i++) { 

            to += groups[i];
            double filler = std::sqrt(groups(i))/2 * 1/arma::norm(W.col(j).subvec(from, to));

            DG.col(j).subvec(from, to).fill(filler); 
            
            from += groups[i];
        }
    }
    //return DG;
}


//calculate the elitist part of the loss function
double elitistLossCpp(const arma::uvec& groups, const arma::vec& elitistlasso, const arma::mat& W, const int& Q){
    if (arma::accu(elitistlasso) == 0) {
        return 0;
    } else {
        double out = 0;
        for (int j = 0; j < Q; j++) {

            unsigned int from = 0;
            int to = -1;

            for (arma::uword i = 0; i < groups.n_elem; i++) { 

                to += groups[i];
                out += std::pow(arma::norm(W.col(j).subvec(from, to), 1), 2) * elitistlasso[j];
                from += groups[i];

            } 
        }
        return out; 
    }
}



/*Calculate the elitist part of the penalty, note this function is a void function,
 also note that if all elements in W are constraints to zero, DE/arma::abs(W) will be 0/0 == NaN,
 in the coordinate descent step NaN gets multiplied with 0 leading to an error, therefore if DE contains exact zeroes
 which I assume can only be the case if the user supplied them, I will change the values into ones, to prevent an error
 in the coordinate descent step */
void elitistPenaltyCpp(const arma::uvec& groups, const arma::mat& W, arma::mat& DE, const int& Q){
    for (int j = 0; j < Q; j++) {

        unsigned int from = 0;
        int to = -1;
 
        for(arma::uword i = 0; i < groups.n_elem; i++){ 

            to += groups[i];
            double filler = arma::norm(W.col(j).subvec(from, to), 1);

            DE.col(j).subvec(from, to).fill(filler); 
            
            from += groups[i];
        }
    }

    DE.elem(find(DE == 0)).ones();
    DE =  DE / arma::abs(W);
}

/* Calculates the loss function */
long double lossFunctionCpp(const arma::mat& X, const arma::mat& W, const arma::mat& P, const arma::uvec& groups, const arma::vec& ridge, const arma::vec& lasso, const arma::vec& grouplasso, const arma::vec& elitistlasso, const double& Q){

        double out = arma::accu(arma::pow((X - X * W * P.t()), 2))
        +  elitistLossCpp(groups, elitistlasso, W, Q)
        +  grouplassoLossCpp(groups, grouplasso, W, Q)
        +  arma::accu(arma::pow(W, 2)*arma::diagmat(ridge))
        +  arma::accu(arma::abs(W)*arma::diagmat(lasso));

        return out;
} 


/*Updates the elements of W per column, could theoretically be done in parallel 
 * This is a void function */
void coordinateDescentStepCpp(const arma::mat& XtXP, const arma::mat& XtX, arma::mat& W, const arma::mat& constraints, const int& Q, const int& J, const arma::mat& D){
    double num, denom;
    for (int q = 0; q < Q; q++) {
        for (int j = 0; j < J; j++) {
            if (constraints(j, q) != 0) {
                if (j == 0) {
                    num = XtXP(j, q) - arma::as_scalar(XtX.row(j).subvec(1, J-1) * W.col(q).subvec(1, J-1));
                } else if (j == J-1) {
                    num = arma::as_scalar(XtXP(j, q) - XtX.row(j).subvec(0, J-2) * W.col(q).subvec(0, J-2));
                } else {
                    num = arma::as_scalar(XtXP(j, q) - XtX.row(j).subvec(0, j-1) * W.col(q).subvec(0, j-1)
                        - XtX.row(j).subvec(j+1, J-1) * W.col(q).subvec(j+1, J-1)); 
                }
                denom = XtX(j, j) + D(j, q);
                W(j, q) = num / denom;
            }
        }
    }
}


/* Find the mimimum of the majorizing function given P, using coordinate descent */
void coordinateDescentCpp(const arma::mat& X, const arma::mat& XtX, const arma::mat& XtXP, const int& itrCoor, arma::mat& W, const arma::mat& constraints, const arma::mat& P, const arma::mat& D, const int& Q, const int& J, const bool& printLoss){

    double tol = std::pow(10.0, -8);
    double a, b;
    arma::vec loss = arma::vec(itrCoor + 1); 
    loss(0) = arma::datum::inf;
    for (int i = 0; i < itrCoor; i++) {
        a = 0;
        for (int q = 0; q < Q; q++) {
            a += arma::as_scalar(arma::pow(W.col(q), 2).t() * D.col(q));
        }

        b = arma::accu(arma::pow(X - X * W * P.t(), 2));
        loss(i + 1) = a + b;
        if (printLoss && (i+1) % 10 == 0 ) {
            Rcpp::Rcout << "iterations coordinate descent step: "  << i+1 << " at loss value: " <<
               loss(i+1) << "\n";
        }
        if (loss(i) - loss(i + 1) < tol) {
            return;
        }
        coordinateDescentStepCpp(XtXP, XtX, W, constraints, Q, J, D);
    }
}

//store in a std::vector, uvec's containing the row indices of the free coefficients in W
//note: with no argument supplied, find will return the indices containing the non-zero elements
std::vector<arma::uvec> makeConstraints(const arma::mat &constraints, const int &Q) {
    
    std::vector<arma::uvec> out;
    out.reserve(Q);

    for (int q = 0; q < Q; q++) {
        arma::uvec indices = arma::find(constraints.col(q));
        out.push_back(indices);

    }
    return out;

}

//' mmsca: Sparse SCA/PCA with and/or: ridge, lasso, group lasso, elitist lasso regularization
//'
//' This function performs PCA/SCA with and/or: ridge, lasso, group lasso, elitist lasso regularization. This function allows for constraining certain weights to zero.
//'  
//' @param X A data matrix of class \code{matrix}
//' @param ncomp The number of components to estimate (an integer)
//' @param ridge A vector containing a ridge parameter for each column of W separately, to set the same ridge penalty for the component weights W, specify: ridge = \code{rep(value, ncomp)}, value is a non-negative double
//' @param lasso A vector containing a ridge parameter for each column of W separately, to set the same lasso penalty for the component weights W, specify: lasso = \code{rep(value, ncomp)}, value is a non-negative double
//' @param grouplasso A vector containing a grouplasso parameter for each column of W separately, to set the same grouplasso penalty for the component weights W, specify: grouplasso = \code{rep(value, ncomp)}, value is a non-negative double
//' @param elitistlasso A vector containing a elitistlasso parameter for each column of W separately, to set the same elitistlasso penalty for the component weights W, specify: elitistlasso = \code{rep(value, ncomp)}, value is a non-negative double
//' @param groups A vector specifying which columns of \code{X} belong to what block. Example: \code{c(10, 100, 1000)}. The first 10 variables belong to the first block, the 100 variables after that belong to the second block etc.
//' @param constraints A matrix of the same dimensions as the component weights matrix W (\code{ncol(X)} x \code{ncomp}). A zero entry corresponds in constraints corresponds to an element in the same location in W that needs to be constraint to zero. A non-zero entry corresponds to an element in the same location in W that needs to be estimated.
//' @param itr The maximum number of iterations (a positive integer)
//' @param Wstart A matrix of \code{ncomp} columns and \code{nrow(X)} rows with starting values for the component weight matrix W, if \code{Wstart} only contains zeros, a warm start is used: the first \code{ncomp} right singular vectors of \code{X}
//' @param tol The convergence is determined by comparing the loss function value after each iteration, if the difference is smaller than \code{tol}, the analysis is converged. Default value is \code{10e-8}
//' @param nStarts The number of random starts the analysis should perform. The first start will be performed with the values given by \code{Wstart}. The consecutive starts will be \code{Wstart} plus a matrix with random uniform values times the current start number (the first start has index zero).
//' @param printLoss A boolean: \code{TRUE} will print the lossfunction value each 1000 iteration.
//' @param coorDes A boolean with the default \code{FALSE}. If coorDes is \code{FALSE} the estimation of the majorizing function to estimate the component weights W conditional on the loadings P will be found using matrix inverses which can be slow. If set to true the marjozing function will be optimized (or partially optimized) using coordinate descent, in many cases coordinate descent will be faster
//' @param coorDesItr An integer specifying the maximum number of iterations for the coordinate descent algorithm, the default is set to 1. You do not have to run this algorithm until convergence before alternating back to the estimation of the loadings. The tolerance for this algorithm is hardcoded and set to \code{10^-8}. 
//' @return A list containing: \cr
//' \code{W} A matrix containing the component weights \cr
//' \code{P} A matrix containing the loadings \cr
//' \code{loss} A numeric variable containing the minimum loss function value of all the \code{nStarts} starts \cr
//' \code{converged} A boolean containing \code{TRUE} if converged \code{FALSE} if not converged.
//' @export
//' @examples
//'
//' J <- 30
//' X <- matrix(rnorm(100*J), 100, J)
//' ncomp <- 3
//' 
//' #An example of sparse SCA with ridge, lasso, and grouplasso regularization, with 2 groups, no constraints, and a "warm" start
//' mmsca(X = X, 
//'        ncomp = ncomp, 
//'        ridge = rep(10e-8, ncomp),
//'        lasso = rep(1, ncomp),
//'        grouplasso = rep(1, ncomp),
//'        elitistlasso = rep(0, ncomp),
//'        groups = c(J/2, J/2), 
//'        constraints = matrix(1, J, ncomp), 
//'        itr = 1000000, 
//'        Wstart = matrix(0, J, ncomp))
//'
//' # Extended example: Perform SCA with group lasso regularization try out all common dinstinctive structures
//' # create sample data, with common and distinctive structure
//' ncomp <- 3 
//' J <- 30
//' comdis <- matrix(1, J, ncomp)
//' comdis[1:15, 1] <- 0 
//' comdis[15:30, 2] <- 0 
//' 
//' comdis <- sparsify(comdis, 0.1) #set 10 percent of the 1's to zero
//' variances <- makeVariance(varianceOfComps = c(100, 80, 90), J = J, error = 0.05) #create realistic eigenvalues
//' dat <- makeDat(n = 100, comdis = comdis, variances = variances)
//' X <- dat$X
//' 
//' results <- mmsca(X = X, 
//'     ncomp = ncomp, 
//'     ridge = rep(10e-8, ncomp),
//'     lasso = rep(0, ncomp),
//'     grouplasso = rep(5, ncomp),
//'     elitistlasso = rep(0, ncomp),
//'     groups = c(J/2, J/2), 
//'     constraints = matrix(1, J, ncomp), 
//'     itr = 1000000, 
//'     Wstart = matrix(0, J, ncomp))
//' 
//' #inspect results
//' results$W
//' dat$P[, 1:ncomp]
//' 
//' #for model selection functions see mmscaModelSelection() and mmscaHyperCubeSelection()
//' 
// [[Rcpp::export]]
Rcpp::List mmsca(const arma::mat& X, const int& ncomp, const arma::vec& ridge, const arma::vec& lasso,  const arma::vec& grouplasso, 
        const arma::vec& elitistlasso, arma::uvec& groups, const arma::mat& constraints, const int& itr,
        arma::mat Wstart, double tol = 10e-8, int nStarts = 1, bool printLoss = true, bool coorDes = false, int coorDesItr = 1) {

    /* some input checking */
    if ((const int) ridge.n_elem != ncomp || (const int) lasso.n_elem != ncomp
            || (const int) grouplasso.n_elem != ncomp || (const int) elitistlasso.n_elem != ncomp) {
        Rcpp::Rcout << "The number of components ncomp: " << ncomp <<"\n";
        Rcpp::Rcout << "Length ridge argument: " << ridge.n_elem << "\n";
        Rcpp::Rcout << "Length lasso argument: " << lasso.n_elem << "\n";
        Rcpp::Rcout << "Length grouplasso argument: " << grouplasso.n_elem << "\n";
        Rcpp::Rcout << "Length elitistlasso argument: " << elitistlasso.n_elem << "\n";

        Rcpp::stop("Either ridge, lasso, grouplasso or elitistlasso is not equal to the number of components ncomp. The penalities are specified per component individually. \n");
    }
    if (arma::sum(groups) != X.n_cols) {
        for (int i = 0; i < (int) groups.n_elem; i++) {
            Rcpp::Rcout << "Group: " << i+1 << ", is of size: " <<  groups[i] << "\n";
        }
        Rcpp::stop("The sum of the groups is not equal to the number of variables: ncol(X)");
    }
    if ((int) Wstart.n_cols != ncomp || Wstart.n_rows != X.n_cols) {
        Rcpp::stop("Wstart should have dimensions ncol(X) x ncomp");
    }
    
    if ((int) constraints.n_cols != ncomp || constraints.n_rows != X.n_cols) {
        Rcpp::stop("constraints should have dimensions ncol(X) x ncomp") ;
    }

    //Object intialization    
    double minLoss, minLossGlobal = arma::datum::inf;
    bool converged = false;
    int J = X.n_cols;
    arma::mat P = arma::randn(J, ncomp);
    arma::mat DG = arma::randn(J, ncomp);
    arma::mat DE = arma::zeros(J, ncomp);
    arma::mat DR = arma::ones(J, ncomp);  

    arma::mat U, U2, V, V2, DL, Dsup;
    arma::vec D, D2;
    arma::svd(U, D, V, X);
    arma::mat W = V.cols(0, ncomp - 1);
    arma::mat XtX = X.t() * X;
    Rcpp::List ret;

    /* needed for the non-contigues submatrix view method: submat*/ 
    std::vector<arma::uvec> indrow = makeConstraints(constraints, ncomp);
    std::vector<arma::uvec> indcol;
    indcol.reserve(ncomp);
    
    for (arma::uword q = 0; q < arma::uword(ncomp); q++) {
        arma::uvec a = {q};
        indcol.push_back(a);
    } 

    for (int j = 0; j < nStarts; j++) {

        arma::vec loss = arma::vec(itr+1); 
        loss.fill(arma::datum::inf);

        
        /* Rcpp cannot handle matrix intialization therefore
         * the user has to give starting values for the algorithm.
         * if the user supplied a matrix that sums to zero normal
         * initialization follows.         
         */

        if (arma::accu(Wstart) == 0) {
            /*the first start of the algorithm will be "warm", the consecutive 
              starts will be random starts */
            //Rcpp::Rcout << "check" << "\n";
            arma::svd(U, D, V, X);
            W = V.cols(0, ncomp - 1);
            W += arma::randu<arma::mat>(size(W)) * j;
            W.elem(find(constraints == 0)).zeros();
        } else {
            /* If the starting value matrix does not sum to zero,
             * meaning the user wants custom starting values, do:
             */
            W = Wstart;
            W += arma::randu<arma::mat>(size(W)) * j;
            W.elem(find(constraints == 0)).zeros();
        }
        /* start algorithm */
        for (int i = 0; i < itr; i++) {
            Rcpp::checkUserInterrupt();

            loss(i + 1) = lossFunctionCpp(X, W, P, groups, ridge, lasso, grouplasso, elitistlasso, ncomp);

            if (printLoss && i % 10 == 0) {
                Rcpp::Rcout << "Start: " << j+1 << " at itr: " << i+1 << " at loss val: " << loss(i+1) << "\n";
            }

            //procruste rotation least squares P given W
            arma::svd(U2, D2, V2, XtX * W);
            P = U2.cols(0, ncomp-1) * V2.t();
            
            //lasso
            DL = arma::pow(arma::abs(W), -1);
            DL.elem(find(DL > std::pow(10.0, 6))).fill(std::pow(10.0, 6));
            //group lasso
            grouplassoPenaltyCpp(groups, W, DG, ncomp); 
            DG.elem(find(DG > std::pow(10.0, 6))).fill(std::pow(10.0, 6));

            //elitist lasso
            if (arma::accu(elitistlasso) != 0) {
                elitistPenaltyCpp(groups, W, DE, ncomp); 
                DE.elem(find(DE > std::pow(10.0, 6))).fill(std::pow(10.0, 6));
            }

            //scale each column of the diagonal matrices by multiplying with penalties in diagonal matrices 
            Dsup = DL*arma::diagmat(lasso/2) + DR*arma::diagmat(ridge) + DG*arma::diagmat(grouplasso/2)
                + DE*arma::diagmat(elitistlasso);

            arma::mat XtXP = XtX * P;

            /* Given P find solution for W with coordinate descent or
            through inverses*/
            if (coorDes) {
                coordinateDescentCpp(X, XtX, XtXP, coorDesItr, W, constraints, P, Dsup, ncomp, J, printLoss);
            } else {
                for (int q = 0; q < ncomp; q++) {
                    W.submat(indrow[q], indcol[q]) =
                        (arma::diagmat(Dsup.submat(indrow[q], indcol[q])) +
                         XtX.submat(indrow[q], indrow[q])).i() * XtXP.submat(indrow[q], indcol[q]);
                }
            }

            //if converged break
            if (loss(i) - loss(i+1) < tol) {
                if (printLoss) {
                    Rcpp::Rcout << "converged" << "\n";
                }
                converged = true;
                break; 
            }
        }

        loss = loss(find(loss != arma::datum::inf));
        minLoss = loss(loss.size() - 1);       

        if (minLoss < minLossGlobal) {

            minLossGlobal = minLoss;
            //set small elements to zero
            W.elem(find(arma::abs(W) < std::pow(10.0, -4))).zeros(); 
            ret["W"] = W;
            ret["P"] = P; 
            ret["loss"] = minLoss;
            ret["converged"] = converged;
        } 
    }

    return ret;
}


