#' allCommonDistinctive: creates all common and distinctive structures for a component weight matrix 
#'
#' A function that generates all common and distinctive structures given the number of variables per block and the number of components 
#'
#' @param vars A vector indicating the number of variables per data block. Example, c(10, 20, 5), there the first 10 variables belong the first block, the next 20 variables belong to the second block, the next 5 variables belong to the third block
#' @param ncomp The number of components that are of interest
#' @param allPermutations A boolean with default \code{TRUE}. If set to \code{TRUE}, all permutations of common and distinctive will be taken into account (a certain component can occur in all positions). If set to \code{FALSE} all permutations will not be taken into account i.e. a certain component can only occur once. If the intent is to use multiple starts use \code{FALSE} else, use \code{TRUE}
#' @param filterZeroSegments A boolean with default \code{TRUE}. If set to \code{TRUE} components were all blocks are set to zero are not possible reducing the number of combinations.  
#' @return 
#'  A list with all common and distinctive structures
#' @export
#' @examples
#'
#' allCommonDistinctive(vars = c(10, 5, 5), ncomp = 3, allPermutations = TRUE, filterZeroSegments = TRUE)
#' allCommonDistinctive(vars = c(10, 5, 5), ncomp = 3, allPermutations = FALSE, filterZeroSegments = FALSE)
#'  
allCommonDistinctive <- function(vars, ncomp, allPermutations = TRUE, 
                                 filterZeroSegments = TRUE) {

    nblocks <- length(vars)
    W  <- matrix(NA, sum(vars), ncomp)

    # Common and distinctive components can only appear once 
    if (allPermutations == FALSE) {

        cd <- as.matrix(expand.grid(rep(list(0:1), nblocks))[-1, ])
        commonSpecific <- gtools::combinations(n = nrow(cd), r = ncomp,
                                      v = c(1:nrow(cd)), repeats.allowed = TRUE)

        allpossibleWmatrices <- rep(list(NA), nrow(commonSpecific))
        b <- 1
        for (i in 1:nrow(commonSpecific)) {
            for (j in 1:ncol(commonSpecific)) {
                W[, j] <-  rep(cd[commonSpecific[i, j], ], times = vars)
            }
            # A zero row or column, corresponds to blocks or components canceled out
            # This might not be wanted if so, filterZeroSegments == TRUE
            if (filterZeroSegments) {
                reject <- sum(apply(W, 1, sum) %in% 0) + sum(apply(W, 2, sum) %in% 0)
                if (!reject) {
                    allpossibleWmatrices[[b]] <- W
                    b <- b + 1
                } else {
                    allpossibleWmatrices[b] <- NULL
                }
            } else {
                allpossibleWmatrices[[b]] <- W
                b <- b + 1
            }
        }

    # Common and distinctive combinations can appear on more locations 
    } else {

        perm <- gtools::permutations(n = 2, r = ncomp*nblocks,
                             v = c(0, 1), repeats.allowed = TRUE) 
        allpossibleWmatrices <- rep(list(NA), nrow(perm))
        cums  <- c(1, cumsum(vars))
        b <- 1
        for (i in 1:nrow(perm)) {
            for (j in 1:ncomp) {
                for (a in 1:nblocks) {
                    W[cums[a]:cums[a+1], j] <- perm[i, (j-1) * nblocks + a]
                }
            }
            # A zero row or column, corresponds to blocks or components canceled out
            # This might not be wanted if so, filterZeroSegments == TRUE
            if (filterZeroSegments) {
                reject <- sum(apply(W, 1, sum) %in% 0) + sum(apply(W, 2, sum) %in% 0)
                if (!reject) {
                    allpossibleWmatrices[[b]] <- W
                    b <- b + 1
                } else {
                    allpossibleWmatrices[b] <- NULL
                }
            } else {
                allpossibleWmatrices[[b]] <- W
                b <- b + 1
            }
        }
    }

    return(allpossibleWmatrices)
}



