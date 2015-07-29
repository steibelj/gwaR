
#' A class for GWA from a GBLUP object.  
#' @title Estimation of SNP effects and their variances 
#' @param gblup An object of class gblup
#' @param x A matrix of standardized genotypes
#' @return an object of the class GWAS: a list with two vectors
#' \itemize{ 
#'    \item {\code{ghat}} {estimated SNP effects} 
#'    \item{\code{varg}} {estimated SNP effect variances}
#'}
#'  @seealso \code{\link{summary.gwas}} \code{\link{plot.gwas}}
#'  @export


gwas <- function(gblup, x) UseMethod("gwas")


#'  @rdname gwas
#'  @export

gwas.default <- function(gblup, x) {
    
    # Check points Check that the class of gwas is gblup
    if (class(gblup) != "gblup") 
        stop("input must be a gblup object")
    
    Vinve <- gblup[["Vinv"]] %*% gblup[["ehat"]]
    
    idx <- Reduce(intersect, list(colnames(x), rownames(gblup[["Vinv"]])))
    
    # calculate the g hat
    ghat <- gblup[["sigma"]][["G"]] * crossprod(t(x[, idx]), Vinve[idx, ])
    
    # Calculate Vinv Q V t(Q) Vinv
    VinQVQtVin <- gblup[["Vinv"]] %*% gblup[["Q"]] %*% gblup[["V"]] %*% t(gblup[["Q"]]) %*% gblup[["Vinv"]]
    colnames(VinQVQtVin) <- rownames(VinQVQtVin)
    
    # apply the function cros prod
    tZViQVQViZ <- diagprod(x[, idx], VinQVQtVin[idx, idx])
    
    
    # Calculate the Variance of the G hat
    var_ghat <- tZViQVQViZ * (gblup[["sigma"]][["G"]])^2
    
    g_res <- data.frame(ghat = ghat, var_ghat = var_ghat)
    class(g_res) <- c("gwas", "data.frame")
    return(g_res)
} 
