
#' A class for GWA from a GBLUP object.  
#' @title Estimation of SNP effects and their variances 
#' @param gblup An object of class gblup
#' @param x A matrix of standardized genotypes
#' @return an object of the class GWAS: a list with two vectors
#' \itemize{ 
#'    \item {\code{ghat}} {estimated SNP effects} 
#'    \item{\code{varg}} {estimated SNP effect variances}
#'}
#'  @seealso \code{\link{summary.GWAS}} \code{\link{plot.GWAS}}
#'  @export


GWAS <- function(gblup, x) UseMethod("GWAS")


#'  @rdname GWAS
#'  @export

GWAS.default <- function(gblup, x) {
    
    # Check points Check that the class of gwas is gblup
    if (class(gblup) != "gblup") 
        stop("input must be a gblup object")
    
    # Filter x so it has the same observations that were previousley filter in GLBUP
    idxa <- rownames(gblup[["Vinv"]]) %in% colnames(x)
    gblup[["Vinv"]] <- gblup[["Vinv"]][idxa, idxa]
    idx <- match(rownames(gblup[["Vinv"]]), colnames(x), )
    x <- x[, idx]
    
    # filter all the other parts of the gblup
    gblup[["ehat"]] <- gblup[["ehat"]][idxa]
    gblup[["Q"]] <- gblup[["Q"]][idxa, idxa]
    gblup[["V"]] <- gblup[["V"]][idxa, idxa]
    
    # calculate the g hat
    ghat <- gblup[["sigma"]][["G"]] * (x %*% (gblup[["Vinv"]] %*% gblup[["ehat"]]))
    
    # Calculate Vinv Q V t(Q) Vinv
    VinvQ <- gblup[["Vinv"]] %*% gblup[["Q"]]
    VinQVQtVin <- VinvQ %*% gblup[["V"]] %*% t(VinvQ)
    
    # apply the function cros prod
    tZViQVQViZ <- diagprod(x, VinQVQtVin)
    
    # Calculate the Variance of the G hat
    var_ghat <- tZViQVQViZ * (gblup[["sigma"]][["G"]])^2
    
    g_res <- data.frame(ghat = ghat, var_ghat = var_ghat)
    class(g_res) <- c("GWAS", "data.frame")
    return(g_res)
} 
