#' A class for GBLUP prediction to implement GWA.  
#' @title prediction of genomic breeding values for use in GWA 
#' @param rsp a string with the name of the response variable if a vector, first element is used. 
#' @param data either a gpData object from synbreed or a dataframe containing records for response variable and classification factors used as fixed or random effects  
#' @param design a list or a vector with formulas for fixed and random effects
#' @param vdata a list with proportional covariance matrices of additional random effects (not included in data)
#' @param G a proportional covariance for genetics effects being predicted. 
#' @param vdata a list with proportional covariance matrices of additional random effects (not included in data)
#' @param wt a vector of weigths (proportional to diagonals of residual covariance matrix), or a matrix of weights with anmimals in rows and response variables in columns
#' @param Res a logical value to indicate if a residual varaiance should be fit. It is T by default and it can be set to F if a proportional residual covariance matrix is provided.
#' @param ... additional parameters for regress function

#' @return an object of the class gblup: a list the following components
#' \itemize{ 
#'    \item {\code{name}} {trait name} 
#'    \item{\code{llik}} {log-likelihood value at parameter estimates}
#'    \item{\code{cycle}} {number of iterations to convergence}
#'    \item{\code{Q}} {a matrix: see regress}
#'    \item{\code{V}} {phenotypic covariance matrix}
#'    \item{\code{Vinv}} {inverse of V}
#'    \item{\code{sigma}} {vector of estimated variance components}
#'    \item{\code{coefm}} {matrix of estimated fixed effects coefficients }
#'    \item{\code{model}} {the fit model}
#'    \item{\code{ehat}} {vector of residuals}
#'    \item{\code{pos}} {see regress}
#'}
#'  @export

gblup <- function(rsp, data, design, G, vdata = NULL, wt = NULL, Res = T, ...) UseMethod("gblup")

#'  @rdname gblup
#'  @export
gblup.default <- function(rsp, data, design, G, vdata = NULL, wt = NULL, Res = T, ...) {
    
    if (!is.matrix(G)) 
        stop("G should be a (relationship) matrix")
    
    if (!is.null(vdata)) 
        if (!is.list(vdata)) 
            stop("vdata should be a list of covariance matrices ")
    vdata$G <- G
    
    if (!is.character(rsp)) 
        stop("resp should be a string with trait name")
    rsp <- rsp[1]
    # design could be a formula (cohersible) or a list of length 1,2
    
    if (is.list(design) & (length(design) > 2)) 
        stop("list design should have only two components")
    
    # if ts is a list of length 2: it gets only the LHS
    if (is.list(design) & (length(design) == 2)) {
        fm2 <- nlme:::getCovariateFormula(as.formula(design[[2]]))
        fm1 <- nlme:::getCovariateFormula(as.formula(design[[1]]))
    } else {
        # if is a list of length 1 or formula, gets the LHS of fixed effects and defines null random effects
        fm2 <- as.formula("~1")
        if (is.list(design)) {
            fm1 <- nlme:::getCovariateFormula(as.formula(design[[1]]))
        } else {
            fm1 <- nlme:::getCovariateFormula(as.formula(design))
        }
    }
    # collapses formulas to use current environment
    fm1 <- as.formula(paste(as.character(fm1), collapse = " "))
    fm2 <- as.formula(paste(as.character(fm2), collapse = " "))
    
    # adds response to fixed effects formula AND genetic effects to random effects formula
    fm1 <- update(fm1, y ~ .)
    fm2 <- update(fm2, ~G + .)
    
    # names of fixed effects
    nms <- all.vars(fm1)[-1]
    # names of random effects, except for genetic effects
    nmr <- all.vars(fm2)
    
    # check if effects and response is present
    if (class(data) == "gpData") {
        fnotgp <- c(nms[!nms %in% colnames(data$covar)], rsp[!rsp %in% colnames(data$pheno)])
        rnotgp <- nmr[!nmr %in% colnames(data$covar)]
    } else {
        fnotgp <- c(nms, rsp)[!c(nms, rsp) %in% colnames(data)]
        rnotgp <- c(nmr)[!c(nmr) %in% colnames(data)]
    }
    
    rnotvd <- nmr[!nmr %in% names(vdata)]
    
    if (rsp %in% fnotgp) 
        stop(paste("response variable", rsp, "not present in data"))
    
    if (length(intersect(rnotgp, rnotvd)) > 0) 
        stop(paste("these random effects are not present neither in data not in vdata:", intersect(rnotgp, rnotvd)))
    
    if (length(fnotgp) > 0) 
        stop(paste("these fixed effects are not present in data", fnotgp))
    
    
    
    if (class(data) == "gpData") {
        y <- na.omit(data$pheno[, rsp, 1])
        ef <- na.omit(data$covar[, c("id", nms[nms %in% colnames(data$covar)], rnotvd), drop = F])
        rownames(ef) <- ef[, "id"]
    } else {
        y <- data[, rsp]
        names(y) <- rownames(data)
        y <- na.omit(y)
        excols <- c(nms, rnotvd)
        ef <- na.omit(data[, excols, drop = F])
    }
  
    Ind <- Res
    
    if (!is.null(wt)) {
        if (length(dim(wt)) >= 2) {
            if (!rsp %in% colnames(wt)) 
                stop(paste("name of response variable", resp, "should also be present in wt matrix"))
            nm <- rownames(wt)
            wt <- wt[, rsp]
            names(wt) <- nm
            
        }
        fm2 <- update(fm2, ~. + wt)
        Ind <- FALSE
    }
    
    
    addnm <- unlist(lapply(vdata, function(x) ifelse(is.null(dim(x)), return(names(x)), return(rownames(x)))))
    addnm <- names(table(addnm))[table(addnm) == length(rnotgp)]
    
    
    # prepare for creating index based on variable number of inputs
    
    ny <- names(y)
    ifelse(is.null(ef), ne <- ny, ne <- rownames(ef))
    ifelse((is.null(addnm)), nar <- ny, nar <- addnm)
    
    idx <- Reduce(intersect, list(ny, ne, nar))
    
    if (is.null(idx)) 
        stop(paste("no observations in common between response, fixed effects and random effects, please check names"))
    
    y <- y[idx]
    # G <- G[idx, idx]
    ifelse(is.null(ef), NA, ef <- ef[idx, , drop = FALSE])
    ifelse(is.null(wt), NA, wt <- diag(wt[idx]))
    
    
    
    for (i in 1:length(vdata)) {
        vdata[[i]] <- vdata[[i]][idx, idx]
    }
    
    
    x <- regress(fm1, fm2, data = c(ef, vdata), identity = Ind, ...)  #possible conflict if the user specifies
    # identity in the ...
    
    coefm <- rbind(cbind(x$beta, sqrt(diag(x$beta.cov))), cbind(x$sigma, sqrt(diag(x$sigma.cov))))
    colnames(coefm) <- c("Estimate", "StdError")
    
    h2 <- x$sigma[["G"]]/sum(x$sigma)
    
    # add an option to compact the output or add the compact version to the summary and withi n the summary add the likihood
    # ratio test save QVQ'
    
    rst <- list(name = rsp, llik = x$llik, cycle = x$cycle, Q = x$Q, V = x$Sigma, Vinv = x$W, sigma = x$sigma, coefm = coefm, 
        model = x$model, ehat = y - x$fitted, pos = x$pos)
    class(rst) <- "gblup"
    return(rst)
} 
