#' Distributed Least square approximation 
#'
#' @param formula a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param K number of workers
#' @param fit_fun user-submitted fitting function
#' @param lasso_fun logical. Should the lars result be included in the final result?
#' @param ind the rule to split the data to workers
#' @param intercept logical. Should an intercept be included in the null model?
#' @param ... optional arguments to fit_fun
#' @return A list of the distributed least squares estimator and covariance matrix
#' @importFrom tibble as.tibble
#' @importFrom pROC roc auc ggroc
#' @importFrom stats contrasts is.empty.model model.matrix model.response
#' @import lars
#' @import ggplot2
#' @export
#' @examples
#' # Not Run
dlsa.fit <- function(formula, data, K, fit_fun, lasso_fun = 1, ind = NULL, intercept = 1, ...)
{
  if (!is.data.frame(data)){
    stop("'data' must be a data.frame, not a matrix or an array")
  }
  
  cal = match.call()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  mlm <- is.matrix(Y)
  N <- if (mlm)
    nrow(Y)
  else length(Y)
  X <- if (!is.empty.model(mt))
    suppressWarnings(model.matrix(mt, mf, contrasts))
  else matrix(0, NROW(Y), 0L)
  X = X[,-1]
  
  if (lasso_fun){
    
    dlsa_result = dlsa(X, Y, K, fit_fun = fit_fun, ind = ind, ...)
    lsa_result = lsa.distribute(dlsa_result$theta_mat, dlsa_result$Sig_inv_list,
                                intercept = intercept, sample_size = N)
    
    z <- list( result = list( theta = dlsa_result$theta, Sig_inv = dlsa_result$Sig_inv,
                              theta_mat = dlsa_result$theta_mat,
                              Sig_inv_list = dlsa_result$Sig_inv_list),
               lsa_result = list( beta.ols = lsa_result$beta.ols, beta.bic = lsa_result$beta.bic,
                                  beta.aic = lsa_result$beta.aic,
                                  aic = lsa_result$aic, bic = lsa_result$bic))
    
  }else{
    
    dlsa_result = dlsa(X, Y, K, fit_fun = fit_fun, ind = ind, ...)

    z <- list( result = list( theta = dlsa_result$theta, Sig_inv = dlsa_result$Sig_inv,
                              theta_mat = dlsa_result$theta_mat,
                              Sig_inv_list = dlsa_result$Sig_inv_list))
    
  }
  
  z
  
}



