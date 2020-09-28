#' Distributed Least square approximation of logistic regression
#'
#' @param formula a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param K number of workers
#' @param lasso_func logical. Should the lars result be included in the final result?
#' @param alpha the critical point to obtain Yhat from Fitted result
#' @param init.beta initial beta values for the logistic regression
#' @param iter.max max iterate times for the logistic regression
#' @param ind the rule to split the data to workers
#' @param intercept logical. Should an intercept be included in the null model?
#' @return A list of the distributed least squares estimator and covariance matrix
#' @importFrom tibble as.tibble
#' @importFrom pROC roc auc ggroc
#' @importFrom stats contrasts is.empty.model model.matrix model.response
#' @import lars
#' @import ggplot2
#' @import MASS
#' @import survival
#' @export
#' @examples
#' # Not Run
dlsa.logistic <- function(formula, data, K, lasso_func = 1, alpha = 0.5,
                         init.beta = NULL, iter.max = NULL, ind = NULL, intercept = 1)
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


  if (lasso_func){

    dlsa_result = dlsa(X, Y, K, init.beta = init.beta, iter.max = iter.max, ind = ind)
    las_result = lsa.distribute(dlsa_result$theta_mat, dlsa_result$Sig_inv_list,
                                intercept = intercept, sample_size = N)

    dlsa_predict = logistic.predict(X, dlsa_result$theta, K, ind = ind, alpha = alpha)
    dlsa_predict_plot = ggplot() +
                          geom_point(mapping = aes(x = dlsa_predict$Xtheta, y = dlsa_predict$p))

    dlsa_roc = roc(dlsa_predict$Yhat, Y)
    dlsa_roc_plot = ggroc(dlsa_roc,size = 1)+theme_minimal() + ggtitle("DLSA ROC curve") +
                geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
    dlsa_auc = auc(dlsa_roc)

    z <- list( dlsa_result = list( dlsa_theta = dlsa_result$theta, dlsa_Sig_inv = dlsa_result$Sig_inv,
                                   dlsa_theta_mat = dlsa_result$theta_mat,
                                   dlsa_Sig_inv_list = dlsa_result$Sig_inv_list),
               las_result = list( beta.ols = las_result$beta.ols, beta.bic = las_result$beta.bic,
                                  beta.aic = las_result$beta.aic,
                                  aic = las_result$aic, bic = las_result$bic),
               dlsa_predict = list( dlsa_predict = dlsa_predict$Yhat,dlsa_predict_plot = dlsa_predict_plot,
                                    dlsa_roc_plot = dlsa_roc_plot, dlsa_auc = dlsa_auc ))

  }else{
    dlsa_result = dlsa(X, Y, K, init.beta = init.beta, iter.max = iter.max, ind = ind)

    dlsa_predict = logistic.predict(X, dlsa_result$theta, K, ind = ind, alpha = alpha)
    dlsa_predict_plot = ggplot() +
                          geom_point(mapping = aes(x = dlsa_predict$Xtheta, y = dlsa_predict$p))

    dlsa_roc = roc(dlsa_predict$Yhat, Y)
    dlsa_roc_plot = ggroc(dlsa_roc,size = 1)+theme_minimal() + ggtitle("DLSA ROC curve") +
               geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
    dlsa_auc = auc(dlsa_roc)

    z <- list( dlsa_result = list( dlsa_theta = dlsa_result$theta, dlsa_Sig_inv = dlsa_result$Sig_inv,
                                   dlsa_theta_mat = dlsa_result$theta_mat,
                                   dlsa_Sig_inv_list = dlsa_result$Sig_inv_list),
               dlsa_predict = list( dlsa_predict = dlsa_predict$Yhat,dlsa_predict_plot = dlsa_predict_plot,
                                    dlsa_roc_plot = dlsa_roc_plot, dlsa_auc = dlsa_auc ) )
  }

  z

}



###################################
##prediect logistic function with theta
logistic.predict <- function(X, theta, K, ind = NULL, alpha = 0.5){

  # split the data
  N = nrow(X)
  if (is.null(ind))
  {
    ind = rep(1:K, each = floor(N/K))
  }

  Xtheta = list()
  p = list()
  Yhat = list()

  for (k in 1:K){

    Xtheta[[k]] = X[ind==k,] %*% theta
    p[[k]] = exp(Xtheta[[k]])/(1 + exp(Xtheta[[k]]))
    Yhat[[k]] = array(0,length(p[[k]]))
    for (i in 1: length(p[[k]])){
      if (p[[k]][i] >= alpha){
        Yhat[[k]][i] = 1
      }else{
        Yhat[[k]][i] = 0
      }
    }
  }
  Xtheta = as.vector(unlist(Xtheta))
  p = as.vector(unlist(p))
  Yhat = as.vector(unlist(Yhat))
  return(list(Xtheta = Xtheta, p = p, Yhat = Yhat))
}

## logistic regression
logistic<-function(X, Y, init.beta = NULL, iter.max = NULL)
{
  N = nrow(X)
  p = ncol(X)
  if (is.null(init.beta))
  {
    beta = rep(0, p)
  }else{
    beta = init.beta
  }
  del = 1
  iter = 0
  if (is.null(iter.max)){
    iter.max = 80
  }
  while(max(abs(del)) > 10^{-4} & iter < iter.max)
  {

    Xbeta = X%*%beta
    prob01 = 1/as.vector(1+exp(Xbeta))*cbind(1, exp(Xbeta))
    obj = sum(Y*log(prob01[,2]))+sum((1-Y)*prob01[,1])

    # obtain gradient and Hessian matrix
    prob = as.vector(exp(Xbeta)/(1+exp(Xbeta)))
    grad = colSums((Y - prob)*X)
    hess = t(X)%*%(prob*(1-prob)*X)
    del = solve(hess)%*%grad
    if (max(abs(del)) > 5){
      del = del * 0.01
    }
    beta = beta + del
    iter = iter + 1
  }
  return(list(theta = beta, Sig_inv = hess, iter = iter))
}
