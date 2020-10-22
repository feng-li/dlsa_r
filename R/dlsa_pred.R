#' Use distributed methods to predict dependent variables based on estimated parameters
#'
#' @param X a symbolic description of the model to be fitted.
#' @param theta a vector containing the estimated parameters of the model.
#' @param K number of workers
#' @param pred_fun user-submitted predicting function
#' @param ind the rule to split the data to workers
#' @param ... optional arguments to pred_fun
#' @return A list of the estimated value of dependent variable and other output of the pred_fun
#' @importFrom tibble as.tibble
#' @export
#' @examples
#' # Not Run
dlsa.pred <- function(X, theta, K, pred_fun, ind = NULL, ...)
{
  
  # split the data
  X = as.matrix(X)
  N = nrow(X)
  if (is.null(ind))
  {
    ind = rep(1:K, each = floor(N/K))
  }
  
  Yhat = list()
  result = list()
  
  for (k in 1:K){
    
    pred_result = pred_fun(X[ind==k,], theta, ...)
    Yhat[[k]] = pred_result$yhat
    result[[k]] = pred_result
    
  }
  
  Yhat = as.vector(unlist(Yhat))
  return(list(Yhat = Yhat, result = result))
}
