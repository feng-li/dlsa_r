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
#' @import foreach
#' @import doParallel
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
    ind = rep(1:(K-1), each = floor(N/K))
    ind = c(ind, rep(K, N-floor(N/K)*(K-1)))
  }
  
  # Parallel Computing
  cl <- makeCluster(K)
  registerDoParallel(cl)
  pred_res <- foreach(i=1:K) %dopar% pred_fun(X = X[ind==i,], theta, ...)
  stopImplicitCluster()
  stopCluster(cl)
  
  # get the result
  Yhat = as.vector(unlist(lapply(res$result,extrac_yhat)))
  
  return(list(Yhat = Yhat, result = pred_res))
}
