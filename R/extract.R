#' extract the yhat part from result list
#' @param x the result list
extrac_yhat <- function(x){
  x$yhat
}
#' extract the Sig_inv part from result list
#' @param x the result list
extrac_Sig_inv <-function(x){
  x$Sig_inv
}
#' extract the theta part from result list
#' @param x the result list
extrac_theta <-function(x){
  x$theta
}
#' calculate the product of Sig_inv and theta from result list
#' @param x the result list
extrac_Sig_inv_theta <-function(x){
  x$Sig_inv%*%x$theta
}


