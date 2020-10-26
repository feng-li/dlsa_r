#' The example of logistic fit function
#'
#' @param X design matrix of dimension n * p.
#' @param Y vector of observations of length n, or a matrix with n rows.
#' @param init.beta initial beta values for the logistic regression.
#' @param iter.max max iterate times for the logistic regression.
#' @return A list of the estimated beta and hess matrix.
logistic_fit <- function(X, Y, init.beta = NULL, iter.max = NULL)
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
  return(list(theta = beta, Sig_inv = hess))
}


#' The example of logistic predict function
#'
#' @param X design matrix of dimension n * p.
#' @param theta the estimated value of beta.
#' @param alpha the critical point to obtain Yhat from Fitted result.
#' @return A list of the product of x and beta, probability and estimated value of observations.
logistic_pred <- function(X, theta, alpha)
{
  xtheta = X %*% theta
  p = exp(xtheta)/(1 + exp(xtheta))
  yhat = array(0,length(p))
  for (i in 1: length(p)){
    if (p[i] >= alpha){
      yhat[i] = 1
    }else{
      yhat[i] = 0
    }
  }
  return(list(xtheta = xtheta, p = p, yhat = yhat))
}


#' The example of plot function
#'
#' @param x the coordinates of points in the plot.
#' @param y the y coordinates of points in the plot.
#' @return the plot result.
plot_fun <- function(x,y){
  plot = ggplot()+
    geom_point(mapping = aes(x = x, y = y))
  return(plot)
}


#' the simulator function for covariates X
#'
#' @param N the row dimension of simulated X.
#' @param p the column dimension of simulated X.
#' @param K number of workers.
#' @param iid logical. whether the simulated X obey independently identically distribution or not.
#' @return the simulated X values.
#' @import MASS
#' @import survival
logistic_simuX<-function(N, p, K, iid = T)
{
  if (iid)
    X = matrix(rnorm(N*p), ncol = p)
  else{
    nks = rep(floor(N/K), K-1)
    nks[K] = N - sum(nks)
    mus = seq(-1, 1, length.out = K)
    mu_mat = rep(1, p)%*%t(mus)
    sigs = seq(0.3, 0.4, length.out = K)
    X_list = lapply(1:K, function(i){
      XSigma = sigs[i]^abs(outer(1:p,1:p,"-"))
      Xi = mvrnorm(n = nks[i], mu = mu_mat[,i], Sigma = XSigma)
    })
    X = do.call(rbind, X_list)
  }
  return(X)
}


#' the simulator function to simulate the Y according to regression model
#'
#' @param X design matrix of dimension n * p.
#' @param beta the beta value of regression.
#' @param reg_type the type of regression.
#' @return the simulated Y values.
logistic_simuY<-function(X, beta, reg_type = "logistic")
{
  if (reg_type == "logistic")
  {
    N = nrow(X)
    prob = exp(X%*%beta)/(1+exp(X%*%beta))
    Y = rbinom(N,size = 1, prob)
  }
  if (reg_type == "cox")
  {
    n = nrow(X)
    xt <- X%*%beta
    yt = rexp(n, exp(xt))
    ut <- rexp(n, 1/(runif(n)*2+1)*exp(-xt))
    Y <- Surv(pmin(yt,ut),yt<=ut)
  }
  if (reg_type == "linear")
  {
    n = nrow(X)
    xbeta <- X%*%beta
    Y <- xbeta + rnorm(n)
  }
  if (reg_type == "poisson")
  {
    n = nrow(X)
    xbeta <- X%*%beta
    lamb = exp(xbeta)
    Y = rpois(n, lamb)
  }
  if (reg_type == "ordered_probit")
  {
    n = nrow(X)
    xbeta <- as.vector(X%*%beta)
    cutoff = c(-1, 0, 0.8)
    cxbeta = t(outer(cutoff, xbeta, "-"))
    cumprob = pnorm(cxbeta)
    cumprob1 = cbind(cumprob,1)
    probs = cumprob1 - cbind(0, cumprob)

    Y = t(apply(probs, 1, function(x) rmultinom(1, size = 1, prob = x)))
  }

  return(Y)
}
