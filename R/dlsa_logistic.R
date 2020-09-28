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



###################################
## dlsa for distributed logistic regression
## K is number of workers

dlsa<-function(X, Y, K, init.beta = NULL, iter.max = NULL, ind = NULL)
{
  # data dimension
  p = ncol(X)
  N = nrow(X)

  # split the data
  if (is.null(ind))
  {
    ind = rep(1:K, each = floor(N/K)) #sample(1:R, N, replace = T)
  }

  # conduct logistic regression in parallel
  Sig_inv_list = list()
  Sig_inv_theta_list = list()

  theta_mat = matrix(0, nrow = p, ncol = K)
  for (k in 1:K)
  {

    lr_k = logistic(X = X[ind==k,], Y = Y[ind==k], init.beta, iter.max)
    Sig_inv_list[[k]] = lr_k$Sig_inv
    Sig_inv_theta_list[[k]] = lr_k$Sig_inv%*%lr_k$theta
    theta_mat[,k] = lr_k$theta

  }
  Sig_inv_sum = Reduce("+", Sig_inv_list)
  Sig_inv_theta_sum = Reduce("+", Sig_inv_theta_list)
  theta = solve(Sig_inv_sum)%*%Sig_inv_theta_sum
  return(list(theta = theta, theta_onehot = rowMeans(theta_mat), Sig_inv = Sig_inv_sum,
              theta_mat = theta_mat, Sig_inv_list = Sig_inv_list))
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



###################################
## Least square approximation. This version May, 2019
## Reference Wang, H. and Leng, C. (2006)
## Rewritten by Xuening Zhu
## Comments and suggestions are welcome
##
## Input
## distributed Sigma Inverse and estiamtors
##
## Output
## beta.ols: the MLE estimate
## beta.bic: the LSA-BIC estimate
## beta.aic: the LSA-AIC estimate

lsa.distribute<-function(theta_mat, Sig_inv_list, intercept = 1, sample_size)
{
  K = ncol(theta_mat)
  Sig_inv_theta_list = list()
  for (k in 1:K)
  {
    Sig_inv_theta_list[[k]] = Sig_inv_list[[k]]%*%theta_mat[,k]
  }
  Sig_inv = Reduce("+", Sig_inv_list)
  beta.ols = solve(Sig_inv)%*%Reduce("+", Sig_inv_theta_list)

  l.fit = lars.lsa(Sig_inv, beta.ols, intercept = intercept, n = sample_size)
  t1 <- sort(l.fit$BIC, ind=T)
  t2 <- sort(l.fit$AIC, ind=T)
  beta <- l.fit$beta
  if(intercept) {
    beta0 <- l.fit$beta0+beta.ols[1]
    beta.bic <- c(beta0[t1$ix[1]],beta[t1$ix[1],])
    beta.aic <- c(beta0[t2$ix[1]],beta[t2$ix[1],])
    bic = min(l.fit$BIC)
    aic = min(l.fit$AIC)
  }
  else {
    beta0 <- l.fit$beta0
    beta.bic <- beta[t1$ix[1],]
    beta.aic <- beta[t2$ix[1],]
    bic = min(l.fit$BIC)
    aic = min(l.fit$AIC)
  }

  obj <- list(beta.ols=beta.ols, beta.bic=beta.bic, beta.aic = beta.aic,
              aic = aic, bic = bic)
  obj

}

###################################
## lars variant for LSA
lars.lsa <- function (Sigma0, b0, intercept,  n,
                      type = c("lasso", "lar"),
                      eps = .Machine$double.eps,max.steps)
{
  type <- match.arg(type)
  b0 = as.vector(b0)
  TYPE <- switch(type, lasso = "LASSO", lar = "LAR")

  #"LASSO"
  n1 <- dim(Sigma0)[1]

  ## handle intercept处理拦截
  if (intercept) {
    a11 <- Sigma0[1,1]
    a12 <- Sigma0[2:n1,1]
    a22 <- Sigma0[2:n1,2:n1]
    Sigma <- a22-outer(a12,a12)/a11
    b <- b0[2:n1]
    beta0 <- crossprod(a12,b)/a11
  }
  else {
    Sigma <- Sigma0
    b <- as.vector(b0)
  }
  Sigma <- diag(abs(b))%*%Sigma%*%diag(abs(b))
  b <- sign(b)

  nm <- dim(Sigma)
  m <- nm[2]
  im <- inactive <- seq(m)

  Cvec <- drop(t(b)%*%Sigma)
  ssy <- sum(Cvec*b)
  if (missing(max.steps))
    max.steps <- 8 * m
  beta <- matrix(0, max.steps + 1, m)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ignores <- NULL

  while ((k < max.steps) & (length(active) < m)) {
    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if (!any(drops)) {
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
        R <- updateR(Sigma[inew, inew], R, drop(Sigma[inew, active]),
                     Gram = TRUE,eps=eps)
        if(attr(R, "rank") == length(active)) {
          ##singularity; back out
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action,  - inew)
        }
        else {
          if(first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
        }
      }
    }
    else action <- -dropid
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1

    if (length(active) >= m) {
      gamhat <- Cmax/A
    }
    else {
      a <- drop(w %*% Sigma[active, -c(active,ignores), drop = FALSE])
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }
    if (type == "lasso") {
      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      # cat('zmin ',zmin, ' gamhat ',gamhat,'\n')
      if (zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
    }
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w

    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w
    Gamrat <- c(Gamrat, gamhat/(Cmax/A))

    arc.length <- c(arc.length, gamhat)
    if (type == "lasso" && any(drops)) {
      dropid <- seq(drops)[drops]
      for (id in rev(dropid)) {
        R <- downdateR(R,id)
      }
      dropid <- active[drops]
      beta[k + 1, dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }

    actions[[k]] <- action
    inactive <- im[-c(active)]
  }
  beta <- beta[seq(k + 1), ]

  dff <- b-t(beta)

  RSS <- diag(t(dff)%*%Sigma%*%dff)

  if(intercept)
    beta <- t(abs(b0[2:n1])*t(beta))
  else
    beta <- t(abs(b0)*t(beta))

  if (intercept) {
    beta0 <- as.vector(beta0)-drop(t(a12)%*%t(beta))/a11
  }
  else {
    beta0 <- rep(0,k+1)
  }
  dof <- apply(abs(beta)>eps,1,sum)
  BIC <- RSS+log(n)*dof
  AIC <- RSS+2*dof
  object <- list(AIC = AIC, BIC = BIC,
                 beta = beta, beta0 = beta0)
  object
}

