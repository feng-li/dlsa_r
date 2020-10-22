#################################################
###test function

library(MASS)
library(survival)


## the simulator for covariates X
simu.X<-function(N, p, K, iid = T)
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

## simulate the Y according to regression model
simu.Y<-function(X, beta, reg_type = "logistic")
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



##create simulate numbers

beta = c(3, 0, 0, 1.5, 0, 0, 2, 0)
p = length(beta)
N = 10000
set.seed(1234)
X = simu.X(N, p, K = 5,iid=T)
Y = simu.Y(X, beta)
XY = data.frame(cbind(X,Y))
X = data.frame(X)
head(XY)

#####################################
##create the fit function 
logistic <- function(X, Y, init.beta = NULL, iter.max = NULL)
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

###################################
##create the predict function
logistic.pred <- function(X, theta, alpha)
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


#################################
#create the plot_fun
plot.exam <- function(x,y){
  plot = ggplot()+
    geom_point(mapping = aes(x = x, y = y))
  return(plot)
}



###############################
###### test the function
####dlsa.fit
dlsa_fit_result <- dlsa.fit(Y~.,XY,K = 5,lasso_fun = 1,logistic,init.beta = NULL)
dlsa_fit_result

####dlsa.pred
dlsa_pred_result <- dlsa.pred(X, dlsa_fit_result$result$theta,K = 5, pred_fun = logistic.pred, alpha = 0.5)
dlsa_pred_result
dlsa_pred_result$result

##提取出result里面的p值和xbeta
xbeta = list()
p = list()
for (i in 1:5){
  xbeta[[i]] = dlsa_pred_result$result[[i]]$xtheta
  p[[i]] = dlsa_pred_result$result[[i]]$p 
}
xbeta = unlist(xbeta)
p = unlist(p)
xbeta
p

####dlsa.plot
dlsa_plot <- dlsa.plot(plot.exam, roc = 1, Y, p, xbeta, p)
dlsa_plot



