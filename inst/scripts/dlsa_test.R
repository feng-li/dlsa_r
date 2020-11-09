
## Load the data: create simulate numbers
beta = c(3, 0, 0, 1.5, 0, 0, 2, 0)
p = length(beta)
N = 10000
K = 5
set.seed(1234)
X = logistic_simuX(N, p, K,iid=T)
Y = logistic_simuY(X, beta)
XY = data.frame(cbind(X,Y))
#X = data.frame(X)

### extract the p value and xbeta form result

extrc <- function(res){
  xbeta = list()
  p = list()
  k = length(res$result)
  for (i in 1:k){
    xbeta[[i]] = res$result[[i]]$xtheta
    p[[i]] = res$result[[i]]$p
  }
  xbeta = unlist(xbeta)
  p = unlist(p)
  return(list(xbeta = xbeta, p = p))
}

### the final result :plot the roc curve

dlsa.fit(Y~., XY, K, lasso_fun = 1, logistic_fit, init.beta = NULL)$result$theta %>% 
  dlsa.pred(X, ., K, pred_fun = logistic_pred, alpha = 0.5) %>%
  extrc() %>% {res = .
               dlsa.plot(logistic_plot, roc = 1, Y, res$p, res$xbeta, res$p)$roc}

### the final result : Yhat

dlsa.fit(Y~., XY, K, lasso_fun = 1, logistic_fit, init.beta = NULL)$result$theta %>% 
  dlsa.pred(X, ., K, pred_fun = logistic_pred, alpha = 0.5) %>% .$Yhat


