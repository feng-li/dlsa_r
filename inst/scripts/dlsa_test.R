#################################################
###test function

library(MASS)
library(survival)

## Load the data: create simulate numbers
beta = c(3, 0, 0, 1.5, 0, 0, 2, 0)
p = length(beta)
N = 10000
set.seed(1234)
X = logistic_simuX(N, p, K = 5,iid=T)
Y = logistic_simuY(X, beta)
XY = data.frame(cbind(X,Y))
X = data.frame(X)
head(XY)


####dlsa.fit
dlsa_fit_result <- dlsa.fit(Y~.,XY,K = 5,lasso_fun = 1,logistic_fit,init.beta = NULL)
dlsa_fit_result

####dlsa.pred
dlsa_pred_result <- dlsa.pred(X, dlsa_fit_result$result$theta,K = 5, pred_fun = logistic_pred, alpha = 0.5)
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
dlsa_plot <- dlsa.plot(plot_fun, roc = 1, Y, p, xbeta, p)
dlsa_plot
