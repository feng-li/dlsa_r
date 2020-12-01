# `dlsa`: Distributed Least Squares Approximation

## Introduction

In this work, we develop a distributed least squares approximation (DLSA) method that is able to solve a large family of regression problems (e.g., linear regression, logistic regression, and Cox's model) on a distributed system. By approximating the local objective function using a local quadratic form, we are able to obtain a combined estimator by taking a weighted average of local estimators. The resulting estimator is proved to be statistically as efficient as the global estimator. Moreover, it requires only one round of communication. We further conduct shrinkage estimation based on the DLSA estimation using an adaptive Lasso approach. The solution can be easily obtained by using the LARS algorithm on the master node. It is theoretically shown that the resulting estimator possesses the oracle property and is selection consistent by using a newly designed distributed Bayesian information criterion (DBIC). The finite sample performance and the computational efficiency are further illustrated by an extensive numerical study and an airline dataset. 

- The entire methodology has been implemented in Python for a de-facto standard Spark system available at [https://github.com/feng-li/dlsa](https://github.com/feng-li/dlsa). 
- This package provides the conceptual demo in R.

## Installation

```r
devtools::install_github('feng-li/dlsa_r')
```

## Quick start 

```r
require("dlsa")
source(system.file("scripts/dlsa_test.R", package = "dlsa"), echo = TRUE)
```
## Implement distributed models

Our method is general and can be used for other models, one only need to specify a fitting function and a predict function. See the [source code of logistic function](https://github.com/feng-li/dlsa_r/blob/master/R/logistic.R) for details.


## References

- [Zhu, X.](https://xueningzhu.github.io/), [Li, F.](http://feng.li/), & [Wang, H.](http://hansheng.gsm.pku.edu.cn/), (2019) Least Squares Approximation for a Distributed System. [_Working Paper_](https://arxiv.org/abs/1908.04904).
