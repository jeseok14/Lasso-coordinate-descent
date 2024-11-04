
### library 및 함수 불러오기
source('functions.R')
library(tidyverse)
library(glmnet)

## lasso 알고리즘 확인
# 정규분포를 따르는 설명변수 생성 후 결과과 확인
set.seed(10)
x <- matrix(runif(2000,-1,1),nrow = 200)
dim(x)
error <- rnorm(200,0,0.1)
yi <- x[,1]+x[,2]-x[,3] + error
lambda <- 2^(-10:10)
result = lasso(yi,x,0.5)
result

## lambda별로 lasso 모델 적용 결과 확인
results <- make_results(x,yi,lambda)
path_plot(results, lambda)

## true beta = (0,1,1,-1,0,...,0)^T 라고 하고 lambda별로 추정된 lasso solution 계산
## MSE(lambda)를 세로축, log(lambda)를 가로축으로 선으로 이어서 그림림
x <- matrix(runif(2000,-1,1),nrow = 200)
error <- rnorm(200,0,0.1)
real_beta <- rep(0,ncol(x))
real_beta[c(2,3,4)] <- c(1,1,-1)
real_beta
yi <- as.vector(x%*%real_beta + error)
lambda <- 2^(-10:10)
results_Q3 <- make_results(x,yi,lambda)
mat_3 <- results_Q3[[1]]$beta
for(i in 2:length(results_Q3)){
  mat_3 <- cbind(mat_3,results_Q3[[i]]$beta)
}
# 각 lambda에 대한 \hat \beta_\lambda
print(mat_3)

MSE_lambda <- apply(mat_3, 2, function(x){ sqrt(sum((x - c(0,real_beta))^2))})
data <- data.frame(log_lambda = log10(lambda), MSE_lambda = MSE_lambda)
data

print(min(data$MSE_lambda))
print(which.min(data$MSE_lambda))
log_lambda_star <- data$log_lambda[which.min(data$MSE_lambda)]
print(10^log_lambda_star)
ggplot(data, aes(log_lambda, MSE_lambda)) + geom_point() + geom_line(linetype = 'dashed') + 
  geom_vline(xintercept = log_lambda_star,linetype = 'dashed', col = 'blue')
data[data$log_lambda == log_lambda_star,]

## solution path 그리기

print(10^(log_lambda_star))
path_plot_Q4(results_Q3,lambda,log_lambda_star)

## 동일한 자료에 대해서 glmnet 패키지를 통한 lasso 수행(glmnet - cv.glmnet())
lambda_star <- 10^(log_lambda_star)
star.result = lasso(yi,x,lambda_star)
star.coef <- star.result$beta
print(star.coef)
cv.out <- cv.glmnet(x, yi, alpha = 1)
bestlam <- cv.out$lambda.min
out <- glmnet(x,yi,alpha = 1)
best_lasso <- glmnet(x, yi, alpha = 1, lambda = bestlam)
lasso.coef <- drop( predict(out, type = "coefficients", s = bestlam))
print(lasso.coef)
## 차이 계산
print(sqrt(sum((lasso.coef - star.coef)^2)))
## 차이가 실제로 작음.

## 10-fold-cross-validation
result_cv <- cvlasso(x,yi,lambda, 10)
mean_mse_values <- apply(result_cv$mse_value, 2, mean)
data_cv <- data.frame(lambda, mean_mse_values)
best_lambda_cv <- lambda[which.min(mean_mse_values)]
data_cv[data_cv$lambda == best_lambda_cv,]

