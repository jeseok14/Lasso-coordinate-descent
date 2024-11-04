soft <- function(t, l){
  if(t > l){
    t - l
  }
  else if(abs(t) <= l){
    0
  }
  else{
    t + l
  }
}

objective_function <- function(y,X,a,l){
  result <- sum((y - X %*% a)^2) + l * sum(abs(a))
  return(result)
}

lasso <- function(y,X,lambda){
  sy <- sqrt(sum((y - mean(y))^2))
  bar_xj <- colMeans(X)
  bar_matrix <- matrix(bar_xj,
                       nrow = nrow(X), ncol = ncol(X),
                       byrow = T)
  X_sub <- X - bar_matrix
  X2_sub <- X_sub^2
  sum_X2_sub <- colSums(X2_sub)
  sx <- sqrt(sum_X2_sub)
  ys <- (y - mean(y))/sy
  Xs <- t(t(X_sub) / sx)
  del <- 10^{-4}
  n <- nrow(Xs)
  p <- ncol(Xs)
  alpha <- matrix(rep(0,p),ncol = p)
  Q <- objective_function(ys,Xs,as.vector(alpha),lambda)
  m <- 0
  while(TRUE){
    m = m + 1
    alpha_n <- as.vector(alpha[m,])
    next_alpha <- alpha_n
    for(j in 1:p){
      xj <- Xs[,j] # n x 1
      Xm <- Xs[,-j] # n x (p-1)
      alpha_m <- alpha_n[-j] # (p - 1)
      t <- xj %*% (ys - Xm %*% alpha_m)
      next_alpha[j] <- soft(t,lambda/2)
    }
    alpha_n <- next_alpha
    alpha <- rbind(alpha, alpha_n)
    Qm <- objective_function(ys, Xs, alpha_n, lambda)
    Q <- append(Q, Qm)
    if(m >= 2){
      if(abs(Q[m-1] - Q[m]) / Q[m-1] < del){
        beta <- (alpha_n * sy)/sx
        beta_0 <- mean(y) - sum(beta * bar_xj)
        break
      }
    }
  }
  return(list(beta = c(beta_0,beta),beta_ = beta,number = m,alpha <- alpha,Q = Q))
}

make_results <- function(x,y,l_vec){
  results <- list()
  for(i in 1:length(l_vec)){
    l <- l_vec[i]
    result <- lasso(y,x,l)
    results[[i]] <- result
  }
  return(results)
}
path_plot <- function(results, lambda){
  mat <- results[[1]]$beta_
  for(i in 2:length(results)){
    mat <- cbind(mat, results[[i]]$beta_)
  }
  tmp <- as.data.frame(as.matrix(mat))
  row.names(tmp) <- paste0('V',1:ncol(x))
  colnames(tmp) <- paste0('s',1:length(lambda))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp)
  tmp$variable <- as.numeric(gsub('s','',tmp$variable))
  tmp$lambda <- lambda[tmp$variable]
  tmp$norm2 <- apply(abs(mat),2, function(x){sqrt(sum(x^2))})[tmp$variable]
  p1 <- ggplot(tmp, aes(norm2, value, color = coef, linetype = coef)) +
    geom_line() + xlab('L2 norm') + ylab('Coefficients') +
    theme_bw()
  return(p1)
}


path_plot_Q4 <- function(results, lambda,log_lambda_star){
  mat <- results[[1]]$beta_
  for(i in 2:length(results)){
    mat <- cbind(mat, results[[i]]$beta_)
  }
  tmp <- as.data.frame(as.matrix(mat))
  row.names(tmp) <- paste0('V',1:ncol(x))
  colnames(tmp) <- paste0('s',1:length(lambda))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp)
  tmp$variable <- as.numeric(gsub('s','',tmp$variable))
  tmp$lambda <- lambda[tmp$variable]
  tmp$norm2 <- apply(abs(mat),2, function(x){sqrt(sum(x^2))})[tmp$variable]
  beta_star <- tmp[tmp$lambda == 10^(log_lambda_star),'value']
  print(beta_star)
  beta_star_norm <- sqrt(sum(beta_star^2))
  print(beta_star_norm)
  p1 <- ggplot(tmp, aes(norm2, value, color = coef, linetype = coef)) +
    geom_line() + geom_vline(xintercept = beta_star_norm, linetype = 'dashed', col ='blue') +
    xlab('L2 norm') + ylab('Coefficinets') +
    theme_bw()
  return(p1)
}


cvlasso <- function(x,y,l_vec,k){
  x_index = sample(1:nrow(x),size = nrow(x))
  split_index = split(x_index, 1:k)
  
  mse_values <- matrix(0, ncol = length(l_vec), nrow = k)
  for(i in 1:k){
    valid_x = x[split_index[[i]],]
    train_x = x[-split_index[[i]],]
    valid_yi = yi[split_index[[i]]]
    train_yi = yi[-split_index[[i]]]
    for(j in 1:length(l_vec)){
      result <- lasso(train_yi, train_x, l_vec[j])
      b <- result$beta
      y_pred <- as.vector(cbind(1, valid_x) %*% b)
      mse_values[i,j] <- mean((y_pred - valid_yi)^2)
    }
  }
  return(list(mse_value = mse_values,k = k))
}

