get_data <- function(tickers, from, to) {
  prices <- xts()
  for (tker in tickers) {
    tmp <- Ad(getSymbols(tker, from = from, to = to, auto.assign = FALSE))
    tmp <- na.approx(tmp, na.rm = FALSE)
    prices <- cbind(prices, tmp)
  }
  colnames(prices) <- tickers
  tclass(prices) <- "Date"
  return(prices)
}

get_return_and_split <- function(ratio, prices){
  
  X_log <- CalculateReturns(prices, "log") %>% na.omit()
  X_lin <- CalculateReturns(prices) %>% na.omit()
  
  N <- ncol(X_log)
  T <- nrow(X_log)
  
  T_trn <- round(ratio * T)
  
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):T, ]
  X_lin_trn <- X_lin[1:T_trn, ]
  X_lin_tst <- X_lin[(T_trn+1):T, ]
  
  return(list(
    log_trn = X_log_trn,
    log_tst = X_log_tst,
    lin_trn = X_lin_trn,
    lin_tst = X_lin_tst,
    X_lin   = X_lin,
    X_log   = X_log
  ))
}

assign_variable_global_env <- function(tickers, from, to, ratio){
  # ratio: train test split ratio
  
  prices <- get_data(tickers, from, to)
  temp <- get_return_and_split(ratio, prices)
  
  
  assign("prices", prices, .GlobalEnv)
  # global assignment
  X_log_trn <<- temp$log_trn
  X_log_tst <<- temp$log_tst
  X_lin_trn <<- temp$lin_trn
  X_lin_tst <<- temp$lin_tst
  X_lin     <<- temp$X_lin
  X_log     <<- temp$X_log
}


# estimate from log returns plus transformation
momentsReturnLog2Lin <- function(mu, Sigma) {
  N <- ncol(Sigma)
  mu_ <- exp(mu + 0.5*diag(Sigma)) - 1
  Sigma_ <- matrix(NA, nrow = N, ncol = N)
  for(ii in 1:N)
    for(jj in 1:N)
      Sigma_[ii, jj] <- exp(mu[ii] + mu[jj] + 0.5*(Sigma[ii, ii]+Sigma[jj, jj])) * (exp(Sigma[ii, jj])-1)
  return( list(mu=mu_, Sigma=Sigma_) )
}

# change periodicity of prices and plot MVP weights allocation (currently only support MVP)


whether_simple_or_log_ret <- function(freq, X_lin_trn, X_log_trn, prices){

  prices_ <- xts()
  for (i in 1:ncol(prices)) {
    prices_ <- cbind(prices_, Cl(to.period(prices[, i], period = freq)))
  }
  colnames(prices_) <- colnames(prices)
  tclass(prices_) <- "Date"
  
  # recompute returns
  X_log <- CalculateReturns(prices_, "log")[-1]
  X_lin <- CalculateReturns(prices_)[-1]
  T <- nrow(X_log)  # number of weeks
  
  # split data into training and set data
  T_trn <- round(0.7*T)  # 70% of data
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):T, ]
  X_lin_trn <- X_lin[1:T_trn, ]
  X_lin_tst <- X_lin[(T_trn+1):T, ]
  
  # estimate mu and Sigma
  mu_lin <- colMeans(X_lin_trn)
  Sigma_lin <- cov(X_lin_trn)
  
  mu_log <- colMeans(X_log_trn)
  Sigma_log <- cov(X_log_trn)
  
  tmp <- momentsReturnLog2Lin(mu_log, Sigma_log)
  mu_log_trans <- tmp$mu
  Sigma_log_trans <- tmp$Sigma
  
  # compute MVPs
  w_MVP_lin <- MVP(mu_lin, Sigma_lin)
  w_MVP_log <- MVP(mu_log, Sigma_log)
  w_MVP_log_trans <- MVP(mu_log_trans, Sigma_log_trans)
  w_MVP_all <- cbind("MVP_lin"       = w_MVP_lin, 
                     "MVP_log"       = w_MVP_log, 
                     "MVP_log_trans" = w_MVP_log_trans)
  
  { par(mfrow = c(1, 1), mar = c(5, 5, 4, 8))
    barplot(
      t(w_MVP_all),
      col = rainbow8equal[1:3],
      legend = colnames(w_MVP_all),
      args.legend = list(
        x = "topright",
        bty = "n",
        inset = c(-0.33, 0.3)
      ),
      main = paste0("MVP allocation (data frequency: ", freq, ")"), 
      xlab = "stocks",
      ylab = "dollars",
      beside = TRUE
    )
  }
  
  ret_MVP_all_trn <- xts(X_lin_trn %*% w_MVP_all, index(X_lin_trn))
  SharpeRatio.annualized(ret_MVP_all_trn)
}

