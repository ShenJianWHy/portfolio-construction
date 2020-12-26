# Equal weighted portfolios
EWP <- function(logRet) {
  N <- ncol(logRet)
  w_EWP <- rep(1/N, N)
  names(w_EWP) <- colnames(logRet)
  return(w_EWP)
}


quintileP <- function(logRet){ #TODO criteria="mu_sigma"
  mu <- colMeans(logRet)
  Sigma <- cov(logRet) 
  N <- ncol(logRet)
  
  # find indices of sorted stocks
  i1 <- sort(mu, decreasing = TRUE, index.return = TRUE)$ix
  i2 <- sort(mu/diag(Sigma), decreasing = TRUE, index.return = TRUE)$ix
  i3 <- sort(mu/sqrt(diag(Sigma)), decreasing = TRUE, index.return = TRUE)$ix
  
  w_QuintP_1 <- w_QuintP_2 <- w_QuintP_3 <- rep(0, N)
  w_QuintP_1[i1[1:round(N/5)]] <- 1/round(N/5)
  w_QuintP_2[i2[1:round(N/5)]] <- 1/round(N/5)
  w_QuintP_3[i3[1:round(N/5)]] <- 1/round(N/5)
  w_QuintP <- cbind("QuintP (mu)"        = w_QuintP_1, 
                    "QuintP (mu/sigma2)" = w_QuintP_2, 
                    "QuintP (mu/sigma)"  = w_QuintP_3)
  rownames(w_QuintP) <- colnames(logRet)
  
  return(w_QuintP)
}

# GMRP global maximum return portfolio
GMRP <- function(mu) {
  w <- Variable(length(mu))
  prob <- Problem(Maximize(t(w) %*% mu), 
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- names(mu)
  return(w)
}
GMVP <- function(Sigma, long_only=FALSE) {
  w <- Variable(nrow(Sigma))
  if(long_only) {
    constraints = list(w >= 0, sum(w) == 1)
  } else {
    constraints = list(sum(w) == 1)
  }
  prob <- Problem(Minimize(quad_form(w, Sigma)), 
                  constraints = constraints)
  result <- solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- colnames(Sigma)
  return(w)
}

MVP <- function(mu, Sigma, lmd = 0.5, long_only=FALSE) {
  
  w <- Variable(nrow(Sigma))
  if(long_only) {
    constraints = list(w >= 0, sum(w) == 1)
  } else {
    constraints = list(sum(w) == 1)
  }
  prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                  constraints = constraints)
  result <- solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- colnames(Sigma)
  return(w)
}

MSRP <- function(mu, Sigma) { #TODO Longonly
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1))
  result <- solve(prob)
  w <- as.vector(result$getValue(w_)/sum(result$getValue(w_)))
  names(w) <- colnames(Sigma)
  return(w)
}

# Inverse volatility portfolio (IVP)
IVP <- function(Sigma) {
  sigma <- sqrt(diag(Sigma))
  w <- 1/sigma
  w <- w/sum(w)
  return(w)
}

# Most diversified portfolio (MDP)
MDP <- function(Sigma) {
  MSRP(mu = sqrt(diag(Sigma)), Sigma=Sigma)
}

# Maximum decorrelation portfolio (MDCP)
MDCP <- function(Sigma) {
  C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
  colnames(C) <- colnames(Sigma)
  return(GMVP(Sigma = C))
}

portfolioDR <- function(X, lmd = 0.5, alpha = 2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  mu <- colMeans(X)
  w <- Variable(N)
  prob <- Problem(Maximize(t(w) %*% mu - (lmd/T) * sum(pos(t(mu) %*% w - X %*% w))^alpha),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

portolioCVaR <- function(X, lmd = 0.5, alpha = 0.95) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  z <- Variable(T)
  zeta <- Variable(1)
  # problem
  prob <- Problem(Maximize(t(w) %*% mu - lmd*zeta - (lmd/(T*(1-alpha))) * sum(z)),
                  constraints = list(z >= 0, z >= -X %*% w - zeta,
                                     w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

portfolioMaxDD <- function(X, c = 0.2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  u <- Variable(T)
  # problem
  # FIXME
  # https://palomar.home.ece.ust.hk/MAFS6010R_lectures/slides_CVaR_portfolio.html#58
  prob <- Problem(Maximize(t(w) %*% mu),
                  constraints = list(w >= 0, sum(w) == 1,
                                     u <= X_cum %*% w + c,
                                     u >= X_cum %*% w,
                                     u[-1] >= u[-T]))
  result <- solve(prob)
  
  tryCatch(
    return(as.vector(result$getValue(w))), 
    error = function(e){message(
      "argument c might be too low so solver can not find solution")
      return(NULL)
    }
  )
}


portfolioAveDD <- function(X, c = 0.2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  u <- Variable(T)
  # problem
  prob <- Problem(Maximize(t(w) %*% mu),
                  constraints = list(w >= 0, sum(w) == 1,
                                     mean(u) <= mean(X_cum %*% w) + c,
                                     u >= X_cum %*% w,
                                     u[-1] >= u[-T]))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}


portfolioCDaR <- function(X, c = 0.1, alpha = 0.95) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  z <- Variable(T)
  zeta <- Variable(1)
  u <- Variable(T)
  # problem
  prob <- Problem(Maximize(t(w) %*% mu),
                  constraints = list(w >= 0, sum(w) == 1,
                                     zeta + (1/(T*(1-alpha))) * sum(z) <= c,
                                     z >= 0, z >= u - X_cum %*% w - zeta,
                                     u >= X_cum %*% w,
                                     u[-1] >= u[-T]))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

