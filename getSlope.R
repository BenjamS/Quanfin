getSlope <- function(in_ts, slope_per = 2, per_ema = NULL, Programatic = T)
{
  #---------------------
  #Make sure required packages are loaded
  #---------------------
  required_packages <- c("pracma", "quantmod")
  lapply(required_packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  #---------------------
  #Rolling slope
  #---------------------
  #in_ts <- xts_cpVWMA[-n_ts, "SPY"]
  #in_ts <- EMA(xts_cp[, 1], per_ema)
  #in_ts <- in_ts[-c(1:(per_ema - 1))]
  if(is.null(per_ema) == F){
    ts <- EMA(in_ts, per_ema)
  }else{
    ts <- in_ts
  }
  n_ts <- length(ts)
  s_start <- slope_per + 1
  s_end <- n_ts
  x <- c(1:(slope_per + 1))
  dydx_mu <- c()
  dydx_sd <- c()
  t <- 0
  #---------------------
  for(s in s_start:s_end)
  {
    datseg <- ts[(s - slope_per):s]
    q <- polyfit(x, datseg, n = 2)
    a2 <- q[1]; a1 <- q[2]; a0 <- q[3]
    this_slope <- 2 * a2 * x + a1
    t <- t + 1
    dydx_mu[t] <- mean(this_slope, na.rm = T)
    dydx_sd[t] <- sd(this_slope, na.rm = T)
  }
  dydx_cv <- dydx_sd / dydx_mu
  mu_dydx_cv <- mean(abs(dydx_cv), na.rm = T)
  ldydx_cv <- log(abs(dydx_cv)) #EMA(log(abs(dydx_cv)), slope_per)
  dydx_sd <- EMA(dydx_sd, slope_per)
  mu_ldydx_cv <- mean(ldydx_cv, na.rm = T)
  mu_dydx_mu <- mean(dydx_mu, na.rm = T)
  sd_dydx_mu <- sd(dydx_mu, na.rm = T)
  cv_dydx_mu <- sd_dydx_mu / mu_dydx_mu
  #-----------
  if(Programatic == T){
    return(dydx_mu)
  }else{
    outlist <- list(dydx_mu, dydx_cv, dydx_sd)
    return()
  }
}
