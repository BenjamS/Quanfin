tsTrends <- function(in_ts, slope_window = 13, before_window = 5, aft_window = 2, quietly = T)
{
  #Captures index of sustained up and down trends of a time series
  #by identifying periods of consecutive time steps in which slope
  #is positve (for up trends) or negative (for down trends)
  #Since this is based on measuring slope in a rolling window,
  #best to use a smoothed series (via EMA for eg.). Consider scaling too.
  #---------------------
  #Make sure required packages are loaded
  #---------------------
  required_packages <- c("pracma", "ggplot2", "plyr", "tidyr", "quantmod")
  lapply(required_packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  #---------------------
  #Rolling slope
  #---------------------
  #per_ema <- 13
  #in_ts <- EMA(xts_cp[, 1], per_ema)
  #in_ts <- in_ts[-c(1:(per_ema - 1))]
  s_start <- slope_window + 1
  s_end <- length(in_ts)
  x <- c(1:(slope_window + 1))
  dydx_mu <- c()
  dydx_sd <- c()
  t <- 0
  #---------------------
  for(s in s_start:s_end)
  {
    datseg <- in_ts[(s - slope_window):s]
    q <- polyfit(x, datseg, n = 2)
    a2 <- q[1]; a1 <- q[2]; a0 <- q[3]
    this_slope <- 2 * a2 * x + a1
    t <- t + 1
    dydx_mu[t] <- mean(this_slope, na.rm = T)
    dydx_sd[t] <- sd(this_slope, na.rm = T)
  }
  dydx_cv <- dydx_sd / dydx_mu
  mu_dydx_cv <- mean(abs(dydx_cv), na.rm = T)
  ldydx_cv <- log(abs(dydx_cv)) #EMA(log(abs(dydx_cv)), slope_window)
  dydx_sd <- EMA(dydx_sd, slope_window)
  mu_ldydx_cv <- mean(ldydx_cv, na.rm = T)
  mu_dydx_mu <- mean(dydx_mu, na.rm = T)
  #---------------------
  #Identify periods of consecutively positive (for up trends)
  #or negative (for down trends) slope
  #---------------------
  df <- fortify(in_ts)
  colnames(df)[2] <- "ts"
  df$Index <- as.Date(df$Index)
  df$dydxmu <- c(rep(NA, slope_window), dydx_mu)
  df$dydxsd <- c(rep(NA, slope_window), dydx_sd)
  df$ldydxcv <- c(rep(NA, slope_window), ldydx_cv)
  #---------------------
  ind_uptrnd <- which(df$dydxmu > mu_dydx_mu)
  ind_dntrnd <- which(df$dydxmu < mu_dydx_mu)
  df$ts_up <- NA
  df$ts_dn <- NA
  df$ts_up[ind_uptrnd] <- df$ts[ind_uptrnd]
  df$ts_dn[ind_dntrnd] <- df$ts[ind_dntrnd]
  #---------------------
  ind_Selsig <- ind_uptrnd[which(diff(ind_uptrnd) != 1)]
  ind_Buysig <- ind_dntrnd[which(diff(ind_dntrnd) != 1)]
  n_b <- length(ind_Buysig)
  n_s <- length(ind_Selsig)
  if(ind_Buysig[1] > ind_Selsig[1]){ind_Selsig <- ind_Selsig[-1]; n_s <- length(ind_Selsig)}
  if(ind_Buysig[n_b] > ind_Selsig[n_s]){ind_Buysig <- ind_Buysig[-n_b]; n_b <- length(ind_Buysig)}
  #---------------------
  df$RetLong_mag <- NA
  df$RetLong_pct <- NA
  df$timeLong <- NA
  df$RetShrt_mag <- NA
  df$RetShrt_pct <- NA
  df$timeShrt <- NA
  df$RetLong_mag[ind_Buysig] <- df$ts[ind_Selsig] - df$ts[ind_Buysig]
  df$RetLong_pct[ind_Buysig] <- df$RetLong_mag[ind_Buysig] / df$ts[ind_Buysig]
  df$timeLong[ind_Buysig] <- ind_Selsig - ind_Buysig
  df$RetShrt_mag[ind_Selsig[-n_s]] <- -(df$ts[ind_Buysig[-1]] - df$ts[ind_Selsig[-n_s]])
  df$RetShrt_pct[ind_Selsig[-n_s]] <- df$RetShrt_mag[ind_Selsig[-n_s]] / df$ts[ind_Selsig[-n_s]]
  df$timeShrt[ind_Selsig[-n_s]] <- ind_Buysig[-1] - ind_Selsig[-n_s]
  #-----------------------
  #Get date of trend start/finish points
  df$BuysigDate <- NA
  df$BuysigDate[ind_Buysig] <- df$Index[ind_Buysig]
  df$BuysigDate <- as.Date(df$BuysigDate)
  df$SelsigDate <- NA
  df$SelsigDate[ind_Selsig] <- df$Index[ind_Selsig]
  df$SelsigDate <- as.Date(df$SelsigDate)
  #-----------------------
  #Create binary buy/sell var for ML model construction
  #I.e. expand single buy/sell signal date to a window of a few days,
  #with more days before the critical date than after it (before_window > aft_window)
  df$Buysig <- 0
  df$Selsig <- 0
  for(i in 1:length(ind_Buysig)){df$Buysig[(ind_Buysig[i] - before_window):(ind_Buysig[i] + aft_window)] <- 1}
  for(i in 1:length(ind_Selsig)){df$Selsig[(ind_Selsig[i] - before_window):(ind_Selsig[i] + aft_window)] <- 1}
  df$Sig <- "Hold"
  df$Sig[which(df$Selsig == 1)] <- "Sell"
  df$Sig[which(df$Buysig == 1)] <- "Buy"
  #-----------------------
  #Interesting... the up trend start (buy) points coincide with max ldydxcv points
  #And these are all at about the same value (for SPY anyway)
  ind_ind_muldydxcvBuy <- which(df$ldydxcv[ind_Buysig] > 0)
  mu_ldydx_cvBuy <- mean(df$ldydxcv[ind_Buysig[ind_ind_muldydxcvBuy]])
  sd_ldydx_cvBuy <- sd(df$ldydxcv[ind_Buysig[ind_ind_muldydxcvBuy]])
  #---------------------------------
  out_RetMagSum <- sum(df$RetLong_mag, na.rm = T)
  out_RetMagMu <- mean(df$RetLong_mag, na.rm = T)
  out_RetMagSd <- sd(df$RetLong_mag, na.rm = T)
  out_RetPctMu <- mean(df$RetLong_pct, na.rm = T)    
  out_RetPctSd <- sd(df$RetLong_pct, na.rm = T)
  sd_dydx_mu <- sd(dydx_mu, na.rm = T)
  outVars <- data.frame(muLdydxCVbuy = mu_ldydx_cvBuy, 
                        sdLdydxCVbuy = sd_ldydx_cvBuy,
                        mudydx = mu_dydx_mu,
                        sddydx = sd_dydx_mu,
                        RetMagSum = out_RetMagSum,
                        RetMagMu = out_RetMagMu,
                        RetMagSd = out_RetMagSd,
                        RetPctMu = out_RetPctMu,
                        RetPctSd = out_RetPctSd,
                        n_uptrends = n_b)
  #---------------------------------
  #Output df with the original ts and its rolling slope series (mean, sd, cv),
  #buy/sell binary var, etc.
  df_ts <- df[, c("Index", "ts", "dydxmu", "dydxsd", "ldydxcv", "Sig")]
  #Output df with only trend start/finish (buy/sell) point data
  ind_trndBegEnd <- unique(c(which(is.na(df$BuysigDate) == F), which(is.na(df$SelsigDate) == F)))
  df_tsCritPoints <- df[ind_trndBegEnd, 
                      c("Index", "RetLong_pct", "RetShrt_pct", "RetLong_mag",
                        "RetShrt_mag", "timeLong", "timeShrt", "dydxmu",
                        "dydxsd", "ldydxcv")]
  df_tsCritPoints <- df_tsCritPoints[order(df_tsCritPoints$Index), ]
  #=================================
  #2 Graphs
  if(quietly == F){
  #Plot of just the ts with trend start/finish (buy/sell) points overlaid
  df_plot <- df[, c("Index", "BuysigDate", "SelsigDate", "RetLong_pct",
                    "RetShrt_pct", "RetLong_mag", "RetShrt_mag",
                    "ts_up", "ts_dn")]
  df_plot <- df_plot %>% gather(trends, value, ts_up:ts_dn)
  gg <- ggplot(df_plot, aes(x = Index, y = value, group = trends, color = trends))
  gg <- gg + geom_line()
  gg <- gg + geom_point(aes(x = SelsigDate, y = value, size = RetShrt_pct), color = "red")
  gg <- gg + geom_point(aes(x = BuysigDate, y = value, size = RetLong_pct), color = "green")
  print(gg)
  #-----------------------
  #Plot roling slope (mean, sd, and cv) with trend start/finish (buy/sell) points overlaid
  #The dydxmu plot includes the mean
  #The ldydxcv plot includes the mean only of buy points > 0
  df_plot <- df[, c("Index", "BuysigDate", "SelsigDate", "RetLong_pct", 
               "RetShrt_pct", "RetLong_mag", "RetShrt_mag", 
               "dydxmu", "dydxsd", "ldydxcv")]
  df_plot <- df_plot %>% gather(ts_type, value, dydxmu:ldydxcv)
  gg <- ggplot(df_plot, aes(x = Index, y = value)) + geom_line()
  gg <- gg + geom_point(aes(x = SelsigDate, y = value, size = RetShrt_pct), color = "red")
  gg <- gg + geom_point(aes(x = BuysigDate, y = value, size = RetLong_pct), color = "green")
  gg <- gg + facet_wrap(~ts_type, ncol = 1, scales = "free")
  gg <- gg + geom_hline(data = data.frame(yint = mu_dydx_mu, ts_type = "dydxmu"), aes(yintercept = yint), linetype = "dotted")
  gg <- gg + geom_hline(data = data.frame(yint = mu_ldydx_cvBuy, ts_type = "ldydxcv"), aes(yintercept = yint), linetype = "dotted")
  print(gg)
  }
  #=================================
  outlist <- list(outVars, df_tsCritPoints, df_ts)
  return(outlist)
}
