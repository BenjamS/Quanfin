tsTrends <- function(in_ts, slope_window = 13,
                     before_window = 5, aft_window = 2,
                     thresh_RetLong_pct = 0.01,
                     thresh_RetShrt_pct = -0.01, quietly = T)
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
  #in_ts <- xts_cpVWMA[, "SPY"]
  #per_ema <- 13
  #in_ts <- EMA(xts_cp[, 1], per_ema)
  #in_ts <- in_ts[-c(1:(per_ema - 1))]
  n_ts <- length(in_ts)
  s_start <- slope_window + 1
  s_end <- n_ts
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
  sd_dydx_mu <- sd(dydx_mu, na.rm = T)
  cv_dydx_mu <- sd_dydx_mu / mu_dydx_mu
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
  #Initial capture based on slope crossings of the mean slope
  ind_uptrnd <- which(df$dydxmu > mu_dydx_mu)
  ind_dntrnd <- which(df$dydxmu <= mu_dydx_mu)
  #For graphing purposes
  df$ts_up <- NA
  df$ts_dn <- NA
  df$ts_up[ind_uptrnd] <- df$ts[ind_uptrnd]
  df$ts_dn[ind_dntrnd] <- df$ts[ind_dntrnd]
  #---------------------
  #Get trend start/finish points by taking differences
  ind_upFin <- ind_uptrnd[which(diff(ind_uptrnd) != 1)]
  ind_upBeg <- ind_dntrnd[which(diff(ind_dntrnd) != 1)]
  #If necessary, remove start/finish points so that first point is a start
  #and last point is a finish, so that you have a set of complete up/down trends
  #Let's adopt convention that we always start with a buy and end with a sell
  n_upBeg <- length(ind_upBeg)
  n_upFin <- length(ind_upFin)
  n_upBeg_raw <- n_upBeg
  n_upFin_raw <- n_upFin
  if(ind_upBeg[1] > ind_upFin[1]){ind_upFin <- ind_upFin[-1]; n_upFin <- length(ind_upFin)}
  if(ind_upBeg[n_upBeg] > ind_upFin[n_upFin]){ind_upBeg <- ind_upBeg[-n_upBeg]; n_upBeg <- length(ind_upBeg)}
  if(sum(ind_upFin - ind_upBeg < 0) > 0){print("Problem with uptrends")}
  ind_dnBeg <- ind_upFin[-n_upFin]
  ind_dnFin <- ind_upBeg[-1]
  n_dnBeg <- length(ind_dnBeg)
  n_dnFin <- length(ind_dnFin)
  n_dnBeg_raw <- n_dnBeg
  n_dnFin_raw <- n_dnFin
  if(sum(ind_dnFin - ind_dnBeg < 0) > 0){print("Problem with downtrends")}
  #Get value of trend increase/decrease (as magnitude and percentage)
  #and trend durations
  RetLong_mag <- df$ts[ind_upFin] - df$ts[ind_upBeg]
  RetLong_pct <- RetLong_mag / df$ts[ind_upBeg]
  timeLong <- ind_upFin - ind_upBeg
  RetShrt_mag <- df$ts[ind_dnFin] - df$ts[ind_dnBeg]
  RetShrt_pct <- RetShrt_mag / df$ts[ind_dnBeg]
  timeShrt <- ind_dnFin - ind_dnBeg
  df$Ret_mag <- NA
  df$Ret_pct <- NA
  df$trendDuration <- NA
  df$Ret_mag[ind_upBeg] <- RetLong_mag
  df$Ret_pct[ind_upBeg] <- RetLong_pct
  df$Ret_mag[ind_dnBeg] <- RetShrt_mag
  df$Ret_pct[ind_dnBeg] <- RetShrt_pct
  df$trendDuration[ind_upBeg] <- timeLong
  df$trendDuration[ind_dnBeg] <- timeShrt
  #---------------------
  #Flag false start/finishes (but don't drop)
  ind_ind_upBeg_keep <- which(RetLong_pct >= thresh_RetLong_pct)
  ind_upBeg_keep <- ind_upBeg[ind_ind_upBeg_keep]
  ind_upFin_keep <- ind_upFin[ind_ind_upBeg_keep]
  n_upBeg_keep <- length(ind_upBeg_keep)
  ind_ind_upBeg_false <- setdiff(1:n_upBeg, ind_ind_upBeg_keep)
  ind_upBeg_false <- ind_upBeg[ind_ind_upBeg_false]
  ind_upFin_false <- ind_upFin[ind_ind_upBeg_false]
  n_upBeg_false <- length(ind_upBeg_false)
  #ind_upFin_keep - ind_upBeg_keep
  #ind_upFin_false - ind_upBeg_false
  RetLong_mag_keep <- RetLong_mag[ind_ind_upBeg_keep]
  RetLong_pct_keep <- RetLong_pct[ind_ind_upBeg_keep]
  timeLong_keep <- ind_upFin_keep - ind_upBeg_keep
  sum_RetLong_mag_keep <- sum(RetLong_mag_keep)
  mu_RetLong_mag_keep <- mean(RetLong_mag_keep)
  sd_RetLong_mag_keep <- sd(RetLong_mag_keep)
  cv_RetLong_mag_keep <- sd_RetLong_mag_keep / mu_RetLong_mag_keep
  mu_RetLong_pct_keep <- mean(RetLong_pct_keep)
  sd_RetLong_pct_keep <- sd(RetLong_pct_keep)
  cv_RetLong_pct_keep <- sd_RetLong_pct_keep / mu_RetLong_pct_keep
  #--------------
  ind_ind_dnBeg_keep <- which(RetShrt_pct <= thresh_RetShrt_pct)
  ind_dnBeg_keep <- ind_dnBeg[ind_ind_dnBeg_keep]
  ind_dnFin_keep <- ind_dnFin[ind_ind_dnBeg_keep]
  n_dnBeg_keep <- length(ind_dnBeg_keep)
  n_dnFin_keep <- n_dnBeg_keep
  ind_ind_dnBeg_false <- setdiff(1:n_dnBeg, ind_ind_dnBeg_keep)
  ind_dnBeg_false <- ind_dnBeg[ind_ind_dnBeg_false]
  ind_dnFin_false <- ind_dnFin[ind_ind_dnBeg_false]
  n_dnBeg_false <- length(ind_dnBeg_false)
  n_dnFin_false <- n_dnBeg_false
  #ind_dnFin_keep - ind_dnBeg_keep
  #ind_dnFin_false - ind_dnBeg_false
  RetShrt_mag_keep <- RetShrt_mag[ind_ind_dnBeg_keep]
  RetShrt_pct_keep <- RetShrt_pct[ind_ind_dnBeg_keep]
  timeShrt_keep <- ind_dnFin_keep - ind_dnBeg_keep
  sum_RetShrt_mag_keep <- sum(RetShrt_mag_keep)
  mu_RetShrt_mag_keep <- mean(RetShrt_mag_keep)
  sd_RetShrt_mag_keep <- sd(RetShrt_mag_keep)
  cv_RetShrt_mag_keep <- sd_RetShrt_mag_keep / mu_RetShrt_mag_keep
  mu_RetShrt_pct_keep <- mean(RetShrt_pct_keep)
  sd_RetShrt_pct_keep <- sd(RetShrt_pct_keep)
  cv_RetShrt_pct_keep <- sd_RetShrt_pct_keep / mu_RetShrt_pct_keep
  #--
  # df$Ret_pct[ind_dnBeg_keep]
  gg <- ggplot(df, aes(x = Index, y = ts)) + geom_line()
  # gg <- ggplot(df, aes(x = Index, y = dydxmu)) + geom_line()
  # gg <- gg + geom_vline(xintercept = df$Index[ind_upBeg_keep], color = "green")
  # gg <- gg + geom_vline(xintercept = df$Index[ind_upFin_keep], color = "red")
  gg <- gg + geom_vline(xintercept = df$Index[ind_dnBeg_keep], color = "blue")
  gg <- gg + geom_vline(xintercept = df$Index[ind_dnFin_keep], color = "orange")
  #gg <- gg + geom_hline(yintercept = mu_dydx_mu, color = "purple")
  gg
  #---------------
  #Signals
  #Create signal date var
  df$SigDate <- NA
  df$SigDate[ind_dnBeg] <- df$Index[ind_dnBeg]
  df$SigDate[ind_upBeg] <- df$Index[ind_upBeg]
  df$SigDate[ind_dnFin] <- df$Index[ind_dnFin]
  df$SigDate[ind_upFin] <- df$Index[ind_upFin]
  #Create signal var over "window of opportunity" determined by params
  #before_window and aft_window. before_window should be greater than aft_window.
  #----------
  #Long Sig
  #----------
  #Hold
  df$SigLong <- "Hold"
  #Buy/Sell
  for(i in 1:n_upFin_keep){df$SigLong[(ind_upFin_keep[i] - before_window):(ind_upFin_keep[i] + aft_window)] <- "Sell"}
  for(i in 1:n_upBeg_keep){df$SigLong[(ind_upBeg_keep[i] - before_window):(ind_upBeg_keep[i] + aft_window)] <- "Buy"}
  #Overwrite buy/sell sigs outside pct return thresholds as false signals
  for(i in 1:n_upFin_false){df$SigLong[(ind_upFin_false[i] - before_window):(ind_upFin_false[i] + aft_window)] <- "False Sell"}
  for(i in 1:n_upBeg_false){df$SigLong[(ind_upBeg_false[i] - before_window):(ind_upBeg_false[i] + aft_window)] <- "False Buy"}
  #----------
  #Short Sig
  #----------
  #Hold
  df$SigShort <- "Hold"
  for(i in 1:n_dnFin_keep){df$SigShort[(ind_dnFin_keep[i] - before_window):(ind_dnFin_keep[i] + aft_window)] <- "Short"}
  for(i in 1:n_dnBeg_keep){df$SigShort[(ind_dnBeg_keep[i] - before_window):(ind_dnBeg_keep[i] + aft_window)] <- "Cover"}
  #Overwrite buy/sell sigs outside thresholds as false signals
  for(i in 1:n_dnFin_false){df$SigShort[(ind_dnFin_false[i] - before_window):(ind_dnFin_false[i] + aft_window)] <- "False Short"}
  for(i in 1:n_dnBeg_false){df$SigShort[(ind_dnBeg_false[i] - before_window):(ind_dnBeg_false[i] + aft_window)] <- "False Cover"}
  #------------------------
  #Trim ts to first uptrend start point and last downtrend stop point
  #so that ML training not messed up
  ind_tsStart <- min(ind_upBeg[1], ind_dnBeg[1]) - before_window
  ind_tsFinish <- max(ind_upFin[n_upFin], ind_dnFin[n_dnFin]) + aft_window
  df <- df[ind_tsStart:ind_tsFinish, ]
  n_ts_trimd <- nrow(df)
  #-----------------------
  #Interesting... the up trend start (buy) points coincide with max ldydxcv points
  #And these are all at about the same value (for SPY anyway)
  ind_ind_muldydxcv_up <- which(df$ldydxcv[ind_upBeg_keep] > 0)
  mu_ldydx_cv_up <- mean(df$ldydxcv[ind_upBeg_keep[ind_ind_muldydxcv_up]])
  sd_ldydx_cv_up <- sd(df$ldydxcv[ind_upBeg[ind_ind_muldydxcv_up]])
  #
  mu_dydx_mu_up <- mean(df$dydxmu[ind_upBeg_keep])
  sd_dydx_mu_up <- sd(df$dydxmu[ind_upBeg_keep])
  cv_dydx_mu_up <- sd_dydx_mu_up / mu_dydx_mu_up
  #
  mu_dydx_mu_dn <- mean(df$dydxmu[ind_dnBeg_keep])
  sd_dydx_mu_dn <- sd(df$dydxmu[ind_dnBeg_keep])
  cv_dydx_mu_dn <- sd_dydx_mu_dn / mu_dydx_mu_dn
  #---------------------------------
  
  outVars <- data.frame(
    n_uptrends_kept = n_upBeg_keep,
    n_downtrends_kept = n_dnBeg_keep,
    n_uptrends_false = n_upBeg_false,
    n_downtrends_false = n_dnBeg_false,

    mudydx = mu_dydx_mu,
    cvdydx = cv_dydx_mu,
    mudydxup = mu_dydx_mu_up,
    cvdydxup = cv_dydx_mu_up,
    mudydxdn = mu_dydx_mu_dn,
    cvdydxdn = cv_dydx_mu_dn,
    
    muLdydxCVup = mu_ldydx_cv_up, 
    sdLdydxCVup = sd_ldydx_cv_up,

    RetMagUp_sum = sum_RetLong_mag_keep,
    RetMagUp_mu = mu_RetLong_mag_keep,
    RetMagUp_cv = cv_RetLong_mag_keep,
    RetPctUp_mu = mu_RetLong_pct_keep,
    RetPctUp_cv = cv_RetLong_pct_keep,
    
    RetMagDown_sum = sum_RetShrt_mag_keep,
    RetMagDown_mu = mu_RetShrt_mag_keep,
    RetMagDown_cv = cv_RetShrt_mag_keep,
    RetPctDown_mu = mu_RetShrt_pct_keep,
    RetPctDown_cv = cv_RetShrt_pct_keep
  )
  
  #---------------------------------
  cat("Raw length of time series: ", n_ts, "\n",
      "Trimmed length of time series: ", n_ts_trimd, "\n",
      "Num. raw uptrend signals: ", n_upBeg_raw, "\n",
      "Num. raw downtrend signals: ", n_upFin_raw, "\n",
      "Num. raw uptrend signals after bookend adjust: ", n_upBeg, "\n",
      "Num. raw downtrend signals after bookend adjust: ", n_dnBeg, "\n",
      "Num. false uptrend signals: ", n_upBeg_false, "\n",
      "Num. false downtrend signals: ", n_dnBeg_false, "\n",
      "Num. uptrend signals within pct change thresholds (kept): ", n_upBeg_keep, "\n",
      "Num. downtrend signals within pct change thresholds (kept): ", n_dnBeg_keep, "\n",
      "Sum of all kept uptrends: ", sum_RetLong_mag_keep, "\n",
      "Sum of all kept downtrends: ", sum_RetShrt_mag_keep, "\n",
      "Mean change kept uptrends: ", mu_RetLong_mag_keep, "\n",
      "Mean change kept downtrends: ", mu_RetShrt_mag_keep, "\n",
      "Mean pct. change kept uptrends: ", mu_RetLong_pct_keep, "\n",
      "Mean pct. change kept downtrends: ", mu_RetShrt_pct_keep, "\n",
      "CV change kept uptrends: ", cv_RetLong_mag_keep, "\n",
      "CV change kept downtrends: ", cv_RetShrt_mag_keep, "\n",
      "Mean slope: ", mu_dydx_mu, "\n",
      "CV slope: ", cv_dydx_mu, "\n",
      "Mean slope at kept uptrend starts: ", mu_dydx_mu_up, "\n",
      "CV slope at kept uptrend starts: ", cv_dydx_mu_up, "\n",
      "Mean slope at kept downtrend starts: ", mu_dydx_mu_dn, "\n",
      "CV slope at kept downtrend starts: ", cv_dydx_mu_dn
  )
  
  #---------------------------------
  #Output df with the original ts and its rolling slope series (mean, sd, cv),
  #buy/sell binary var, etc.
  df_ts <- df[, c("Index", "ts", "dydxmu", "dydxsd", "ldydxcv", "SigDate", "SigLong", "SigShort")]
  #Output df with only trend start/finish (buy/sell) point data
  ind_trndBegEnd <- unique(c(which(is.na(df$BuysigDate) == F), which(is.na(df$SelsigDate) == F)))
  df_tsCritPoints <- df[ind_trndBegEnd, 
                        c("Index", "SigDate", "SigLong", "SigShort", "ts", "dydxmu", "dydxsd", "ldydxcv",
                          "Ret_mag", "Ret_pct", "trendDuration")]
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
    gg <- gg + geom_hline(data = data.frame(yint = mu_ldydx_cv_up, ts_type = "ldydxcv"), aes(yintercept = yint), linetype = "dotted")
    print(gg)
  }
  #=================================
  outlist <- list(outVars, df_tsCritPoints, df_ts)
  return(outlist)
}
