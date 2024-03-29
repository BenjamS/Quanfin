tsTrends <- function(in_ts, slope_window = 5,
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
  required_packages <- c("pracma", "ggplot2", "plyr", "tidyr", "quantmod", "ggpubr")
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
  ind_dnFin <- ind_upBeg
  ind_dnBeg <- ind_upFin
  #If necessary, remove start/finish points so that first point is a start
  #and last point is a finish, so that you have a set of complete up/down trends
  #Let's adopt convention that we always start with a buy and end with a sell
  n_upBeg <- length(ind_upBeg)
  n_upFin <- length(ind_upFin)
  n_upBeg_raw <- n_upBeg
  n_upFin_raw <- n_upFin
  n_dnBeg <- length(ind_dnBeg)
  n_dnFin <- length(ind_dnFin)
  n_dnBeg_raw <- n_dnBeg
  n_dnFin_raw <- n_dnFin
  if(ind_upBeg[1] > ind_upFin[1]){ind_upFin <- ind_upFin[-1]; n_upFin <- length(ind_upFin)}
  if(ind_upBeg[n_upBeg] > ind_upFin[n_upFin]){ind_upBeg <- ind_upBeg[-n_upBeg]; n_upBeg <- length(ind_upBeg)}
  if(sum(ind_upFin - ind_upBeg < 0) > 0){print("Problem with uptrends")}
  if(ind_dnBeg[1] > ind_dnFin[1]){ind_dnFin <- ind_dnFin[-1]; n_dnFin <- length(ind_dnFin)}
  if(ind_dnBeg[n_dnBeg] > ind_dnFin[n_dnFin]){ind_dnBeg <- ind_dnBeg[-n_dnBeg]; n_dnBeg <- length(ind_dnBeg)}
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
  n_upFin_keep <- n_upBeg_keep
  ind_ind_upBeg_false <- setdiff(1:n_upBeg, ind_ind_upBeg_keep)
  ind_upBeg_false <- ind_upBeg[ind_ind_upBeg_false]
  ind_upFin_false <- ind_upFin[ind_ind_upBeg_false]
  n_upBeg_false <- length(ind_upBeg_false)
  n_upFin_false <- n_upBeg_false
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
  #-----------------------
  #--
  # df$Ret_pct[ind_dnBeg_keep]
  # gg <- ggplot(df, aes(x = Index, y = ts)) + geom_line()
  # # gg <- ggplot(df, aes(x = Index, y = dydxmu)) + geom_line()
  # gg <- gg + geom_vline(xintercept = df$Index[ind_upBeg_keep], color = "green")
  # gg <- gg + geom_vline(xintercept = df$Index[ind_upFin_keep], color = "red")
  # # gg <- gg + geom_vline(xintercept = df$Index[ind_dnBeg_keep], color = "blue")
  # # gg <- gg + geom_vline(xintercept = df$Index[ind_dnFin_keep], color = "orange")
  # #gg <- gg + geom_hline(yintercept = mu_dydx_mu, color = "purple")
  # gg
  #---------------
  #Signals
  #Create signal type var for subsetting
  df$SigUp <- NA
  df$SigUp[ind_upBeg_keep] <- "Uptrend Start"
  df$SigUp[ind_upFin_keep] <- "Uptrend Stop"
  df$SigUp[ind_upBeg_false] <- "Uptrend False Start"
  df$SigUp[ind_upFin_false] <- "Uptrend False Stop"
  df$SigDown <- NA
  df$SigDown[ind_dnBeg_keep] <- "Downtrend Start"
  df$SigDown[ind_dnFin_keep] <- "Downtrend Stop"
  df$SigDown[ind_dnBeg_false] <- "Downtrend False Start"
  df$SigDown[ind_dnFin_false] <- "Downtrend False Stop"
  #Create signal var over "window of opportunity" determined by params
  #before_window and aft_window. before_window should be greater than aft_window.
  #----------
  #Long Sig
  #----------
  #Hold
  df$SigUpWin <- "Hold"
  #Buy/Sell
  for(i in 1:n_upFin_keep){df$SigUpWin[(ind_upFin_keep[i] - before_window):(ind_upFin_keep[i] + aft_window)] <- "Uptrend Stop"}
  for(i in 1:n_upBeg_keep){df$SigUpWin[(ind_upBeg_keep[i] - before_window):(ind_upBeg_keep[i] + aft_window)] <- "Uptrend Start"}
  #Overwrite buy/sell sigs outside pct return thresholds as false signals
  for(i in 1:n_upFin_false){df$SigUpWin[(ind_upFin_false[i] - before_window):(ind_upFin_false[i] + aft_window)] <- "Uptrend False Stop"}
  for(i in 1:n_upBeg_false){df$SigUpWin[(ind_upBeg_false[i] - before_window):(ind_upBeg_false[i] + aft_window)] <- "Uptrend False Start"}
  #----------
  #Short Sig
  #----------
  #Hold
  df$SigDownWin <- "Hold"
  for(i in 1:n_dnFin_keep){df$SigDownWin[(ind_dnFin_keep[i] - before_window):(ind_dnFin_keep[i] + aft_window)] <- "Downtrend Start"}
  for(i in 1:n_dnBeg_keep){df$SigDownWin[(ind_dnBeg_keep[i] - before_window):(ind_dnBeg_keep[i] + aft_window)] <- "Downtrend Stop"}
  #Overwrite buy/sell sigs outside thresholds as false signals
  for(i in 1:n_dnFin_false){df$SigDownWin[(ind_dnFin_false[i] - before_window):(ind_dnFin_false[i] + aft_window)] <- "Downtrend False Start"}
  for(i in 1:n_dnBeg_false){df$SigDownWin[(ind_dnBeg_false[i] - before_window):(ind_dnBeg_false[i] + aft_window)] <- "Downtrend False Stop"}
  #-----------------------
  #Interesting... the up trend start (buy) points coincide with max ldydxcv points
  #And these are all at about the same value (for SPY anyway)
  ind_ind_muldydxcv_up <- which(df$ldydxcv[ind_upBeg_keep] > 0)
  mu_ldydx_cv_up <- mean(df$ldydxcv[ind_upBeg_keep[ind_ind_muldydxcv_up]])
  sd_ldydx_cv_up <- sd(df$ldydxcv[ind_upBeg[ind_ind_muldydxcv_up]])
  cv_ldydx_cv_up <- sd_ldydx_cv_up / mu_ldydx_cv_up
  #
  mu_dydx_mu_up <- mean(df$dydxmu[ind_upBeg_keep])
  sd_dydx_mu_up <- sd(df$dydxmu[ind_upBeg_keep])
  cv_dydx_mu_up <- sd_dydx_mu_up / mu_dydx_mu_up
  #
  mu_dydx_mu_dn <- mean(df$dydxmu[ind_dnBeg_keep])
  sd_dydx_mu_dn <- sd(df$dydxmu[ind_dnBeg_keep])
  cv_dydx_mu_dn <- sd_dydx_mu_dn / mu_dydx_mu_dn
  #---------------------------------
  #Get trimmed ts index to first uptrend start point and last downtrend stop point
  #so that ML training not messed up
  ind_tsStart <- min(ind_upBeg[1], ind_dnBeg[1]) - before_window
  ind_tsFinish <- max(ind_upFin[n_upFin], ind_dnFin[n_dnFin]) + aft_window
  # df <- df[ind_tsStart:ind_tsFinish, ]
  n_ts_trimd <- nrow(df)
  #---------------------------------
  outVars <- data.frame(
    ind_tsStart = ind_tsStart,
    ind_tsFinish = ind_tsFinish,
    
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
    cvLdydxCVup = cv_ldydx_cv_up,
    
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
  #Output df with the original ts and its rolling slope series (mean, sd, cv),
  #and signal vars.
  df_ts <- df[, c("Index", "ts", "dydxmu", "dydxsd", "ldydxcv", "SigUpWin", "SigDownWin")]
  #Output dfs with only trend start/finish (buy/sell) point data
  df_ts_upEvents <- subset(df, SigUp %in% c("Uptrend Start"))
  df_ts_upEvents <- df_ts_upEvents[, c("Index", "Ret_pct", "ts", "dydxmu", "dydxsd", "ldydxcv", "SigUp")]
  df_ts_dnEvents <- subset(df, SigDown %in% c("Downtrend Start"))
  df_ts_dnEvents <- df_ts_dnEvents[, c("Index", "Ret_pct", "ts", "dydxmu", "dydxsd", "ldydxcv", "SigDown")]
  #===================================================
  if(quietly == F){
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
        "Mean rolling mean slope: ", mu_dydx_mu, "\n",
        "CV roling mean slope: ", cv_dydx_mu, "\n",
        "Mean rolling mean slope at kept uptrend starts: ", mu_dydx_mu_up, "\n",
        "CV rolling mean slope at kept uptrend starts: ", cv_dydx_mu_up, "\n",
        "Mean rolling mean slope at kept downtrend starts: ", mu_dydx_mu_dn, "\n",
        "CV rolling mean slope at kept downtrend starts: ", cv_dydx_mu_dn, "\n",
        "Mean logged rolling cv of (abs val of) slope at kept uptrend starts: ", mu_ldydx_cv_up, "\n",
        "CV logged rolling cv of (abs val of) slope at kept uptrend starts: ", cv_ldydx_cv_up
    )
    #===================================================
    #4 Graphs
    #===================================================
    #===================================================
    #1
    df_plotUp <- df[, c("Index", "SigUp", "SigDown", "ts", "Ret_pct", "trendDuration", "dydxmu", "dydxsd", "ldydxcv")]
    df_plotUp$Ret_pct <- df_plotUp$Ret_pct * 100
    df_plotUp$SigUpStartDate <- df_plotUp$Index
    ind_UptrendStart <- which(df$SigUp == "Uptrend Start")
    df_plotUp$SigUpStartDate[-ind_UptrendStart] <- NA
    df_plotUp$SigUpStopDate <- df_plotUp$Index
    ind_UptrendStop <- which(df$SigUp == "Uptrend Stop")
    df_plotUp$SigUpStopDate[-ind_UptrendStop] <- NA
    df_plotUp$Ret_pct[-ind_UptrendStart] <- NA
    df_plotUp$trendDuration[-ind_UptrendStart] <- NA
    df_plotUp$`Pct Change/Time` <- df_plotUp$Ret_pct / df_plotUp$trendDuration
    #--
    #Plot of just the ts with trend start/finish (buy/sell) points overlaid
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotUp, aes(x = Index, y = ts))
    gg <- gg + geom_point(data = df_plotUp, aes(x = SigUpStartDate, y = ts), color = "green")
    gg <- gg + geom_point(data = df_plotUp, aes(x = SigUpStopDate, y = ts), color = "red")
    gg <- gg + geom_rect(aes(xmin = df_plotUp$Index[ind_UptrendStart], xmax = df_plotUp$Index[ind_UptrendStop],
                             ymin = -Inf, ymax = Inf, fill = "Uptrends"), alpha = 0.3)
    gg <- gg + scale_fill_manual("", values = "green", guide = guide_legend(override.aes = list(alpha = 1)))
    gg <- gg + theme(axis.title.x = element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
    gg <- gg + ylab("Time Series")
    #    gg <- gg + theme_light()
    gg_Up1 <- gg
    #--
    #Accompanying plot of false uptrends
    ind_UptrendStart_false <- which(df$SigUp == "Uptrend False Start")
    ind_UptrendStop_false <- which(df$SigUp == "Uptrend False Stop")
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotUp, aes(x = Index, y = ts))
    gg <- gg + geom_point(data = df_plotUp, aes(x = SigUpStartDate, y = ts), color = "green")
    gg <- gg + geom_point(data = df_plotUp, aes(x = SigUpStopDate, y = ts), color = "red")
    gg <- gg + geom_rect(aes(xmin = df_plotUp$Index[ind_UptrendStart_false], xmax = df_plotUp$Index[ind_UptrendStop_false],
                             ymin = -Inf, ymax = Inf, fill = "False Uptrends"), alpha = 0.3)
    gg <- gg + scale_fill_manual("", values = "violet", guide = guide_legend(override.aes = list(alpha = 1)))
    gg <- gg + ylab("Time Series")
    gg <- gg + theme(axis.title.x = element_blank())#,
                     # axis.text.x=element_blank(),
                     # axis.ticks.x=element_blank())
    #    gg <- gg + theme_light()
    gg_UpFalse1 <- gg
    #--
    #Accompanying plot of the pct change, shaded to reflect pct change/time also
    df_plotUp_ret <- data.frame(xmin = df_plotUp$SigUpStartDate[ind_UptrendStart],
                                xmax = df_plotUp$SigUpStopDate[ind_UptrendStop],
                                ymin = 0, ymax = df_plotUp$Ret_pct[ind_UptrendStart], 
                                pctch_per_time = df_plotUp$`Pct Change/Time`[ind_UptrendStart])
    colnames(df_plotUp_ret)[ncol(df_plotUp_ret)] <- "Pct Change/Time"
    df_plotUp$zero <- 0
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotUp, aes(x = Index, y = zero))
    gg <- gg + geom_rect(data = df_plotUp_ret,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = `Pct Change/Time`), alpha = 0.9)
    gg <- gg + scale_fill_gradient(low = "white", high = "green")
    gg <- gg + ylab("Pct Change")
    gg <- gg + theme(axis.title.x = element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
    #    gg <- gg + theme_light()
    gg_Up2 <- gg
    #--
    #Arrange these plots into one display
    gg1 <- ggarrange(gg_Up1, gg_Up2, gg_UpFalse1, ncol = 1, nrow = 3, align = "v")
    print(gg1)
    #===================================================
    #2
    #Do same for down trends
    df_plotDn <- df[, c("Index", "SigUp", "SigDown", "ts", "Ret_pct", "trendDuration", "dydxmu", "dydxsd", "ldydxcv")]
    df_plotDn$Ret_pct <- df_plotDn$Ret_pct * 100
    df_plotDn$SigDownStartDate <- df_plotDn$Index
    ind_DowntrendStart <- which(df$SigDown == "Downtrend Start")
    df_plotDn$SigDownStartDate[-ind_DowntrendStart] <- NA
    df_plotDn$SigDownStopDate <- df_plotDn$Index
    ind_DowntrendStop <- which(df$SigDown == "Downtrend Stop")
    df_plotDn$SigDownStopDate[-ind_DowntrendStop] <- NA
    df_plotDn$Ret_pct[-ind_DowntrendStart] <- NA
    df_plotDn$trendDuration[-ind_DowntrendStart] <- NA
    df_plotDn$`Pct Change/Time` <- df_plotDn$Ret_pct / df_plotDn$trendDuration
    #--
    #Plot of just the ts with trend start/finish (buy/sell) points overlaid
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotDn, aes(x = Index, y = ts))
    gg <- gg + geom_point(data = df_plotDn, aes(x = SigDownStartDate, y = ts), color = "green")
    gg <- gg + geom_point(data = df_plotDn, aes(x = SigDownStopDate, y = ts), color = "red")
    gg <- gg + geom_rect(aes(xmin = df_plotDn$Index[ind_DowntrendStart], xmax = df_plotDn$Index[ind_DowntrendStop],
                             ymin = -Inf, ymax = Inf, fill = "Downtrends"), alpha = 0.3)
    gg <- gg + theme(axis.title.x = element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
    gg <- gg + ylab("Time Series")
    gg <- gg + scale_fill_manual("", values = "green", guide = guide_legend(override.aes = list(alpha = 1)))
    #gg <- gg + theme_light()
    gg_Dn1 <- gg
    #--
    #Accompanying plot of false downtrends
    ind_DowntrendStart_false <- which(df$SigDown == "Downtrend False Start")
    ind_DowntrendStop_false <- which(df$SigDown == "Downtrend False Stop")
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotDn, aes(x = Index, y = ts))
    gg <- gg + geom_point(data = df_plotDn, aes(x = SigDownStartDate, y = ts), color = "green")
    gg <- gg + geom_point(data = df_plotDn, aes(x = SigDownStopDate, y = ts), color = "red")
    gg <- gg + geom_rect(aes(xmin = df_plotDn$Index[ind_DowntrendStart_false], xmax = df_plotDn$Index[ind_DowntrendStop_false],
                             ymin = -Inf, ymax = Inf, fill = "False Downtrends"), alpha = 0.3)
    gg <- gg + ylab("Time Series")
    gg <- gg + scale_fill_manual("", values = "violet", guide = guide_legend(override.aes = list(alpha = 1)))
    gg <- gg + theme(axis.title.x = element_blank())
    #gg <- gg + theme_light()
    gg_DnFalse1 <- gg
    #--
    #Accompanying plot of the pct change, shaded to reflect pct change/time also
    df_plotDn_ret <- data.frame(xmin = df_plotDn$SigDownStartDate[ind_DowntrendStart],
                                xmax = df_plotDn$SigDownStopDate[ind_DowntrendStop],
                                ymin = 0, ymax = df_plotDn$Ret_pct[ind_DowntrendStart], 
                                pctch_per_time = df_plotDn$`Pct Change/Time`[ind_DowntrendStart])
    colnames(df_plotDn_ret)[ncol(df_plotDn_ret)] <- "Pct Change/Time"
    df_plotDn$zero <- 0
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plotDn, aes(x = Index, y = zero))
    gg <- gg + geom_rect(data = df_plotDn_ret,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = `Pct Change/Time`), alpha = 0.9)
    gg <- gg + scale_fill_gradient(low = "green", high = "white")
    gg <- gg + ylab("Pct Change")
    gg <- gg + theme(axis.title.x = element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
    #gg <- gg + theme_light()
    gg_Dn2 <- gg
    #--
    #Arrange these plots into one display
    gg2 <- ggarrange(gg_Dn1, gg_Dn2, gg_DnFalse1, ncol = 1, nrow = 3, align = "v")
    print(gg2)
    #===================================================
    #3
    #Plot roling slope (mean, sd, and cv) with trend start/finish (buy/sell) points overlaid
    #The dydxmu plot includes the mean
    #The ldydxcv plot includes the mean only of buy points > 0
    df_plot <- df_plotUp %>% gather(ts_type, Value, dydxmu:ldydxcv)
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plot, aes(x = Index, y = zero), color = "red", linetype = "dotted")
    gg <- gg + geom_rect(data = df_plotUp_ret,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = `Pct Change/Time`), alpha = 0.9)
    gg <- gg + scale_fill_gradient(low = "white", high = "green")
    gg <- gg + geom_line(data = df_plot, aes(x = Index, y = Value))
    # gg <- gg + theme(#axis.title.x = element_blank(),
    #                  #axis.text.x = element_blank(),
    #                  #axis.ticks.x = element_blank(),
    #                  #axis.title.y = element_blank()
    #                  )
    gg <- gg + facet_wrap(~ ts_type, ncol = 1, scales = "free")
    gg <- gg + xlab("Date") + ylab("Value")
    gg <- gg + geom_hline(data = data.frame(yint = mu_dydx_mu_up, ts_type = "dydxmu"), aes(yintercept = yint), color = "violet", linetype = "dashed")
    gg <- gg + geom_hline(data = data.frame(yint = mu_ldydx_cv_up, ts_type = "ldydxcv"), aes(yintercept = yint), color = "violet",linetype = "dashed")
    gg3 <- gg
    print(gg3)
    #===================================================
    #4
    #Now same for downtrends
    df_plot <- df_plotDn %>% gather(ts_type, Value, dydxmu:ldydxcv)
    gg <- ggplot()
    gg <- gg + geom_line(data = df_plot, aes(x = Index, y = zero), color = "red", linetype = "dotted")
    gg <- gg + geom_rect(data = df_plotDn_ret,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = `Pct Change/Time`), alpha = 0.9)
    gg <- gg + scale_fill_gradient(low = "green", high = "white")
    gg <- gg + geom_line(data = df_plot, aes(x = Index, y = Value))
    # gg <- gg + theme(#axis.text.x = element_blank(),
    #                  #axis.ticks.x = element_blank(),
    #                  #axis.title.y = element_blank()
    #                 )
    gg <- gg + facet_wrap(~ ts_type, ncol = 1, scales = "free")
    gg <- gg + xlab("Date") + ylab("Value")
    gg <- gg + geom_hline(data = data.frame(yint = mu_dydx_mu_up, ts_type = "dydxmu"), aes(yintercept = yint), color = "violet", linetype = "dashed")
    gg <- gg + geom_hline(data = data.frame(yint = mu_ldydx_cv_up, ts_type = "ldydxcv"), aes(yintercept = yint), color = "violet",linetype = "dashed")
    gg4 <- gg
    print(gg4)
    #---------------------------------------------------
    gg_finalUp <- ggarrange(gg1, gg3, ncol = 1, nrow = 2)
    gg_finalDn <- ggarrange(gg2, gg4, ncol = 1, nrow = 2)
    print(gg_finalUp)
    print(gg_finalDn)
  }
  #=================================
  outlist <- list(df_ts, df_ts_upEvents, df_ts_dnEvents, outVars)
  return(outlist)
}
