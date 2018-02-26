tsSlope <- function(in_ts, wind = 17, graph = T)
{
  s_start <- wind + 1
  s_end <- length(in_ts) - wind
  x <- c(1:(2 * wind + 1))
  dydx_mu <- c()
  dydx_sd <- c()
  t <- 0
  for(s in s_start:s_end)
  {
    datseg <- in_ts[(s - wind):(s + wind)]
    q <- polyfit(x, datseg, n = 2)
    a2 <- q[1]; a1 <- q[2]; a0 <- q[3]
    this_slope <- 2 * a2 * x + a1
    t <- t + 1
    dydx_mu[t] <- mean(this_slope)
    dydx_sd[t] <- sd(this_slope)
  }
  dydx_cv <- dydx_sd / dydx_mu
  #-----------------------------
  ts_in <- fortify(in_ts[s_start:s_end, ])
  colnames(ts_in)[2] <- "ts"
  n_diff <- 4
  ts_in$tsDiff <- c(rep(NA, n_diff), diff(ts_in$ts, differences = n_diff))
  df_out <- data.frame(ts_in, dydx_mu, dydx_sd, dydx_cv)
  df_out$Index <- as.Date(df_out$Index)
  colnames(df_out)[1] <- "Date"
  #-----------------------------
  if(graph == T)
  { 
    ind_zeros <- which(abs(df_out$dydx_mu) < 10^-3)
    sd_diff <- sd(df_out$tsDiff, na.rm = T)
    ind_diffs <- which(abs(df_out$tsDiff) > 3 * sd_diff)
    df_plot <- df_out %>% gather(Series, Value, ts:dydx_sd)
    gg <- ggplot(df_plot, aes(x = Date, y = Value))
    gg <- gg + geom_line()
    gg <- gg + geom_vline(xintercept = df_out$Date[ind_diffs], color = "green")
    gg <- gg + facet_wrap(~Series, ncol = 1, scales = "free")
    gg
  }
  #-----------------------------
  return(df_out)
}
