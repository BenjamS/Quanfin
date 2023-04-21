#https://online.stat.psu.edu/stat510/lesson/8/8.2
library(tidyverse)
library(tidyquant)
library(patchwork)
library(lubridate)
#=============================================================================
av_api_key("HQBAWJK4Y3YW81VG")
#=============================================================================
# Define functions
get_S_and_corrXS <- function(mat_X_in){
  # mat_P = eigenvectors of the data correlation matrix
  # mat_G = corresponding eigenvalues
  mat_X_centered <- scale(mat_X_in, scale = F)
  # out_svd <- svd(mat_X_centered)
  # sing_values <- out_svd$d
  # n_obs <- nrow(mat_X_centered)
  # eig_values <- sing_values^2 / (n_obs - 1)
  # mat_P <- out_svd$v
  
  #mat_X_in[which(is.na(mat_X_in[, 15])), 15]
  
  
  mat_P <- eigen(cov(mat_X_in))$vectors
  if(mean(mat_P[, 1]) < 0){mat_P <- -mat_P}
  eig_values <- eigen(cov(mat_X_in))$values
  mat_G <- diag(eig_values)
  
  #mat_P_sigs <- mat_P[, 1:n_signals]
  # eig_values[1:n_signals] / eigen(cov(mat_X_centered))$values[1:n_signals] #check
  # mat_P / eigen(cov(mat_X_in))$vectors #check
  
  #mat_G <- diag(eig_values)
  
  #mat_G_sigs <- matU[, 1:n_signals]
  #---------------------------------------------
  sd_X <- apply(mat_X_in, 2, sd)
  D_sdX_inv <- diag(1 / sd_X)
  cormat_XS <- D_sdX_inv %*% mat_P %*% sqrt(mat_G)
  row.names(cormat_XS) <- colnames(mat_X_in)
  mat_L <- cormat_XS
  #mat_L <- diag(1 / apply(mat_X_in, 2, sd)) %*% mat_P %*% sqrt(mat_G)
  #---------------------------------------------------------
  # Set sign of eigenvectors such that signals best conform to their most highly correlated items
  # First have to get average of highest correlated items for each signal
  corrThresh <- 0.55
  n_items <- ncol(mat_L)
  list_X_hiCorr_avg <- list()
  for(i in 1:n_items){
    this_loadvec <- mat_L[, i]
    ind_tracks <- which(abs(this_loadvec) >= corrThresh)
    if(length(ind_tracks) == 0){
      ind_tracks <- which(abs(this_loadvec) == max(abs(this_loadvec)))
    }
    if(length(ind_tracks) == 1){
      list_X_hiCorr_avg[[i]] <- mat_X_centered[, ind_tracks]
    }else{
      loadvec_kept <- this_loadvec[ind_tracks]
      list_X_hiCorr_avg[[i]] <- rowMeans(mat_X_centered[, ind_tracks])
      
    }
  }
  mat_X_hiCorr_avg <- do.call(cbind, list_X_hiCorr_avg)
  mat_S_all <- mat_X_centered %*% mat_P
  #mat_S_all <- mat_X_in %*% mat_P
  for(i in 1:n_items){
    this_S <- mat_S_all[, i]
    this_X_hiCorr_avg <- mat_X_hiCorr_avg[, i]
    mse <- mean((this_S - this_X_hiCorr_avg)^2)
    mse_neg <- mean((-this_S - this_X_hiCorr_avg)^2)
    if(mse_neg < mse){
      mat_P[, i] <- -mat_P[, i]
    }
  }
  cormat_XS <- D_sdX_inv %*% mat_P %*% sqrt(mat_G)
  row.names(cormat_XS) <- colnames(mat_X_in)
  mat_L <- cormat_XS
  mat_S_all <- mat_X_centered %*% mat_P
  #---------------------------------------------
  # res <- FactoMineR::PCA(mat_pctDiff_in, scale.unit = F, ncp = ncol(mat_pctDiff_in), graph = F)
  # mat_L_FactoMiner <- res$var$coord
  # mat_L / mat_L_FactoMiner
  
  list_out <- list(mat_S_all, cormat_XS, eig_values, mat_P)
  return(list_out)
}
#-----------------------------------------------------------------------------
# Function to plot signal-item correlations (loadings)
plot_corrXS_barchart <- function(mat_L, group_info = NULL, xAxis_title = NULL, sigNames = NULL){
  
  n_signals <- ncol(mat_L)
  df_plot <- data.frame(Item = row.names(mat_L), mat_L)
  df_plot$Item <- as.character(df_plot$Item)
  #-------------------------------------------------------
  if(is.null(sigNames)){
    signal_id <- paste("Signal", 1:n_signals)
  }else{
    #signal_id <- paste("Signal", 1:n_signals, "\n", sigNames)
    signal_id <- sigNames
  }
  colnames(df_plot)[2:(n_signals + 1)] <- signal_id
  #-------------------------------------------------------
  gathercols <- as.character(signal_id) 
  df_plot <- gather_(df_plot, "Signal", "Correlation", gathercols)
  df_plot <- transform(df_plot,
                       Signal = factor(Signal, levels = gathercols))
  
  if(!is.null(group_info)){
    outlist <- group_fn(group_info)
    cols_ordered_by_group <- outlist[[1]]
    group_color_vec <- outlist[[2]]
    group_vec_ordered <- outlist[[3]]
    df_match_group <- data.frame(Item = cols_ordered_by_group, Group = group_vec_ordered)
    df_plot <- merge(df_plot, df_match_group, by = "Item")
    df_plot <- df_plot[order(df_plot$Group), ]
    df_plot$Item <- factor(df_plot$Item, levels = unique(df_plot$Item))
    gg <- ggplot(df_plot, aes(x = Item, y = Correlation, fill = Group))
    gg <- gg + scale_fill_manual(values = unique(group_color_vec))
  }else{
    df_plot$Item <- factor(df_plot$Item,
                           levels = rev(unique(df_plot$Item)))
    gg <- ggplot(df_plot, aes(x = Item, y = Correlation))
  }
  gg <- gg + geom_bar(stat = "identity", color = "black", position = "dodge")
  gg <- gg + ylim(limits = c(-1, 1))
  gg <- gg + facet_wrap(~ Signal, nrow = 1)
  if(!is.null(xAxis_title)){
    gg <- gg + labs(y = xAxis_title)
  }
  gg <- gg + theme(axis.text = element_text(size = 7),
                   axis.title.x = element_text(size = 7),
                   axis.title.y = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 7),
                   strip.text = element_text(size = 7))
  gg <- gg + coord_equal()
  gg <- gg + coord_flip()
  gg
  
}
#-----------------------------------------------------------------------------
# For percentile oscillator
pctileFun <- function(x){
  out <- ecdf(x)(x[length(x)])
  return(out)
}
#----------------------------------------------------------------------------
# Extract seasonality
# getSeasons <- function(df_in,
#                        freq = 52,
#                        mod_type = "multiplicative",
#                        show_graphs = T){
#   # daily_freq <- 104
#   # weekly_freq <- 52
#   # monthly_freq <- 12
#   # freq <- daily_freq
#   tsStart <- c(year(df_in$date)[1], week(df_in$date)[1])
#   #tsStart <- c(year(df_in$date)[1], month(df_in$date)[1])
#   item_vec <- unique(df_in$Item)
#   n_items <- length(item_vec)
#   ratio_vec <- c()
#   ratioCV_vec <- c()
#   list_df <- list()
#   for(i in 1:n_items){
#     this_item <- item_vec[i]
#     print(this_item)
#     this_dateVec <- subset(df_in, Item == this_item)$date
#     this_series <- subset(df_in, Item == this_item)$Value
#     names(this_series) <- this_dateVec
#     this_series <- na.approx(this_series)
#     this_series <- na.trim(this_series)
#     this_dateVec <- names(this_series)
#     this_ts <- ts(this_series, start = tsStart, frequency = freq)
#     #plot.ts(this_ts)
#     ts_decomp <- this_ts %>%
#       decompose(type = mod_type) #"additive" or "multiplicative"
#     df_out <- data.frame(date = this_dateVec, Item = this_item, Value = as.numeric(ts_decomp$seasonal))
#     list_df[[i]] <- df_out
#     mu_ratio <- mean(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T)
#     ratio_vec[i] <- mu_ratio
#     ratioCV_vec[i] <- sd(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T) / mu_ratio
#     if(show_graphs){
#       plot(ts_decomp)
#       Sys.sleep(2)
#     }
#     
#   }
#   # names(ratio_vec) <- item_vec
#   # names(ratioCV_vec) <- item_vec
#   # hist(ratio_vec)
#   # hist(ratioCV_vec)
#   df_s <- as.data.frame(do.call(rbind, list_df))
#   list_out <- list(df_s, ratio_vec, ratioCV_vec)
#   return(list_out)
# }
#=============================================================================
#=============================================================================
# End function definition
#=============================================================================
#=============================================================================
#=============================================================================
fx_vec1 <- c("eur/gbp", "eur/usd",
            "aud/usd", "usd/jpy",
            "aud/jpy")
fx_vec2 <- c("eur/jpy", "usd/cad", "usd/chf")
stock_vec <- c("ES=F", "GC=F", "NG=F", "CL=F", "CT=F", "KC=F", "CC=F",
               "SB=F", "ZB=F", "ZN=F", "ZF=F", "XLF", "XLI", "XLK",
               "IYR", "XLV", "XLY", "XLP")
fromdate <- "2011-01-01"
#-----------------------------------------------------------------------------
# fx_emaPer_vec <- rep(34, length(fx_vec))
# stock_emaPer_vec <- rep(55, length(stock_vec))
# emaPer_vec <- c(fx_emaPer_vec, stock_emaPer_vec)
#-----------------------------------------------------------------------------
time_step <- "weekly" #1min, 5min, 15min, 30min, 60min, daily, weekly
#=============================================================================
df_ohlcv <- stock_vec %>% 
  tq_get(get = "stock.prices",
         from = fromdate, periodicity = time_step) %>%
  as.data.frame()
df_ohlcv$p <- df_ohlcv$adjusted
df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
dfStk <- df_ohlcv[, c("symbol", "date", "p", "volume", "diffHiLo")]
# AlphaVantage allows only 5 queries per minute.
# So split queries into two batches with minute wait in between.
df_ohlcv1 <- fx_vec1 %>%  tq_get(get = "alphavantager",
                               av_fun = "FX_WEEKLY",
                               outputsize = "full") %>% as.data.frame()
Sys.sleep(61)
df_ohlcv2 <- fx_vec2 %>%  tq_get(get = "alphavantager",
                               av_fun = "FX_WEEKLY",
                               outputsize = "full") %>% as.data.frame()
df_ohlcv <- as.data.frame(rbind(df_ohlcv1, df_ohlcv2))
df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
colnames(df_ohlcv)[2] <- "date"
dfFx <- df_ohlcv[, c("symbol", "date", "p", "diffHiLo")]
#-----------------------------------------------------------------------------
o <- apply(dfFx, 2, function(x) length(which(is.na(x))))
#table(o)
ind_na <- which(is.na(df$p))
if(length(ind_na) != 0){
  df$p <- na.approx(df$p)
}
#-----------------------------------------------------------------------------
# Get percentile oscillator series
rollWind <- 34
dfStk <- dfStk %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
dfFx <- dfFx %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
#-----------------------------------------------------------------------------
# Get p - ema oscillator
dfStk <- dfStk %>% group_by(symbol) %>%
  mutate(ema = EMA(p, per_ema_for_detrend)) %>%
  mutate(difEMA = p - ema) %>%
  as.data.frame()
dfFx <- dfFx %>% group_by(symbol) %>%
  mutate(ema = EMA(p, per_ema_for_detrend)) %>%
  mutate(difEMA = p - ema) %>%
  as.data.frame()
#-----------------------------------------------------------------------------
gg <- ggplot(dfFx, aes(x = date, y = pctlOsc))
gg <- gg + geom_line()
gg <- gg + scale_x_date(breaks = scales::breaks_pretty(n = 4), labels = scales::date_format("%b\n%Y"))
gg <- gg + facet_wrap(~symbol, scales = "free_y")
gg
#-----------------------------------------------------------------------------
#df_pca <- dfFx[which(year(dfFx$date) > 2018), c("symbol", "date", "pctlOsc")]
df <- as.data.frame(rbind(dfFx[, c("symbol", "date", "pctlOsc")], dfStk[, c("symbol", "date", "pctlOsc")]))
df_pca <- df[which(year(df$date) > 2018), ]
df_pca <- df_pca %>% spread(symbol, pctlOsc)
# Remove leading NAs if using percentile oscillator
ind_rm <- 1:(rollWind - 1)
df_pca <- df_pca[-ind_rm, ]
o <- apply(df_pca, 2, function(x) length(which(is.na(x))))
#table(o)
#ind_na <- which(is.na(df$p))
ind_na <- which(o > 0)
if(length(ind_na) != 0){
  df_pca[, ind_na] <- na.approx(df_pca[, ind_na])
}
# Shorten names
# colnames(df_pca) <- gsub("ICE FUTURES U.S.", "ICE", colnames(df_pca))
# colnames(df_pca) <- gsub("CHICAGO BOARD OF TRADE", "CBoT", colnames(df_pca))
# colnames(df_pca) <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", colnames(df_pca))
# colnames(df_pca) <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", colnames(df_pca))
# colnames(df_pca) <- gsub("COMMODITY EXCHANGE INC.", "CE", colnames(df_pca))
mat_X_in <- as.matrix(df_pca[, -1])
o <- apply(mat_X_in, 1, function(x) length(which(is.na(x))))
#table(o)
row.names(mat_X_in) <- df_pca$date
mat_X_in <- mat_X_in[-which(o > 0), ]
out <- get_S_and_corrXS(mat_X_in)
mat_L <- out[[2]]
eigVals <- out[[3]]
pctExplnd <- cumsum(eigVals) / sum(eigVals)
ind90pctExplnd <- which(pctExplnd >= 0.90)[1]
pctExplndVec <- round(eigVals[1:ind90pctExplnd], 3)
mat_L <- mat_L[, 1:ind90pctExplnd]
plot_corrXS_barchart(mat_L, group_info = NULL, xAxis_title = NULL,
                     sigNames = as.character(pctExplndVec))
mat_L <- out[[2]]
mat_Lrot <- varimax(mat_L)[[1]]
mat_Lrot <- matrix(as.numeric(mat_Lrot),
                   attributes(mat_Lrot)$dim,
                   dimnames = attributes(mat_Lrot)$dimnames)
mat_R <- varimax(mat_L)[[2]]
mat_R <- matrix(as.numeric(mat_R),
                attributes(mat_R)$dim,
                dimnames = attributes(mat_R)$dimnames)
xAxis_title <- "Varimax Rotated Correlation"
mat_Lrot <- mat_Lrot[, 1:5]
plot_corrXS_barchart(mat_Lrot, group_info = NULL, xAxis_title, sigNames = NULL)





























# Read in the Commitment of Traders (CoT) data
# Financials CoT data
this_folder <- "D:/OneDrive - CGIAR/Documents 1/CIAT 2/finAnalysis/data/"
this_filename <- "Finance CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_finCoT <- read.csv(this_filepath, stringsAsFactors = F)
# Commodities CoT data
this_filename <- "Commodity CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_comCoT <- read.csv(this_filepath, stringsAsFactors = F)
#-----------------------------------------------------------------------------
# Get rolling percentile series
# Set rolling percentile window size for weekly CoT data
rollWind_CoT <- 52
#---
df_pctlFinCoT <- subset(df_finCoT, Element == "Smart money net position (% of OI)")
df_pctlFinCoT$Element <- NULL
df_pctlFinCoT <- df_pctlFinCoT %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind_CoT, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
df_pctlComCoT <- subset(df_comCoT, Element == "Smart money net position (% of OI)")
df_pctlComCoT$Element <- NULL
df_pctlComCoT <- df_pctlComCoT %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind_CoT, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
#-----------------------------------------------------------------------------
df_in <- df_pctlComCoT
df_in <- df_in[-which(is.na(df_in$Value)), ]
out_list <- getSeasons(df_in,
                       freq = 52,
                       mod_type = "additive",
                       show_graphs = T)
df_s <- out_list[[1]]
ratio1 <- out_list[[2]]
ratio2 <- out_list[[3]]