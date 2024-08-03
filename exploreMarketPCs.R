#https://online.stat.psu.edu/stat510/lesson/8/8.2
library(tidyverse)
library(tidyquant)
library(patchwork)
library(lubridate)
library(scales)
library(WaveletComp)
library(kableExtra)
library(xtable)
library(flextable)
library(randomcoloR)
#=============================================================================
# Define functions
get_S_and_corrXS <- function(mat_X_in){
  # mat_P = eigenvectors of the data correlation matrix
  # mat_G = corresponding eigenvalues
  #mat_X_centered <- scale(mat_X_in, scale = F)
  # out_svd <- svd(mat_X_centered)
  # sing_values <- out_svd$d
  # n_obs <- nrow(mat_X_centered)
  # eig_values <- sing_values^2 / (n_obs - 1)
  # mat_P <- out_svd$v
  
  #mat_X_in[which(is.na(mat_X_in[, 15])), 15]
  
  mat_P <- eigen(cov(mat_X_in))$vectors
  if(mean(mat_P[, 1]) < 0){mat_P <- -mat_P}
  eig_values <- eigen(cov(mat_X_in))$values
  eig_values[which(eig_values < 0)] <- 0
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
      list_X_hiCorr_avg[[i]] <- mat_X_in[, ind_tracks]
    }else{
      loadvec_kept <- this_loadvec[ind_tracks]
      list_X_hiCorr_avg[[i]] <- rowMeans(mat_X_in[, ind_tracks])
      
    }
  }
  mat_X_hiCorr_avg <- do.call(cbind, list_X_hiCorr_avg)
  mat_S_all <- mat_X_in %*% mat_P %*% diag(1 / sqrt(eig_values)) * 1 / 2
  #mat_S_all <- (mat_S_all + 1) * 1 / 2
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
  # Recalculate signals based on new P with changed signs
  mat_S_all <- mat_X_in %*% mat_P %*% diag(1 / sqrt(eig_values)) * 1 / 2
  #mat_S_all <- (mat_S_all + 1) * 1 / 2
  cormat_XS <- D_sdX_inv %*% mat_P %*% sqrt(mat_G)
  row.names(cormat_XS) <- colnames(mat_X_in)
  mat_L <- cormat_XS
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
  gg <- gg + geom_hline(yintercept = 0, color = "black")
  gg <- gg + geom_hline(yintercept = c(0.5, -0.5), color = "red")
  gg <- gg + geom_bar(stat = "identity", position = "dodge")
  gg <- gg + ylim(limits = c(-1, 1))
  gg <- gg + facet_wrap(~Signal, ncol = 1, strip.position = "right")
  # if(!is.null(xAxis_title)){
  #   gg <- gg + labs(x = xAxis_title)
  # }
  gg <- gg + theme_classic()
  gg <- gg + theme(axis.text.y = element_text(size = 7),
                   axis.text.x = element_text(size = 7, angle = 60, hjust = 1),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 7),
                   legend.key.size = unit(0.3, "cm"),
                   legend.position = "top",
                   panel.background = element_rect(fill = "black"),
                   strip.text = element_text(size = 7))
  gg <- gg + guides(fill = guide_legend(nrow = 2, byrow = T))
  #gg <- gg + coord_equal()
  #gg <- gg + coord_flip()
  gg
  
}
#--------------------------------------------------------------
# Define function to order data by group
group_fn <- function(groupInfo){
  listGroups <- groupInfo[[1]]
  groupNames <- groupInfo[[2]]
  groupColors <- groupInfo[[3]]
  varNames_ordered <- do.call(c, listGroups)
  nGroups <- length(groupNames)
  n_items <- length(varNames_ordered)
  if(is.null(groupColors)){
    bag_of_colors <- distinctColorPalette(k = 5 * nGroups)
    groupColors <- sample(bag_of_colors, nGroups)
    #group_colors <- viridis::viridis_pal(option = "D")(length(group_names))
  }
  #if(reverse_order){group_colors <- rev(group_colors)}
  #varNames_ordered <- colnames(mat_pctDiff)
  group_vec <- rep(NA, n_items)
  group_color_vec <- rep(NA, n_items)
  for(i in 1:nGroups){
    this_group_vec <- listGroups[[i]]
    this_group_name <- groupNames[i]
    this_group_color <- groupColors[i]
    group_vec[which(varNames_ordered %in% this_group_vec)] <- this_group_name
    group_color_vec[which(varNames_ordered %in% this_group_vec)] <- this_group_color
  }
  ind_ordered_cols <- order(factor(group_vec))
  cols_ordered_by_group <- as.character(varNames_ordered[ind_ordered_cols])
  group_color_vec <- group_color_vec[ind_ordered_cols]
  group_vec_ordered <- group_vec[ind_ordered_cols]
  out_list <- list(cols_ordered_by_group, group_color_vec, group_vec_ordered, ind_ordered_cols, group_vec)
  return(out_list)
}
#-----------------------------------------------------------------------------
# For percentile oscillator
pctileFun <- function(x){
  out <- ecdf(x)(x[length(x)])
  return(out)
}
#========================================================================
fitWave <- function(ts, per_vec, pval_thresh = 0.01, n_lookAhead){
  n_t <- length(ts)
  t <- 1:n_t
  regrsrs_sin <- paste0("sin(2 * pi / per_vec[", c(1:length(per_vec)), "] * t)", collapse = " + ")
  regrsrs_cos <- paste0("cos(2 * pi / per_vec[", c(1:length(per_vec)), "] * t)", collapse = " + ")
  regrsrs <- paste(regrsrs_sin, regrsrs_cos, sep = " + ")
  this_formula <- as.formula(paste("ts ~ -1 +", regrsrs)[1])
  linmod <- lm(this_formula)
  summod <- summary(linmod)
  #summod
  pvals <- as.numeric(summod$coefficients[, 4])
  ind_rm <- which(pvals > pval_thresh)
  round <- 0
  while(length(ind_rm) > 0){
    round <- round + 1
    print(paste("round ", round))
    if(round == 7){
      print("Too many rounds, aborting.")
      break
    }
    regrsrs_char <- strsplit(as.character(regrsrs), " \\+ ")[[1]]
    regrsrs_char <- regrsrs_char[-ind_rm]
    regrsrs <- paste(regrsrs_char, collapse = " + ")
    this_formula <- as.formula(paste("ts ~ -1 +", regrsrs)[1])
    linmod <- lm(this_formula)
    #---------------------
    summod <- summary(linmod)
    #print(summod)
    pvals <- as.numeric(summod$coefficients[, 4])
    ind_rm <- which(pvals > pval_thresh)
  }
  fitted_wave <- fitted(linmod)
  t <- (n_t + 1):(n_t + n_lookAhead)
  regrsrs_fut <- names(linmod$coefficients)
  #regrsrs_fut <- gsub("t", "t_fut", regrsrs_fut)
  list_fut <- list()
  for(i in 1:length(regrsrs_fut)){
    list_fut[[i]] <- eval(parse(text = regrsrs_fut[i]))
  }
  df_fut <- as.data.frame(do.call(cbind, list_fut))
  predict_wave <- predict(linmod, df_fut, interval = "prediction")  
  print(summod)
  
  list_out <- list(fitted_wave, predict_wave, summod)
  return(list_out)
}
#========================================================================
plot_validation <- function(dfYhat, dfYpred, dfPlotTs, colorVec = NULL){
  symbVec <- unique(dfPlotTs$symbol)
  nTs <- ncol(dfYhat)
  if(is.null(colorVec)){
    bag_of_colors <- distinctColorPalette(k = 5 * nTs)
    theseColors <- sample(bag_of_colors, nTs)
  }else{
    theseColors <- colorVec
  }
  ind_divide <- nrow(dfYhat)
  dfPlotMod <- as.data.frame(rbind(dfYhat, dfYpred))
  dateVec <- subset(dfPlotTs, symbol == symbVec[1])$date; dateChrVec <- as.character(dateVec); lenTs <- length(dateVec)
  dfPlotTs$date_chr <- as.character(dfPlotTs$date)
  dfPlotMod$date <- dateVec; dfPlotMod$date_chr <- dateChrVec
  my_breaks <- dateChrVec[seq.int(1, lenTs, length.out = 30)]
  dfPlotMod <- dfPlotMod %>% gather_("symbol", "tsMod", symbVec)
  #dfPlotMod$date_chr <- as.factor(dfPlotMod$date_chr)
  thisTitle <- "Model validation"
  gg <- ggplot()
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", lwd = 1)
  gg <- gg + geom_vline(xintercept = ind_divide, lwd = 1, color = "blue")
  gg <- gg + geom_line(data = dfPlotTs, aes(x = date_chr, y = ts, color = symbol, group = symbol))#, color = colors_dtFit[2], lwd = 1.1)
  gg <- gg + geom_line(data = dfPlotMod, aes(x = date_chr, y = tsMod, color = symbol, group = symbol), lwd = 1.1)#, color = colors_dtFit[2], lwd = 1.1)
  gg <- gg + scale_color_manual(values = theseColors)
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + facet_wrap(~symbol, ncol = 1, strip.position = "right")
  gg <- gg + labs(title = thisTitle)
  gg <- gg + theme_dark()
  gg <- gg + theme(axis.title = element_blank(),
                   #legend.title = element_blank(),
                   legend.position = "none",
                   axis.text.x = element_text(angle = 60, hjust = 1),
                   plot.title = element_text(size = 10))
  print(gg)
  # gg <- ggplot()
  # gg <- gg + geom_line(data = dfPlotMod, aes(x = date_chr, y = tsMod, color = Type, group = Type), lwd = 1.1)#, color = colors_dtFit[2], lwd = 1.1)
  # gg <- gg + geom_line(data = dfPlotTs, aes(x = date_chr, y = ts, color = Type, group = Type))#, color = colors_dtFit[1])
  # gg <- gg + scale_color_manual(values = colors_dtFit)
  # gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", lwd = 1)
  # gg <- gg + geom_vline(aes(xintercept = ind_divide), lwd = 1, color = "violet")
  # gg <- gg + theme_bw()
  # gg <- gg + scale_x_discrete(breaks = my_breaks)
  # gg <- gg + theme(axis.title = element_blank(),
  #                  legend.title = element_blank(),
  #                  axis.text.x = element_text(angle = 60, hjust = 1))
  # gg_dtFit <- gg
  
  # dfPlot_ts$date_chr <- as.factor(dfPlot_ts$date_chr)
  # dfPlot_yhat_p <- subset(dfPlot_ts, Type == "yhat_p")
  # dfPlot_p_ema <- subset(dfPlot_ts, Type != "yhat_p")
  # gg <- ggplot()
  # gg <- gg + geom_line(data = dfPlot_p_ema, aes(x = date_chr, y = Value, group = Type, color = Type))
  # gg <- gg + scale_color_manual(values = colors_ts[1:2])
  # gg <- gg + geom_line(data = dfPlot_yhat_p, aes(x = date_chr, y = Value, group = 1), color = colors_ts[3], lwd = 1.1)
  # gg <- gg + geom_vline(aes(xintercept = ind_divide), color = "violet", lwd = 1)
  # gg <- gg + theme_bw()
  # gg <- gg + scale_x_discrete(breaks = my_breaks)
  # gg <- gg + theme(axis.title = element_blank(),
  #                  axis.text.x = element_text(angle = 60, hjust = 1),
  #                  #legend.position = "bottom",
  #                  legend.title = element_blank())
  # gg_ts <- gg
  # 
  # gg_together <- gg_dtFit + gg_ts + plot_layout(ncol = 1, heights = c(2, 1))
  
  # print(gg_together)
}
#========================================================================
plot_prediction <- function(dfYhat, dfYpred, dfPlotTs, time_step, n_lookAhead = 34, n_lookAhead_zoom = 21, n_lookBack_zoom = 21,
                            colorVec = NULL, showCloseUp = F){
  symbVec <- unique(dfPlotTs$symbol)
  nTs <- ncol(dfYhat)
  if(is.null(colorVec)){
    bag_of_colors <- distinctColorPalette(k = 5 * nTs)
    theseColors <- sample(bag_of_colors, nTs)
  }else{
    theseColors <- colorVec
  }
  ind_end <- nrow(dfYhat)
  dateVec <- subset(dfPlotTs, symbol == symbVec[1])$date; dateChrVec <- as.character(dateVec); lenTs <- length(dateVec)
  dfPlotMod <- dfYhat
  dfPlotMod$date_chr <- dateChrVec
  dfPlotMod$set <- "fit"
  #dfPlot$t <- NULL
  #---------------------------------------
  time_step_unit <- as.character(stringr::str_extract_all(time_step, "[a-z]+")[[1]])
  if(time_step_unit == "min"){
    time_step_num <- as.numeric(stringr::str_extract_all(time_step, "[0-9]+")[[1]])
    #time_step_num <- 60 * time_step_num
  }
  
  if(time_step_unit == "daily"){
    time_step_num <- 1
  }
  
  if(time_step_unit == "weekly"){
    time_step_num <- 1
  }
  # #---------------------------------------
  date_fut <- 1:n_lookAhead * time_step_num
  if(time_step_unit == "min"){
    date_fut <- round(date_fut / 60, 2)
    time_step_unit <- "hrs"
  }
  if(time_step_unit == "weekly"){
    # date_fut <- round(date_fut * 7, 2)
    # time_step_unit <- "days"
    time_step_unit <- "weeks"
  }
  if(time_step_unit == "daily"){
    time_step_unit <- "days"
  }
  date_fut_chr <- as.character(paste("+", date_fut, time_step_unit))
  #----------------------------------------------------------------------
  #  df_add <- data.frame(date = date_fut, date_chr = date_fut_chr, slope = NA, yhat = ypredict[, 1], set = "predict")
  #df_add <- data.frame(date_chr = date_fut_chr, slope = NA, yhat = ypredict[, 1], set = "predict")
  dfAdd <- dfYpred
  dfAdd$date_chr <- date_fut_chr; dfAdd$set <- "predict"
  dfPlotMod <- as.data.frame(rbind(dfPlotMod, dfAdd))
  dateChrVec <- dfPlotMod$date_chr; lenTs <- length(dateChrVec)
  my_breaks <- dateChrVec[seq.int(1, lenTs, length.out = 30)]
  dfPlotMod$date_chr <- factor(dfPlotMod$date_chr, levels = dateChrVec)
  dfPlotMod <- dfPlotMod %>% gather_("symbol", "tsMod", symbVec)
  dfPlotTs$date_chr <- as.character(dfPlotTs$date)
  # this_title <- paste(c(paste(time_step, "chart")),
  #                       #paste(per_ema_for_detrend, "step detrend")),
  #                     collapse = ", ")
  thisTitle <- "Model prediction"
  gg <- ggplot()
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", size = 1)
  gg <- gg + geom_vline(xintercept = ind_end, color = "blue", size = 1)
  gg <- gg + geom_line(data = dfPlotMod, aes(x = date_chr, y = tsMod, group = symbol, color = symbol), lwd = 1.1)
  gg <- gg + geom_line(data = dfPlotTs, aes(x = date_chr, y = ts, group = symbol, color = symbol))
  gg <- gg + scale_color_manual(values = theseColors)
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + facet_wrap(~symbol, ncol = 1, strip.position = "right")
  gg <- gg + labs(title = thisTitle)
  gg <- gg + theme_dark()
  gg <- gg + theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1),
                   #legend.title = element_blank(),
                   legend.position = "none",
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(size = 10))
  #gg <- gg + scale_color_brewer(palette = "Dark2")
  print(gg)
  
  if(showCloseUp){
    # Zoom in
    ind_end_new <- n_lookBack_zoom
    dateChrZoomVec <- dateChrVec[(ind_end - n_lookBack_zoom):(ind_end + n_lookAhead_zoom)]
    lenTsZoom <- length(dateChrZoomVec)
    dfPlotModZoom <- dfPlotMod %>% mutate(date_chr = as.character(date_chr)) %>%
      subset(date_chr %in% dateChrZoomVec)
    dfPlotTsZoom <- dfPlotTs %>% subset(date_chr >= dateChrZoomVec[1])
    my_breaks_zoom <- dateChrZoomVec[seq.int(1, lenTsZoom, length.out = 30)]
    dfPlotModZoom <- dfPlotModZoom %>% spread(symbol, tsMod) %>%
      mutate(date_chr = factor(date_chr, levels = dateChrZoomVec)) %>%
      gather_("symbol", "tsMod", symbVec)
    gg <- ggplot()
    gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", size = 1)
    gg <- gg + geom_vline(xintercept = ind_end_new, color = "blue", size = 1)
    gg <- gg + geom_line(data = dfPlotModZoom, aes(x = date_chr, y = tsMod, group = symbol, color = symbol))#, color = colors_dtFit[1], lwd = 1.1)
    gg <- gg + geom_line(data = dfPlotTsZoom, aes(x = date_chr, y = ts, group = symbol, color = symbol))#, color = colors_dtFit[2])
    gg <- gg + scale_color_manual(values = theseColors)
    gg <- gg + scale_x_discrete(breaks = my_breaks_zoom)
    gg <- gg + facet_wrap(~symbol, ncol = 1, strip.position = "right")
    gg <- gg + labs(title = paste(thisTitle, "close up"))
    gg <- gg + theme_bw()
    gg <- gg + theme(axis.title = element_blank(),
                     axis.text.x = element_text(angle = 60, hjust = 1),
                     legend.position = "none",
                     axis.ticks.x = element_blank(),
                     plot.title = element_text(size = 10))
    
    print(gg)
  }
  
  
}
#========================================================================
getCycles <- function(waveAnalysis, plotPeriodogram = T){
  df_periodogram <- data.frame(period = waveAnalysis$Period, power = waveAnalysis$Power.avg)
  ind_critPoints <- findPeaks(df_periodogram$power)
  critPers <- df_periodogram$period[ind_critPoints]
  critPers
  #-------------------
  u <- df_periodogram$period
  if(max(u) > 1.5 * max(critPers)){
    ind_rm <- which(u > 1.5 * max(critPers))
    df_plot <- df_periodogram[-ind_rm, ]
  }else{
    df_plot <- df_periodogram
  }
  if(plotPeriodogram){
    gg <- ggplot(df_plot, aes(x = period, y = power))
    gg <- gg + geom_line()
    gg <- gg + geom_vline(xintercept = critPers, color = "cyan", size = 1.2)
    gg <- gg + theme_bw()
    print(gg)
  }
  #-------------------
  # Get periods, ordered from most power to least
  critPwrs <- df_periodogram$power[ind_critPoints]
  ind_order <- order(critPwrs, decreasing = T)
  per_vec <- critPers[ind_order]
  pwr_vec <- critPwrs[ind_order]
  dfMainCycles <- data.frame(Num = c(1:length(per_vec)), Period = per_vec, Power = pwr_vec)
  #-------------------
  return(list(dfMainCycles, df_plot))
  
}
#=============================================================================
#=============================================================================
#=============================================================================
# End function definition
#=============================================================================
#=============================================================================
dfDrivers <- NULL
#=============================================================================
#=============================================================================
# Manual global overview
spy_sector_symbs <- c("XLF", "XLC", "XLY", "XLP", "XLV", "XLK", "RWR",
                      "XLU", "XLI", "XBI", "IYT", "XLE", "SPY", "IJR", "IJH", "URTH") #"TTEK"
spy_sector_detail <- c("Financials", "Communications", "Luxury goods", "Consumer\ngoods",
                       "Healthcare", "Technology", "Real estate", "Utilities", "Industrial",
                       "Biotechnology", "Transportation", "Energy", "S&P 500 ETF",
                       "Small Cap ETF", "Mid Cap ETF", "Devd Mkts ETF") #"Gov. foreign aid"
# semicondctr_symbs <- c("TSM", "NVDA", "AVGO", "ASML", "AMD", "QCOM", "TXN",
#                        "AMAT", "ADI", "MU", "LRCX", "KLAC", "INTC")
# semicondctr_detail <- c("Taiwan Semiconductor Manufacturing Co. Ltd.",
#                         "Nvidia", "Broadcom", "ASML (Netherlands semicndctr firm)",
#                         "AMD", "Qualcomm", "Texas Instruments", "Applied Materials",
#                         "Analog Devices", "Micron Technology",
#                         "Lam Research", "KLA", "Intel")
minerals_symbs <- c("GLD", "IAU", "SLV", "PPLT", "CPER", "DBB") #"XME"
minerals_detail <- c("Gold a", "Gold b", "Silver", "Platinum", "Copper", "Industrial metals") #"US metals and mining"
agriculture_symbs <- c("WEAT", "CORN", "KROP", "SOYB", "CANE", "DBA", "TAGS")
agriculture_detail <- c("Wheat", "Maize", "AgTech", "Soybean", "Sugar", "General Ag a", "General Ag b")
energy_symbs <- c("WTI", "BNO", "USO", "WOOD", "ICLN", "UNG")
energy_detail <- c("Oil (W&T)", "Oil (Brent)", "Oil (US ETF)", "Timber", "Clean energy", "US natural gas")
currency_symbs <- c("EMLC", "UUP", "FXE", "FXY", "FXF", "FXC", "FXB", "FXA")
# currency_detail <- c("Emerging mkt currencies", "USD", "EUR", "JPY", "CHF", "CND", "GBP", "AUD")
currency_detail <- c("EUR/USD", "USD/JPY",
                     "USD/CHF", "USD/CAD", "GBP/USD", "AUD/USD", "USD/INR")
currency_symbs <- c("EURUSD=X", "JPY=X", "CHF=X", "CAD=X",
                    "GBPUSD=X", "AUDUSD=X", "INR=X")
emerg_mkt_symbs <- c("ELD", "BKF", "VWOB", "FXI")
emerg_mkt_detail <- c("Emerg mkts debt", "BRIC countries", "Emerg mkts gov. bonds", "China ETF")
crypto_symbs <- c("BLOK", "LEGR", "BITQ", "BTC=F")
crypto_detail <- c("Blockchain tech.", "Blockchain companies",
                   "Crypto Industry Innovators", "Bitcoin futures")
Tbond_symbs <- c("IEI", "IEF", "BIL")#"TLT")
Tbond_detail <- c("T-bond 3-7 yrs", "T-bond 7-10 yrs", "T-bond 1-3 months")#"T-bond 20+ yrs")

ts_symb_vec <- c(spy_sector_symbs, minerals_symbs, agriculture_symbs, energy_symbs,
                 currency_symbs, emerg_mkt_symbs, crypto_symbs, Tbond_symbs)
ts_detail_vec <- c(spy_sector_detail, minerals_detail, agriculture_detail, energy_detail,
                   currency_detail, emerg_mkt_detail, crypto_detail, Tbond_detail)
listGroupsDrivers <- list(spy_sector_symbs, minerals_symbs, agriculture_symbs, energy_symbs,
                          currency_symbs, emerg_mkt_symbs, crypto_symbs, Tbond_symbs)
groupNamesDrivers <- c("US Sectors", "Minerals", "Agriculture", "Energy", "Major Currencies",
                       "Emerging Markets", "Crypto", "T-Bonds")
names(listGroupsDrivers) <- groupNamesDrivers
nGroupsDrivers <- length(listGroupsDrivers)
bag_of_colors <- distinctColorPalette(k = 5 * nGroupsDrivers)
groupColorsDrivers <- sample(bag_of_colors, nGroupsDrivers)
groupInfoDrivers <- list(listGroupsDrivers, groupNamesDrivers, groupColorsDrivers)
dfDrivers <- data.frame(Ticker = ts_symb_vec, Name = ts_detail_vec)
dfDrivers$Sector <- NA;
for(i in 1:nGroupsDrivers){
  dfDrivers$Sector[which(dfDrivers$Ticker %in% listGroupsDrivers[[i]])] <- names(listGroupsDrivers)[i]
}
driversVec <- dfDrivers$Ticker
fromdate <- Sys.Date() - 1000
#=============================================================================
#=============================================================================
#=============================================================================
#=============================================================================
# Or use stocks from stock screener
# Tiingo has good free stock screener with option to export to csv:
# https://app.tiingo.com/screener/overview
#workFolder <- "D:/OneDrive - CGIAR/Documents 1/Personal stuff/quanFin/"
workFolder <- "/home/ben/Documents/finAnalysis/"
thisFile <-"LargeCap.csv"
#thisFile <-"midCap.csv"
thisFilepath <- paste0(workFolder, thisFile)
#list.files(workFolder)
dfThese <- read.csv(thisFilepath, stringsAsFactors = F) #Tiingo screener
theseStks <- unique(dfThese$Ticker)
unique(dfThese$Sector)
dfAvailRaw <- riingo::supported_tickers() %>% as.data.frame() #Get list of all available stocks
colnames(dfAvailRaw)[1] <- "Ticker"
dfAvail <- dfAvailRaw %>% merge(dfThese)
theseExchngs <- c("NYSE", "NASDAQ", "AMEX")
dfAvail <- dfAvail %>% subset(startDate < "2019-01-01" &
                                exchange %in% theseExchngs & 
                                endDate == (Sys.Date()))
allStks <- unique(dfAvail$Ticker)
length(allStks)
#allStks <- allStks[-which(allStks == "AHL-P-C")]
fromdate <- Sys.Date() - round(length(allStks) * 2)
length(fromdate:Sys.Date())
stkVec <- allStks
#=============================================================================
# Sort out sector info, especially for graphing purposes
unique(dfAvail$Sector)
dfAvail$Sector[grep("Materials", dfAvail$Sector)] <- "Materials"
dfAvail$Sector[grep("Tech", dfAvail$Sector)] <- "Technology"
dfAvail$Sector[grep("Unknown", dfAvail$Sector)] <- "Unknown"
sctrVec <- unique(dfAvail$Sector)
listGroups <- list()
for(i in 1:length(sctrVec)){
  listGroups[[sctrVec[i]]] <- dfAvail[, c("Ticker", "Sector")] %>%
    subset(Sector == sctrVec[i]) %>% .$Ticker
}
groupNames <- sctrVec
nGroups <- length(listGroups)
bag_of_colors <- distinctColorPalette(k = 5 * nGroups)
groupColors <- sample(bag_of_colors, nGroups)
groupInfo <- list(listGroups, groupNames, groupColors)
#----------------------------------------------------------------------------
# Add in drivers
colnames(dfAvail)
colnames(dfDrivers)
dfDrivers <- dfDrivers %>% mutate(exchange = NA,
                                  assetType = NA,
                                  priceCurrency = NA,
                                  startDate = NA,
                                  endDate = NA,
                                  Market.Cap = NA,
                                  Industry = NA,
                                  `X2.Year.Beta.to.S.P500` = NA,
                                  `P.E.Ratio` = NA,
                                  `Debt.Equity..D.E..Ratio` = NA)
dfDrivers <- dfDrivers[, colnames(dfAvail)]
dfAvail <- rbind(dfAvail, dfDrivers) %>% as.data.frame()
stkVec <- dfAvail$Ticker %>% unique()
# outlist <- group_fn(groupInfo)
# cols_ordered_by_group <- outlist[[1]]
# group_color_vec <- outlist[[2]]
# group_vec_ordered <- outlist[[3]]
# ind_ordered_cols <- outlist[[4]]
# df_match_group <- data.frame(Item = cols_ordered_by_group, Group = group_vec_ordered)
#=============================================================================
# Download price series
df_ohlcv <- stkVec %>% tq_get(get = "stock.prices", from = fromdate) %>% as.data.frame()
length(unique(df_ohlcv$symbol))
o <- apply(df_ohlcv, 2, function(x) sum(is.na(x))); o; table(o)
unique(df_ohlcv$symbol[which(is.na(df_ohlcv$low))])
df_ohlcv$p <- df_ohlcv$adjusted
#df_ohlcv[c(28305, 28306), ]
df_ohlcv <- df_ohlcv %>% group_by(symbol) %>%
  mutate(dup = duplicated(date)) %>% subset(dup == F)
dfx <- df_ohlcv[, c("symbol", "date", "p")] %>% spread(symbol, p)
o <- apply(dfx, 2, function(x) length(which(is.na(x)))); o; table(o)
which(o == max(o))
dfx[, -1] <- na.spline(dfx[, -1]) %>% as.data.frame()
dfx$LTM <- NULL
gathercols <- colnames(dfx)[-1]
df_ohlcv <- dfx %>% gather_("symbol", "p", gathercols)
#df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
#dfStk <- df_ohlcv[, c("symbol", "date", "p", "volume", "diffHiLo")]
#o <- apply(df_ohlcv, 2, function(x) length(which(is.na(x)))); o
# saveToFolder <- "D:/OneDrive - CGIAR/Documents 1/Personal stuff/quanFin/"
# fileName <- "largeCap.rds"
# saveFilePath <- paste0(saveToFolder, fileName)
# saveRDS(dfStk, saveFilePath)
#=============================================================================
# Get percentile oscillator series
dfStk <- df_ohlcv[, c("symbol", "date", "p")]
rollWind <- 55
dfStk <- dfStk %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  mutate(pctlOsc = 2 * pctlOsc - 1) %>% # For signal forcast/validation to work pctlOsc has to straddle y=0
  as.data.frame()
dfx <- dfStk[, c("symbol", "date", "pctlOsc")] %>% spread(symbol, pctlOsc)
o <- apply(dfx, 2, function(x) sum(is.na(x))); table(o)
notThese <- names(which(o >= rollWind))
dfStk <- dfStk %>% subset(!(symbol %in% notThese))
dfStk <- dfStk[-which(is.na(dfStk$pctlOsc)), c("symbol", "date", "pctlOsc")]
#=============================================================================
df_pca <- dfStk %>% spread(symbol, pctlOsc)
apply(df_pca, 2, function(x) length(which(is.na(x)))) %>% table()
mat_X_in <- as.matrix(df_pca[, -1]) %>% scale(scale = F)
row.names(mat_X_in) <- df_pca$date
apply(mat_X_in, 2, function(x) length(which(is.na(x)))) %>% table()
out <- get_S_and_corrXS(mat_X_in)
matS <- out[[1]]
mat_L <- out[[2]]
eigVals <- out[[3]]
matP <- out[[4]]
mat_Lrot <- varimax(mat_L)[[1]]
#matS <- mat_X_in %*% matP
pctExplnd <- cumsum(eigVals) / sum(eigVals)
cutOff1 <- which(pctExplnd >= 0.80)[1]
#cutOff2 <- which(colMeans(abs(mat_Lrot)) > 0.15)
cutOff2 <- which(eigVals / sum(eigVals) < 0.05)[1]
cutOff <- min(max(cutOff1), max(cutOff2)); cutOff
#cutOff <- 5
pctExplndVec <- round(100 * eigVals[1:cutOff] / sum(eigVals), 2)
mat_L <- mat_L[, 1:cutOff]
mat_Lrot <- mat_Lrot[, 1:cutOff]
sigNames <- paste0("PC ", 1:cutOff, "\n", pctExplndVec)
plot_corrXS_barchart(mat_Lrot, groupInfo, xAxis_title, sigNames)
#==========================================================================
# Plot signals (PCs) in time domain together with highest correlated stocks
stkSymbs <- row.names(mat_Lrot)
listDrivers <- list()
for(i in 1:cutOff){
  theseCorrs <- abs(mat_Lrot[, i])
  #ind <- which(theseCorrs > 0.6)
  ind <- which(theseCorrs > quantile(theseCorrs, probs = 0.98))
  #dfx <- dfAvail[, c("Ticker", "Name", "Sector", "Industry", "Market.Cap", "exchange")] %>% subset(Ticker %in% stkSymbs[ind])
  dfx <- dfAvail[, c("Ticker", "Name", "Sector")] %>% subset(Ticker %in% stkSymbs[ind])
  dfx$PC <- i; dfx$Lcorr <- mat_Lrot[ind, i]
  listDrivers[[i]] <- dfx
}
dfHiLcorStks <- as.data.frame(do.call(rbind, listDrivers)) %>% subset(abs(Lcorr) > 0.5) #merge(dfAvail) %>%
dfHiLcorStks$PC <- as.integer(dfHiLcorStks$PC)
#matSadj <- (mat_X_in %*% matP[, 1:cutOff] %*% diag(1 / sqrt(eigVals[1:cutOff]))) * 1 / 2
#dfS <- matSadj %>% as.data.frame()
dfS <- matS[, 1:cutOff] %>% as.data.frame()
#dfS$V3 <- -dfS$V3
dfS$date <- df_pca$date
dfS <- dfS %>% gather_("PC", "pctlOsc", colnames(dfS)[-ncol(dfS)])
dfS$PC <- as.integer(gsub("V", "", dfS$PC))
#-----------------------------------------------------------------------
listGg <- list()
for(i in 1:cutOff){
  theseHiCor <- dfHiLcorStks$Ticker[which(dfHiLcorStks$PC == i)]
  nStks <- length(theseHiCor)
  bag_of_colors <- distinctColorPalette(k = 5 * nStks)
  theseColors <- sample(bag_of_colors, nStks)
  dfTracks <- df_pca[, c("date", theseHiCor)] %>% gather_("ticker", "pctlOsc", theseHiCor)
  dfPlotS <- dfS %>% subset(PC == i)
  gg <- ggplot()
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red")
  gg <- gg + geom_line(data = dfPlotS, aes(x = date, y = pctlOsc), color = "grey", lwd = 1.3)
  gg <- gg + geom_line(data = dfTracks, aes(x = date, y = pctlOsc, group = ticker, color = ticker))
  gg <- gg + scale_color_manual(values = theseColors)
  #gg <- gg + ylim(-0.1, 1.1)
  gg <- gg + scale_y_continuous(breaks = c(-1, 0, 1))
  gg <- gg + theme_bw()
  gg <- gg + theme(axis.title = element_blank(),
                   legend.title = element_blank(),
                   axis.text = element_text(size = 7),
                   legend.text = element_text(size = 7))
  listGg[[i]] <- gg
}
wrap_plots(listGg, ncol = 1)
#=======================================================================
# Fit cyclic forecast model
# Validate and predict
#-----------------------------------------------------------------------
power_threshold <- "percentile log 0.1"#"drop lowest" #0.12
#-----------------------------------------------------------------------
# max_length <- 1000 # max length of time series
# rollWind <- 110
# per_ema_for_detrend <- rollWind
backtest_fraction <- 1 / 8
lookAhead_fraction <- 1 / 8
lookBack_fraction <- lookAhead_fraction * 0.5
#waveletComp params:
lowerPeriod <- 2^2 #2 * slope, #2^1
upperPeriod <- 2^10 #floor(nrow(df_wave)/3)*slope, #2^7
#-----------------------------------------------------------------------
slope <- 1 #(if daily data)
dj <- 1 / 250
this_period_label <- "period (days)"
# if(time_step_unit == "min"){
#   time_step_num <- as.numeric(stringr::str_extract_all(time_step, "[0-9]+")[[1]])
#   slope <- time_step_num / 60
#   dj <- 1 / 20
#   this_period_label <- "period (minutes)"
# }else{
#   slope <- 1 #(if daily data)
#   dj <- 1 / 250
#   this_period_label <- "period (days)"
# }
tsLength <- nrow(df_pca)
n_backtest = round(backtest_fraction * tsLength)
ind_fit <- 1:(tsLength - n_backtest)
ind_test <- setdiff(1:tsLength, ind_fit)
n_lookAhead <- round(lookAhead_fraction * tsLength)
##=========================================================================
# Start with look at all PCs (up to cutOff)
dfTs <- dfS; colnames(dfTs)[2:3] <- c("symbol", "ts"); dfTs$symbol <- paste("PC", dfTs$symbol)
##=========================================================================
# After identifying PCs of interest based on validation and prediction of them all,
# look at position/movement of highly correlated securities relative to selected PCs
PCfocus <- c(5)
dfTs <- dfS; colnames(dfTs)[2:3] <- c("symbol", "ts"); dfTs$symbol <- paste("PC", dfTs$symbol)
dfTs <- dfTs %>% subset(symbol %in% paste("PC", PCfocus))
dfTheseHiCor <- dfHiLcorStks[which(dfHiLcorStks$PC %in% PCfocus), c("Ticker", "PC")] %>%
  rename(symbol = Ticker)
theseHiCor <- dfTheseHiCor$symbol
dfTsAdd <- df_pca[, c("date", theseHiCor)]
if(length(theseHiCor) > 1){
  dfTsAdd <- dfTsAdd %>% gather_("symbol", "ts", theseHiCor)
}else{
  dfTsAdd$symbol <- theseHiCor; colnames(dfTsAdd)[2] <- "ts"
  dfTsAdd <- dfTsAdd[, c("date", "symbol", "ts")]
}
dfTsAdd <- dfTsAdd %>% merge(dfTheseHiCor) %>% rename(CorrWithPC = PC)
dfTs$CorrWithPC <- gsub("PC ", "", dfTs$symbol) %>% as.integer()
dfTs <- dfTs %>% rbind(dfTsAdd) %>% as.data.frame()
ind <- setdiff(c(1:nrow(dfTs)), grep("PC", dfTs$symbol))
dfTs$symbol[ind] <- paste(dfTs$symbol[ind], dfTs$CorrWithPC[ind])
#=========================================================================
# Fit model
symbVec <- unique(dfTs$symbol)
nTs <- length(symbVec)
listPCwaves <- list()
listPCwtImages <- list()
listDfMainCycles <- list()
listDfPeriodograms <- list()
listYhatValid <- list()
listYpredValid <- list()
listYhatPred <- list()
listYpredPred <- list()
for(i in 1:nTs){
  thisSymb <- symbVec[i]
  df_wave <- dfTs %>% subset(symbol == thisSymb)
  #------------------------------------------------------------------------
  waveAnalysis <- analyze.wavelet(df_wave, "ts",
                                  loess.span = 0,
                                  slope, #slope = time_step_num / 60 #if time_unit=min
                                  dj,
                                  lowerPeriod,
                                  upperPeriod, 
                                  make.pval = TRUE, n.sim = 10,
                                  verbose = F)
  #listPCwaves[[i]] <- waveAnalysis
  #------------------------------------------------------------
  # Plot beautiful waveComp analysis
  # wt.image(waveAnalysis, n.levels = 250, periodlab = "period (active minutes)",legend.params = list(lab = "wavelet power levels"), spec.time.axis = list(at = ind, labels = df_wave$date[ind]))
  wtImage <- wt.image(waveAnalysis, n.levels = 250, periodlab = this_period_label, legend.params = list(lab = "wavelet power levels", mar = 4.7))
  # my.rec <- reconstruct(my.w)
  # x.rec <- my.rec$series$x.r  # x: name of original series
  #listPCwtImages[[i]] <- wtImage
  #------------------------------------------------------------
  # Get periods and plot periodogram
  listOut <- getCycles(waveAnalysis, plotPeriodogram = F)
  dfMainCycles <- listOut[[1]]; dfMainCycles$symbol <- thisSymb
  dfPeriodograms <- listOut[[2]]; dfPeriodograms$symbol <- thisSymb
  listDfMainCycles[[i]] <- dfMainCycles
  listDfPeriodograms[[i]] <- dfPeriodograms
  #kable(round(dfMainCycles, 2)) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  #------------------------------------------------------------
  # Get info necessary for wave model validation (backtest)
  tsFit <- df_wave$ts[ind_fit]
  per_vec <- dfMainCycles$Period #if daily
  # if(time_step_unit == "min"){
  #   per_vec <- df_mainCycles$Period * 60 / time_step_num
  # }else{
  #   per_vec <- df_mainCycles$Period
  # }
  if(!is.null(power_threshold)){
    # Keep only periods above a certain power threshold
    power_threshold_input <- power_threshold
    power_threshold_str <- paste(as.character(stringr::str_extract_all(power_threshold, "[a-z]+")[[1]]), collapse = " ")
    if(power_threshold_str == "drop lowest"){
      power_threshold <- min(dfMainCycles$Power)
    }
    if(power_threshold_str == "percentile log"){
      power_threshold_num <- as.numeric(stringr::str_extract_all(power_threshold, "[0-9.]+")[[1]])
      x <- quantile(log(dfMainCycles$Power), power_threshold_num)
      power_threshold <- exp(x)
    }
    ind_keep <- which(dfMainCycles$Power > power_threshold)
    per_vec <- per_vec[ind_keep]
    
  }
  out_fitWave <- fitWave(tsFit, per_vec, pval_thresh = 0.01, n_backtest)
  yhat_validate <- out_fitWave[[1]]
  ypredict_validate <- out_fitWave[[2]]; ypredict_validate <- ypredict_validate[, 1]
  listYhatValid[[i]] <- yhat_validate
  listYpredValid[[i]] <- ypredict_validate
  #------------------------------------------------------------
  # Get info necessary for wave model prediction
  ts <- df_wave$ts
  out_fitWave <- fitWave(ts, per_vec, pval_thresh = 0.01, n_lookAhead)
  yhat_pred <- out_fitWave[[1]]
  ypredict_pred <- out_fitWave[[2]]; ypredict_pred <- ypredict_pred[, 1]
  summod <- out_fitWave[[3]]
  #kable(round(xtable(summod), 4)) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  listYhatPred[[i]] <- yhat_pred
  listYpredPred[[i]] <- ypredict_pred
  
}
#------------------------------------------------------------
dfYhatValid <- as.data.frame(do.call(cbind, listYhatValid))
dfYpredValid <- as.data.frame(do.call(cbind, listYpredValid))
dfYhatPred <- as.data.frame(do.call(cbind, listYhatPred))
dfYpredPred <- as.data.frame(do.call(cbind, listYpredPred))
colnames(dfYhatValid) <- symbVec; colnames(dfYpredValid) <- symbVec
colnames(dfYhatPred) <- symbVec; colnames(dfYpredPred) <- symbVec
#------------------------------------------------------------
# Display only ts with good validation
dateVec <- dfTs$date; lenTs <- length(dateVec)
lenTsValid <- nrow(dfYpredValid); dateVecValid <- dateVec[(lenTs - lenTsValid + 1):lenTs]
symbValidScore <- dfYpredValid %>% mutate(date = dateVecValid) %>%
  gather_("symbol", "tsPredValid", symbVec) %>%
  merge(dfTs[, c("date", "symbol", "ts")]) %>%
  group_by(symbol) %>% mutate(x = (ts - tsPredValid)^2) %>%
  summarise(sse = sum(x)) %>% subset(sse < quantile(sse, probs = 1)) %>%
  .$symbol;symbValidScore
dfYhatValidSel <- dfYhatValid[, symbValidScore];dfYpredValidSel <- dfYpredValid[, symbValidScore]
dfYhatPredSel <- dfYhatPred[, symbValidScore];dfYpredPredSel <- dfYpredPred[, symbValidScore]
dfTsSelect <- dfTs %>% subset(symbol %in% symbValidScore)
nTsSelect <- length(symbValidScore)
#------------------------------------------------------------
# Backtest
bag_of_colors <- distinctColorPalette(k = 5 * nTsSelect)
theseColors <- sample(bag_of_colors, nTsSelect)
plot_validation(dfYhatValidSel, dfYpredValidSel, dfTsSelect, colorVec = theseColors)
bag_of_colors <- distinctColorPalette(k = 5 * nTs)
theseColors <- sample(bag_of_colors, nTs)
plot_validation(dfYhatValid, dfYpredValid, dfTs, colorVec = theseColors)
#------------------------------------------------------------
# Prediction
n_lookAhead_zoom <- n_lookAhead
n_lookBack_zoom <- round(lookBack_fraction * nrow(df_wave))
time_step <- "daily"
plot_prediction(dfYhatPredSel, dfYpredPredSel, dfTsSelect, time_step, n_lookAhead, n_lookAhead_zoom, n_lookBack_zoom,
                colorVec = theseColors, showCloseUp = F)
plot_prediction(dfYhatPred, dfYpredPred, dfTs, time_step, n_lookAhead, n_lookAhead_zoom, n_lookBack_zoom,
                colorVec = theseColors, showCloseUp = F)
#------------------------------------------------------------
# Check fundamentals
dfLook <- dfAvail %>% subset(Ticker %in% gsub(" .*", "", symbValidScore))
View(dfLook)
#------------------------------------------------------------
# Periodograms
dfPeriodograms <- as.data.frame(do.call(rbind, listDfPeriodograms))
gg <- ggplot(dfPeriodograms, aes(x = period, y = power))
gg <- gg + geom_line()
gg <- gg + facet_wrap(~symbol, ncol = 1, strip.position = "right")
#gg <- gg + geom_vline(xintercept = critPers, color = "cyan", size = 1.2)
gg <- gg + theme_bw()
print(gg)
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
# End
#==============================================================




#

































# Examine individual signals
thisPC <- 5
dfPlotS <-dfS %>% subset(PC == thisPC)
gg <- ggplot(dfPlotS, aes(x = date, y = val))
gg <- gg + geom_hline(yintercept = c(0, 1 / 2, 1), color = "red")
gg <- gg + geom_line()
gg
#
theseHiCor <- dfHiLcorStks$Ticker[which(dfHiLcorStks$PC == thisPC)]
nStks <- length(theseHiCor)
bag_of_colors <- distinctColorPalette(k = 5 * nStks)
theseColors <- sample(bag_of_colors, nStks)
dfTracks <- df_pca[, c("date", theseHiCor)] %>% gather_("ts", "val", theseHiCor)
gg <- ggplot()
gg <- gg + geom_hline(yintercept = c(0, 1 / 2, 1), color = "red")
gg <- gg + geom_line(data = dfPlotS, aes(x = date, y = Signal, group = 1), color = "grey", lwd = 1.3)
gg <- gg + geom_line(data = dfTracks, aes(x = date, y = val, group = ts, color = ts))
gg <- gg + scale_color_manual(values = theseColors)
#gg <- gg + scale_x_discrete(breaks = xAxis_labels)
gg <- gg + theme_bw()
gg
#--------------------------------------------------------------------------
# Plot stock divergence from signal
dfPCDiv <- dfTracks %>% merge(dfPlotS) %>%
  mutate(PCdiv = Signal - val)
#gg <- ggplot(dfPCDiv, aes(x = date, y = PCdiv, group = ts, color = ts))
gg <- ggplot(dfPCDiv, aes(x = date, y = PCdiv))
gg <- gg + geom_hline(yintercept = 0, color = "red") + geom_line()
gg <- gg + geom_hline(yintercept = c(-1, 1), color = "red") + geom_line()
gg <- gg + facet_wrap(~ts)
gg <- gg + theme_bw()
gg
#==========================================================================
# Save to file
workFolder <- "/home/ben/Documents/finAnalysis/"
saveFile <-"hiLcor.csv"
saveFilepath <- paste0(workFolder, saveFile)
saveRDS(dfHiLcorStks, saveFilepath)
#write.csv(dfHiLcorStks, thisFilepath)

#==========================================================================


#

































#==========================================================================
#==========================================================================
#==========================================================================
# End
#==========================================================================
#==========================================================================
#==========================================================================
#==========================================================================
#plot_corrXS_barchart(mat_L, group_info = NULL, xAxis_title = NULL,
#                     sigNames = as.character(pctExplndVec))
# mat_L <- out[[2]]
# mat_Lrot <- varimax(mat_L)[[1]]
# mat_Lrot <- matrix(as.numeric(mat_Lrot),
#                    attributes(mat_Lrot)$dim,
#                    dimnames = attributes(mat_Lrot)$dimnames)
# mat_R <- varimax(mat_L)[[2]]
# mat_R <- matrix(as.numeric(mat_R),
#                 attributes(mat_R)$dim,
#                 dimnames = attributes(mat_R)$dimnames)

# dfBar <- mat_Lrot %>% as.data.frame()
# colnames(dfBar) <- sigNames
# gatherCols <- colnames(dfBar)
# dfBar$Ticker <- colnames(mat_X_in)
# dfBar <- dfBar %>% gather_("PC", "val", gatherCols)
# dfPlot <- dfBar %>% merge(dfAvail[, c("Ticker", "Sector", "Industry")])
# dfPlot$Sector <- factor(dfPlot$Sector, levels = unique(dfPlot$Sector))
# dfPlot$Ticker <- factor(dfPlot$Ticker, levels = unique(dfPlot$Ticker))
# 
# gg <- ggplot(dfPlot, aes(x = val, y = Ticker, fill = Sector))
# gg <- gg + geom_bar(stat = "identity", position = position_dodge(width = 0.8))
# gg <- gg + geom_vline(xintercept = 0, color = "red")
# gg <- gg + facet_wrap(~PC, nrow = 1)
# gg <- gg + theme_bw()
# gg


#df_pca <- dfFx[which(year(dfFx$date) > 2018), c("symbol", "date", "pctlOsc")]
#df <- as.data.frame(rbind(dfFx[, c("symbol", "date", "pctlOsc")], dfStk[, c("symbol", "date", "pctlOsc")]))
df <- dfStk[, c("symbol", "date", "pctlOsc")]
#df_pca <- df[which(year(df$date) > 2018), ]
df_pca <- df %>% spread(symbol, pctlOsc)
# Remove leading NAs if using percentile oscillator
# ind_rm <- 1:(rollWind - 1)
# df_pca <- df_pca[-ind_rm, ]
#---
o <- apply(df_pca, 2, function(x) length(which(is.na(x)))); o
ind_na <- which(is.na(df$p))
ind_na <- which(o > 0)
if(length(ind_na) != 0){
  df_pca[, ind_na] <- na.approx(df_pca[, ind_na])
}
#---
# Shorten names
# colnames(df_pca) <- gsub("ICE FUTURES U.S.", "ICE", colnames(df_pca))
# colnames(df_pca) <- gsub("CHICAGO BOARD OF TRADE", "CBoT", colnames(df_pca))
# colnames(df_pca) <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", colnames(df_pca))
# colnames(df_pca) <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", colnames(df_pca))
# colnames(df_pca) <- gsub("COMMODITY EXCHANGE INC.", "CE", colnames(df_pca))







































fxVec <- c("eurgbp", "eurusd",
           "audusd", "usdjpy",
           "audjpy", "btcusd")
#fx_vec2 <- c("eurjpy", "usdcad", "usdchf")
#stock_vec <- NULL
stkVec <- allStks
stock_vec <- c("XLF", "XLI", "XLK",
               "IYR", "XLV", "XLY", "XLP")
# stock_vec <- c("ES=F", "GC=F", "NG=F", "CL=F", "CT=F", "KC=F", "CC=F",
#                "SB=F", "ZB=F", "ZN=F", "ZF=F", "XLF", "XLI", "XLK",
#                "IYR", "XLV", "XLY", "XLP")
#max_length <- 2000
max_length <- round(length(stkVec) * 1.75)
#-----------------------------------------------------------------------------
# fx_emaPer_vec <- rep(34, length(fx_vec))
# stock_emaPer_vec <- rep(55, length(stock_vec))
# emaPer_vec <- c(fx_emaPer_vec, stock_emaPer_vec)
#-----------------------------------------------------------------------------
time_step <- "daily" #1min, 5min, 15min, 30min, 60min, daily, weekly
#=============================================================================
time_step_unit <- as.character(stringr::str_extract_all(time_step, "[a-z]+")[[1]])
time_step_num <- as.numeric(stringr::str_extract_all(time_step, "[0-9]+")[[1]])
if(time_step_unit == "min"){
  stkGetField <- "tiingo.iex"
  fxFreqField <- time_step
  fromdate <- Sys.Date() - max_length / (7.5 * (60 / time_step_num))
}else{
  stkGetField <- "tiingo"
  fxFreqField <- "1day"
  fromdate <- Sys.Date() - max_length
}
#=============================================================================
av_api_key("HQBAWJK4Y3YW81VG")
tiingo_api_key("36ed3f9ac9c2ba969b80c8389bc9a1e1fcdfcc42")
#=============================================================================
# Stocks
if(!is.null(stock_vec)){
  df_ohlcv <- stkVec %>%
    tq_get(get = stkGetField,
           from = fromdate, resample_frequency = time_step) %>%
    as.data.frame()
  #df_ohlcv$p <- df_ohlcv$adjusted
  df_ohlcv$p <- rowSums(df_ohlcv[, c(5:7)]) / 3
  df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
  dfStk <- df_ohlcv[, c("symbol", "date", "p", "volume", "diffHiLo")]
}
#----------------------------------------------------------------------------
# Forex or crypto pairs
if(!is.null(fxVec)){
  df_ohlcv <- fxVec %>%
    tq_get(get = "tiingo.crypto",
           from = fromdate, resample_frequency = fxFreqField) %>%
    as.data.frame()
  df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
  df_ohlcv$p <- rowSums(df_ohlcv[, c(5:7)]) / 3
  dfFx <- df_ohlcv[, c("symbol", "date", "p", "volume", "diffHiLo",
                       "volumeNotional", "tradesDone")]
}
#=============================================================================
# Check if there are any NAs
o <- apply(dfFx, 2, function(x) length(which(is.na(x)))); o
ind_na <- which(is.na(dfFx$p))
if(length(ind_na) != 0){
  df$p <- na.approx(df$p)
}
o <- apply(dfStk, 2, function(x) length(which(is.na(x)))); o
ind_na <- which(is.na(dfFx$p))
if(length(ind_na) != 0){
  df$p <- na.approx(df$p)
}
#-----------------------------------------------------------------------------
# Get percentile oscillator series
rollWind <- 89
dfStk <- dfStk %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
dfFx <- dfFx %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
dfStk <- dfStk[-which(is.na(dfStk$pctlOsc)), ]
dfFx <- dfFx[-which(is.na(dfFx$pctlOsc)), ]
#-----------------------------------------------------------------------------
# # Get p - ema oscillator
# dfStk <- dfStk %>% group_by(symbol) %>%
#   mutate(ema = EMA(p, per_ema_for_detrend)) %>%
#   mutate(difEMA = p - ema) %>%
#   as.data.frame()
# dfFx <- dfFx %>% group_by(symbol) %>%
#   mutate(ema = EMA(p, per_ema_for_detrend)) %>%
#   mutate(difEMA = p - ema) %>%
#   as.data.frame()
#-----------------------------------------------------------------------------
gg <- ggplot(dfFx, aes(x = date, y = pctlOsc))
gg <- gg + geom_line()
if(time_step_unit != "min"){
  dfFx$date <- as.Date(dfFx$date)
  gg <- gg + scale_x_date(breaks = scales::breaks_pretty(n = 4), labels = scales::date_format("%b\n%Y"))
}
gg <- gg + facet_wrap(~symbol, scales = "free_y")
gg

gg <- ggplot(dfStk, aes(x = date, y = pctlOsc))
gg <- gg + geom_line()
if(time_step_unit != "min"){
  gg <- gg + scale_x_date(breaks = scales::breaks_pretty(n = 4), labels = scales::date_format("%b\n%Y"))
}
gg <- gg + facet_wrap(~symbol, scales = "free_y")
gg
#-----------------------------------------------------------------------------




# # Read in the Commitment of Traders (CoT) data
# # Financials CoT data
# this_folder <- "D:/OneDrive - CGIAR/Documents 1/CIAT 2/finAnalysis/data/"
# this_filename <- "Finance CoT_processed.csv"
# this_filepath <- paste0(this_folder, this_filename)
# df_finCoT <- read.csv(this_filepath, stringsAsFactors = F)
# # Commodities CoT data
# this_filename <- "Commodity CoT_processed.csv"
# this_filepath <- paste0(this_folder, this_filename)
# df_comCoT <- read.csv(this_filepath, stringsAsFactors = F)
# #-----------------------------------------------------------------------------
# # Get rolling percentile series
# # Set rolling percentile window size for weekly CoT data
# rollWind_CoT <- 52
# #---
# df_pctlFinCoT <- subset(df_finCoT, Element == "Smart money net position (% of OI)")
# df_pctlFinCoT$Element <- NULL
# df_pctlFinCoT <- df_pctlFinCoT %>% group_by(Item) %>%
#   mutate(Value = rollapply(Value, rollWind_CoT, pctileFun, fill = NA, align = "right")) %>%
#   as.data.frame()
# df_pctlComCoT <- subset(df_comCoT, Element == "Smart money net position (% of OI)")
# df_pctlComCoT$Element <- NULL
# df_pctlComCoT <- df_pctlComCoT %>% group_by(Item) %>%
#   mutate(Value = rollapply(Value, rollWind_CoT, pctileFun, fill = NA, align = "right")) %>%
#   as.data.frame()
# #-----------------------------------------------------------------------------
# df_in <- df_pctlComCoT
# df_in <- df_in[-which(is.na(df_in$Value)), ]
# out_list <- getSeasons(df_in,
#                        freq = 52,
#                        mod_type = "additive",
#                        show_graphs = T)
# df_s <- out_list[[1]]
# ratio1 <- out_list[[2]]
# ratio2 <- out_list[[3]]

















# if(!is.null(stock_vec)){
#   if(time_step_unit == "min"){
#     # If intraday
#     # df_ohlcv <- stock_symbol %>%
#     #   tq_get(get = "alphavantage", av_fun = "TIME_SERIES_INTRADAY", interval = time_step, outputsize = "full") %>% as.data.frame()
#     # df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#     
#     fromdate <- Sys.Date() - max_length / (7.5 * (60 / time_step_num))
#     df_ohlcv <- stock_vec[1] %>% tq_get(get = "tiingo.iex",
#                                         from   = fromdate,
#                                         resample_frequency = time_step) %>% as.data.frame()
#     df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#     
#     
#   }
#   if(time_step_unit == "daily"){
#     # If daily
#     fromdate <- Sys.Date() - max_length
#     df_ohlcv <- stock_vec %>% tq_get(get = "stock.prices",
#                                      from = fromdate) %>% as.data.frame()
#     df_ohlcv$p <- df_ohlcv$adjusted
#     # df_ohlcv <- stock_symbol %>%
#     #   tq_get(get = "alphavantager", av_fun = "TIME_SERIES_DAILY_ADJUSTED", outputsize = "full") %>% as.data.frame(tbl_ohlcv)
#     # df_ohlcv$p <- df_ohlcv$adjusted_close
#   }
#   if(time_step_unit == "weekly"){
#     # If weekly
#     fromdate <- Sys.Date() - max_length * 5
#     df_ohlcv <- stock_vec %>% tq_get(get = "stock.prices",
#                                      from = fromdate,
#                                      periodicity = time_step_unit) %>% as.data.frame()
#     df_ohlcv$p <- df_ohlcv$adjusted
#     
#     # df_ohlcv <- stock_symbol %>%
#     #   tq_get(get = "alphavantager", av_fun = "TIME_SERIES_WEEKLY_ADJUSTED", outputsize = "full") %>% as.data.frame()
#     # df_ohlcv$p <- df_ohlcv$adjusted_close
#   }
# }
# # Forex
# # AlphaVantage allows only 5 queries per minute.
# # So split queries into two batches with minute wait in between.
# if(!is.null(fx_vec1)){
#   
#   # df_ohlcv1 <- fx_vec1 %>%  tq_get(get = "alphavantager",
#   #                                av_fun = "FX_WEEKLY",
#   #                                outputsize = "full") %>% as.data.frame()
#   # Sys.sleep(61)
#   # df_ohlcv2 <- fx_vec2 %>%  tq_get(get = "alphavantager",
#   #                                av_fun = "FX_WEEKLY",
#   #                                outputsize = "full") %>% as.data.frame()
#   # df_ohlcv <- as.data.frame(rbind(df_ohlcv1, df_ohlcv2))
#   # df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#   # df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
#   # colnames(df_ohlcv)[2] <- "date"
#   # dfFx <- df_ohlcv[, c("symbol", "date", "p", "diffHiLo")]
#   if(time_step_unit == "min"){
#     # If intraday
#     df_ohlcv1 <- fx_vec1 %>% tq_get(get = "alphavantage",
#                                     av_fun     = "TIME_SERIES_INTRADAY",
#                                     interval   = time_step,
#                                     outputsize = "full") %>% as.data.frame()
#     Sys.sleep(61)
#     df_ohlcv2 <- fx_vec2 %>% tq_get(get = "alphavantage",
#                                     av_fun     = "TIME_SERIES_INTRADAY",
#                                     interval   = time_step,
#                                     outputsize = "full") %>% as.data.frame()
#     df_ohlcv <- as.data.frame(rbind(df_ohlcv1, df_ohlcv2))
#     df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#     df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
#     colnames(df_ohlcv)[2] <- "date"
#     dfFx <- df_ohlcv[, c("symbol", "date", "p", "diffHiLo")]
#   }
#   if(time_step_unit == "daily"){
#     # If daily
#     df_ohlcv1 <- fx_vec1 %>% tq_get(get = "alphavantager", av_fun = "FX_DAILY", from_symbol = symb_currency_from, outputsize = "full") %>% as.data.frame()
#     df_ohlcv2 <- fx_vec2 %>% tq_get(get = "alphavantager", av_fun = "FX_DAILY", from_symbol = symb_currency_from, outputsize = "full") %>% as.data.frame()
#     df_ohlcv <- as.data.frame(rbind(df_ohlcv1, df_ohlcv2))
#     df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#     df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
#     colnames(df_ohlcv)[2] <- "date"
#     dfFx <- df_ohlcv[, c("symbol", "date", "p", "diffHiLo")]
#   }
#   if(time_step_unit == "weekly"){
#     # If weekly
#     df_ohlcv <- fx_vec1 %>% tq_get("", get = "alphavantager", av_fun = "FX_WEEKLY", from_symbol = symb_currency_from, to_symbol = symb_currency_to, outputsize = "full") %>% as.data.frame()
#     df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
#   }
#   df_ohlcv$symbol <- NULL
# }


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