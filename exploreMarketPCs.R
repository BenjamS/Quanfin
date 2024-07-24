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
  mat_S_all <- (mat_S_all + 1) * 1 / 2
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
                   panel.background = element_blank(),
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
  n_groups <- length(groupNames)
  n_items <- length(varNames_ordered)
  if(is.null(groupColors)){
    bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * n_groups)
    groupColors <- sample(bag_of_colors, n_groups)
    #group_colors <- viridis::viridis_pal(option = "D")(length(group_names))
  }
  #if(reverse_order){group_colors <- rev(group_colors)}
  #varNames_ordered <- colnames(mat_pctDiff)
  group_vec <- rep(NA, n_items)
  group_color_vec <- rep(NA, n_items)
  for(i in 1:n_groups){
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
plot_validation <- function(yhat, ypredict, df_wave){
  df_plot <- df_wave
  #df_plot <- df_wave[, c("date", "date_chr", "p", "ema", "slope")]
  # df_plot <- df_wave[, c("date", "date_chr", "p", "ema", "pctlOsc")]
  df_plot$yhat <- c(yhat, ypredict[, 1])
  ind_divide <- length(yhat)
  # # df_plot$set <- NA
  # # df_plot$set[ind_fit] <- "fit"
  # # df_plot$set[ind_test] <- "test"
  # df_plot$yhat_p <- df_plot$p + df_plot$yhat
  my_breaks <- df_plot$date_chr[seq.int(1, length(df_plot$date_chr), length.out = 30)]
  df_plot$date_chr <- NULL
  df_plot <- df_plot %>% gather(Type, Value, pctlOsc:yhat)
  df_plot$Value <- as.numeric(df_plot$Value)
  # #unique(df_plot$Type)
  # df_plot_ts <- subset(df_plot, Type %in% c("p", "ema", "yhat_p"))
  # #df_plot_ts_dtFit <- subset(df_plot, Type %in% c("slope", "yhat"))
  df_plot_ts_dtFit <- df_plot
  
  n_types_fit <- length(unique(df_plot_ts_dtFit$Type))
  #n_types_ts <- length(unique(df_plot_ts$Type))
  #distinct_colors <- randomcoloR::distinctColorPalette(k = n_types_fit + n_types_ts)
  distinct_colors <- randomcoloR::distinctColorPalette(k = 4 * n_types_fit)
  colors_dtFit <- distinct_colors[1:n_types_fit]
  #colors_ts <- distinct_colors[(n_types_fit + 1):(n_types_fit + n_types_ts)]
  
  df_plot_ts_dtFit$date_chr <- as.factor(as.character(df_plot_ts_dtFit$date))
  #df_plot_dt <- subset(df_plot_ts_dtFit, Type == "slope")
  df_plot_dt <- subset(df_plot_ts_dtFit, Type == "pctlOsc")
  df_plot_dtFit <- subset(df_plot_ts_dtFit, Type == "yhat")
  
  gg <- ggplot()
  gg <- gg + geom_line(data = df_plot_dtFit, aes(x = date_chr, y = Value, color = Type, group = 1), lwd = 1.1)#, color = colors_dtFit[2], lwd = 1.1)
  gg <- gg + geom_line(data = df_plot_dt, aes(x = date_chr, y = Value, color = Type, group = 1))#, color = colors_dtFit[1])
  gg <- gg + scale_color_manual(values = colors_dtFit)
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", lwd = 1)
  gg <- gg + geom_vline(aes(xintercept = ind_divide), lwd = 1, color = "violet")
  gg <- gg + theme_bw()
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + theme(axis.title = element_blank(),
                   legend.title = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1))
  #gg <- gg + scale_color_brewer(palette = "Dark2")
  gg_dtFit <- gg
  
  # df_plot_ts$date_chr <- as.factor(df_plot_ts$date_chr)
  # df_plot_yhat_p <- subset(df_plot_ts, Type == "yhat_p")
  # df_plot_p_ema <- subset(df_plot_ts, Type != "yhat_p")
  # gg <- ggplot()
  # gg <- gg + geom_line(data = df_plot_p_ema, aes(x = date_chr, y = Value, group = Type, color = Type))
  # gg <- gg + scale_color_manual(values = colors_ts[1:2])
  # gg <- gg + geom_line(data = df_plot_yhat_p, aes(x = date_chr, y = Value, group = 1), color = colors_ts[3], lwd = 1.1)
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
  print(gg_dtFit)
}
#========================================================================
plot_prediction <- function(yhat, ypredict, df_wave, time_step, n_lookAhead = 34, n_lookAhead_zoom = 21, n_lookBack_zoom = 21){
  ind_end <- nrow(df_wave)
  #df_plot <- df_wave[, c("date", "date_chr", "slope")]
  df_plot <- df_wave[, c("date", "date_chr", "pctlOsc")]
  df_plot$yhat <- yhat
  df_plot$set <- "fit"
  df_plot$t <- NULL
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
  # future_start <- df_plot$date[nrow(df_plot)] + time_step_num
  # future_stop <- df_plot$date[nrow(df_plot)] + 2 * time_step_num * n_lookAhead
  # date_fut <- seq(future_start, future_stop, by = time_step_num)
  # date_fut <- date_fut[!weekdays(date_fut) %in% c('Saturday','Sunday')]
  # date_fut <- date_fut[1:n_lookAhead]
  # date_fut_chr <- as.character(date_fut)
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
  #---------------------------------------
  #  df_add <- data.frame(date = date_fut, date_chr = date_fut_chr, slope = NA, yhat = ypredict[, 1], set = "predict")
  #df_add <- data.frame(date_chr = date_fut_chr, slope = NA, yhat = ypredict[, 1], set = "predict")
  df_add <- data.frame(date_chr = date_fut_chr, pctlOsc = NA, yhat = ypredict[, 1], set = "predict")
  df_plot$date <- NULL
  df_plot <- rbind(df_plot, df_add)
  df_plot$date_chr <- factor(df_plot$date_chr, levels = df_plot$date_chr)
  my_breaks <- df_plot$date_chr[seq.int(1, length(df_plot$date_chr), length.out = 30)]
  n <- 30
  distinct_colors <- randomcoloR::distinctColorPalette(k = n)
  colors_dtFit <- distinct_colors[sample(1:n, 2)]
  # df_plot_dt <- subset(df_plot, Type == "slope")
  # df_plot_dtFit <- subset(df_plot, Type == "yhat")
  this_title <- paste(c(paste(time_step, "chart")),
                        #paste(per_ema_for_detrend, "step detrend")),
                      collapse = ", ")  
  gg <- ggplot()
  gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = yhat, group = 1), color = colors_dtFit[1], lwd = 1.1)
  # gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = slope, group = 1), color = colors_dtFit[2])
  gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = pctlOsc, group = 1), color = colors_dtFit[2])
  gg <- gg + geom_vline(xintercept = ind_end, color = "blue", size = 1)
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", size = 1)
  gg <- gg + labs(title = this_title)
  gg <- gg + theme_bw()
  gg <- gg + theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1),
                   legend.title = element_blank(),
                   #legend.position = "top",
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(size = 10))
  #gg <- gg + scale_color_brewer(palette = "Dark2")
  print(gg)
  
  
  # Zoom in
  ind_end_new <- n_lookBack_zoom
  df_plot_zoom <- df_plot[(ind_end - n_lookBack_zoom):(ind_end + n_lookAhead_zoom), ]
  my_breaks_zoom <- df_plot_zoom$date_chr[seq.int(1, length(df_plot_zoom$date_chr), length.out = 30)]
  gg <- ggplot()
  gg <- gg + geom_line(data = df_plot_zoom, aes(x = date_chr, y = yhat, group = 1), color = colors_dtFit[1], lwd = 1.1)
  # gg <- gg + geom_line(data = df_plot_zoom, aes(x = date_chr, y = slope, group = 1), color = colors_dtFit[2])
  gg <- gg + geom_line(data = df_plot_zoom, aes(x = date_chr, y = pctlOsc, group = 1), color = colors_dtFit[2])
  gg <- gg + geom_vline(xintercept = ind_end_new, color = "blue", size = 1)
  gg <- gg + scale_x_discrete(breaks = my_breaks_zoom)
  gg <- gg + geom_hline(yintercept = c(-1, 0, 1), color = "red", size = 1)
  gg <- gg + labs(title = paste(this_title, "close up"))
  gg <- gg + theme_bw()
  gg <- gg + theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1),
                   legend.title = element_blank(),
                   #legend.position = "top",
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(size = 10))
  #gg <- gg + scale_color_brewer(palette = "Dark2")
  print(gg)
  
}
#========================================================================
get_cycles <- function(waveAnalysis){
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
  gg <- ggplot(df_plot, aes(x = period, y = power))
  gg <- gg + geom_line()
  gg <- gg + geom_vline(xintercept = critPers, color = "cyan", size = 1.2)
  gg <- gg + theme_bw()
  print(gg)
  #-------------------
  # Get periods, ordered from most power to least
  critPwrs <- df_periodogram$power[ind_critPoints]
  ind_order <- order(critPwrs, decreasing = T)
  per_vec <- critPers[ind_order]
  pwr_vec <- critPwrs[ind_order]
  df_mainCycles <- data.frame(Num = c(1:length(per_vec)), Period = per_vec, Power = pwr_vec)
  #-------------------
  return(df_mainCycles)
  
}
#========================================================================
#=============================================================================
#=============================================================================
# End function definition
#=============================================================================
#=============================================================================
# Tiingo has good free stock screener with option to export to csv:
# https://app.tiingo.com/screener/overview
#workFolder <- "D:/OneDrive - CGIAR/Documents 1/Personal stuff/quanFin/"
workFolder <- "/home/ben/Documents/finAnalysis/"
#thisFile <-"soundFundmntls.csv"
thisFile <-"overLeveraged.csv"
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
                                endDate == (Sys.Date() - 0))
#unique(dfAvailRaw$exchange)
#theseExchngs <- c("NYSE", "NASDAQ", "AMEX")
# dfAvail <- dfAvailRaw %>% subset(exchange %in% theseExchngs &
#                                    startDate < "2019-01-01" &
#                                    endDate == (Sys.Date() - 1) &
#                                    ticker %in% theseStks)
allStks <- unique(dfAvail$Ticker)
length(allStks)
#allStks <- allStks[-which(allStks == "AHL-P-C")]
fromdate <- Sys.Date() - round(length(allStks) * 12)
stkVec <- allStks
#=============================================================================
# Sort out sector info, especially for graphing purposes
unique(dfAvail$Sector)
dfAvail$Sector[grep("Materials", dfAvail$Sector)] <- "Materials"
dfAvail$Sector[grep("Tech", dfAvail$Sector)] <- "Technology"
dfAvail$Sector[grep("Unknown", dfAvail$Sector)] <- "Unknown"
sctrVec <- unique(dfAvail$Sector)
# list_groups <- list(spy_sector_detail, minerals_detail, agriculture_detail, energy_detail,
#                     currency_detail, emerg_mkt_detail, Tbond_detail) #crypto_detail,
# group_names <- c("US Sectors", "Minerals", "Agriculture", "Energy", "Major Currency Pairs",
#                  "Emerging Markets", "T-Bonds") #"Cryptocurrencies/\nBlockchain"
listGroups <- list()
for(i in 1:length(sctrVec)){
  listGroups[[sctrVec[i]]] <- dfAvail[, c("Ticker", "Sector")] %>%
    subset(Sector == sctrVec[i]) %>% .$Ticker
}
groupNames <- sctrVec
n_groups <- length(listGroups)
#group_colors <- RColorBrewer::brewer.pal(n = n_groups, "Dark2")
# "Darjeeling"
# group_colors <- wesanderson::wes_palette("Darjeeling1", n = n_groups, type = "continuous")
bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * n_groups)
groupColors <- sample(bag_of_colors, n_groups)
groupInfo <- list(listGroups, groupNames, groupColors)
# outlist <- group_fn(groupInfo)
# cols_ordered_by_group <- outlist[[1]]
# group_color_vec <- outlist[[2]]
# group_vec_ordered <- outlist[[3]]
# ind_ordered_cols <- outlist[[4]]
# df_match_group <- data.frame(Item = cols_ordered_by_group, Group = group_vec_ordered)
#=============================================================================
# Download price series
df_ohlcv  <- stkVec %>% tq_get(get = "stock.prices", from = fromdate) %>% as.data.frame()
length(unique(df_ohlcv$symbol))
o <- apply(df_ohlcv, 2, function(x) length(which(is.na(x)))); o
df_ohlcv$p <- df_ohlcv$adjusted
# dfx <- df_ohlcv[, c("symbol", "date", "p")] %>% spread(symbol, p)
# o <- apply(dfx, 2, function(x) length(which(is.na(x)))); o
df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
#dfStk <- df_ohlcv[, c("symbol", "date", "p", "volume", "diffHiLo")]
#o <- apply(df_ohlcv, 2, function(x) length(which(is.na(x)))); o
# saveToFolder <- "D:/OneDrive - CGIAR/Documents 1/Personal stuff/quanFin/"
# fileName <- "largeCap.rds"
# saveFilePath <- paste0(saveToFolder, fileName)
# saveRDS(dfStk, saveFilePath)
#=============================================================================
# Get percentile oscillator series
dfStk <- df_ohlcv[, c("symbol", "date", "p")]
rollWind <- 89
dfStk <- dfStk %>% group_by(symbol) %>%
  mutate(pctlOsc = rollapply(p, rollWind, pctileFun, fill = NA, align = "right")) %>%
  mutate(pctlOsc = 2 * pctlOsc - 1) %>%
  as.data.frame()
dfx <- dfStk[, c("symbol", "date", "pctlOsc")] %>% spread(symbol, pctlOsc)
o <- apply(dfx, 2, function(x) length(which(is.na(x)))); o;table(o)
max(o);which(o == max(o))
dfStk <- dfStk %>% subset(symbol != "ARKO")
dfStk <- dfStk %>% subset(symbol != "TECX")
dfStk <- dfStk[-which(is.na(dfStk$pctlOsc)), c("symbol", "date", "pctlOsc")]
o <- apply(dfStk, 2, function(x) length(which(is.na(x))));o
#=============================================================================
# df <- dfStk[, c("symbol", "date", "pctlOsc")]
# o <- apply(df, 2, function(x) length(which(is.nan(x))));o
df_pca <- dfStk %>% spread(symbol, pctlOsc)
o <- apply(df_pca, 2, function(x) length(which(is.na(x))));o;table(o)
mat_X_in <- as.matrix(df_pca[, -1]) %>% scale(scale = F)
row.names(mat_X_in) <- df_pca$date
o <- apply(mat_X_in, 2, function(x) length(which(is.na(x))));o;table(o)
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
cutOff <- min(max(cutOff1), max(cutOff2))
#if(cutOff > 7){cutOff <- 7}
#cutOff <- 7
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
  dfx <- dfAvail[, c("Ticker", "Name", "Sector", "Industry", "Market.Cap", "exchange")] %>% subset(Ticker %in% stkSymbs[ind])
  dfx$PC <- i; dfx$Lcorr <- mat_Lrot[ind, i]
  listDrivers[[i]] <- dfx
}
dfHiLcorStks <- as.data.frame(do.call(rbind, listDrivers)) %>% merge(dfAvail) %>%
  subset(Lcorr > 0.6)
dfHiLcorStks$PC <- as.integer(dfHiLcorStks$PC)
#matSadj <- (mat_X_in %*% matP[, 1:cutOff] %*% diag(1 / sqrt(eigVals[1:cutOff]))) * 1 / 2
#dfS <- matSadj %>% as.data.frame()
dfS <- matS[, 1:cutOff] %>% as.data.frame()
dfS$date <- df_pca$date
dfS <- dfS %>% gather_("PC", "pctlOsc", colnames(dfS)[-ncol(dfS)])
dfS$PC <- as.integer(gsub("V", "", dfS$PC))
#-----------------------------------------------------------------------
listGg <- list()
for(i in 1:cutOff){
  theseHiCor <- dfHiLcorStks$Ticker[which(dfHiLcorStks$PC == i)]
  nStks <- length(theseHiCor)
  bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * nStks)
  theseColors <- sample(bag_of_colors, nStks)
  dfTracks <- df_pca[, c("date", theseHiCor)] %>% gather_("ticker", "pctlOsc", theseHiCor)
  dfPlotS <-dfS %>% subset(PC == i)
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
# Forecast
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
i <- 3
df_wave <- dfS %>% subset(PC == i)
df_wave$date_chr <- as.character(df_wave$date)
#df_wave <- df[-c(1:(rollWind + 1)), c("date", "date_chr", "p", "ema", "pctlOsc")]
#------------------------------------------------------------------------
waveAnalysis <- analyze.wavelet(df_wave, "pctlOsc",
                                # waveAnalysis <- analyze.wavelet(df_wave, "slope",
                                loess.span = 0,
                                slope, #slope = time_step_num / 60 #if time_unit=min
                                dj,
                                lowerPeriod,
                                upperPeriod, 
                                make.pval = TRUE, n.sim = 10,
                                verbose = F)
#------------------------------------------------------------
# Plot beautiful waveComp analysis
# wt.image(waveAnalysis, n.levels = 250, periodlab = "period (active minutes)",legend.params = list(lab = "wavelet power levels"), spec.time.axis = list(at = ind, labels = df_wave$date[ind]))
wtImage <- wt.image(waveAnalysis, n.levels = 250, periodlab = this_period_label, legend.params = list(lab = "wavelet power levels", mar = 4.7))
# my.rec <- reconstruct(my.w)
# x.rec <- my.rec$series$x.r  # x: name of original series
#------------------------------------------------------------
# Get periods and plot periodogram
df_mainCycles <- get_cycles(waveAnalysis)
kable(round(df_mainCycles, 2)) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
#------------------------------------------------------------
# Get validate (backtest) info
n_backtest = round(backtest_fraction * nrow(df_wave))
ind_fit <- 1:(nrow(df_wave) - n_backtest)
ind_test <- setdiff(1:nrow(df_wave), ind_fit)
#ts <- df_wave$slope[ind_fit]
ts <- df_wave$pctlOsc[ind_fit]
per_vec <- df_mainCycles$Period #if daily
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
    power_threshold <- min(df_mainCycles$Power)
  }
  if(power_threshold_str == "percentile log"){
    power_threshold_num <- as.numeric(stringr::str_extract_all(power_threshold, "[0-9.]+")[[1]])
    x <- quantile(log(df_mainCycles$Power), power_threshold_num)
    power_threshold <- exp(x)
  }
  ind_keep <- which(df_mainCycles$Power > power_threshold)
  per_vec <- per_vec[ind_keep]
  
}
out_fitWave <- fitWave(ts, per_vec, pval_thresh = 0.01, n_backtest)
yhat_validate <- out_fitWave[[1]]
ypredict_validate <- out_fitWave[[2]]
#------------------------------------------------------------
# Get predict info
#ts <- df_wave$slope
ts <- df_wave$pctlOsc
n_lookAhead <- round(lookAhead_fraction * length(ts))
out_fitWave <- fitWave(ts, per_vec, pval_thresh = 0.01, n_lookAhead)
yhat_pred <- out_fitWave[[1]]
ypredict_pred <- out_fitWave[[2]]
summod <- out_fitWave[[3]]
kable(round(xtable(summod), 4)) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
#------------------------------------------------------------
# Backtest
plot_validation(yhat_validate, ypredict_validate, df_wave)
#------------------------------------------------------------
# Prediction
n_lookAhead_zoom <- n_lookAhead
n_lookBack_zoom <- round(lookBack_fraction * nrow(df_wave))
time_step <- "daily"
plot_prediction(yhat_pred, ypredict_pred, df_wave, time_step, n_lookAhead, n_lookAhead_zoom, n_lookBack_zoom)
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
bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * nStks)
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
