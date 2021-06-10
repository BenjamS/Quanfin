library(tidyverse)
library(tidyquant)
library(patchwork)
library(lubridate)
#--------------------------------------------------------------------------------
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
#=============================================================================
#=============================================================================
# Define fns
#=============================================================================
#=============================================================================
getTsTrendInds <- function(in_df, thresh_pct_uptrend = 0.7, thresh_pct_dntrend = -0.7){
  # Captures index of sustained up and down trends of a time series
  # by identifying periods of consecutive time steps in which slope
  # (i.e. the detrended series) is positve (for up trends) 
  # or negative (for down trends).
  # in_df = data.frame(date = etc., p = numeric, dt = numeric)
  # where dt equals the detrended series, usually made by subtracting a short period
  # (eg. 3-13) ema from the raw ts. This is essentially the slope (first derivative)
  # of the ts.
  #=======================================================
  #Initial capture based on slope crossings of the mean slope
  ind_uptrnd <- which(in_df$dt > 0)
  ind_dntrnd <- which(in_df$dt <= 0)
  #---------------------
  #Get trend start/finish points by taking differences
  ind_upFin <- ind_dntrnd[which(diff(ind_dntrnd) != 1) + 1]
  ind_upBeg <- ind_uptrnd[which(diff(ind_uptrnd) != 1) + 1]
  ind_dnFin <- ind_upBeg
  ind_dnBeg <- ind_upFin
  #If necessary, remove start/finish points so that you have one finish point
  #for every start point (i.e. a complete set)
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
  #=================================
  #Uptrends
  df_a <- in_df
  colnames(df_a)[2:ncol(df_a)] <- paste("Start", colnames(df_a)[2:ncol(df_a)])
  df_a$`Start date` <- NA
  df_a$`Start date`[ind_upBeg] <- as.character(df_a$date[ind_upBeg])
  df_a <- subset(df_a, is.na(`Start date`) != T)
  df_a$Date <- NULL
  df_b <- in_df
  colnames(df_b)[2:ncol(df_b)] <- paste("Stop", colnames(df_b)[2:ncol(df_b)])
  df_b$`Stop date` <- NA
  df_b$`Stop date`[ind_upFin] <- as.character(df_b$date[ind_upFin])
  df_b <- subset(df_b, is.na(`Stop date`) != T)
  df_b$Date <- NULL
  #--
  df_a$date <- NULL
  df_b$date <- NULL
  df_upTrends <- cbind(df_a, df_b)
  #--
  thresh_pct_uptrend <- 0
  df_upTrends$Change <- df_upTrends$`Stop p` - df_upTrends$`Start p`
  df_upTrends$`Pct. Change` <- 100 * df_upTrends$Change / df_upTrends$`Start p`
  df_upTrends$Duration <- ind_upFin - ind_upBeg
  #df_upTrends$`Pct. Change/Time` <- df_upTrends$`Pct. Change` / df_upTrends$Duration
  df_upTrends$`Pct. Change/Time` <- df_upTrends$`Pct. Change` * exp(-0.05 * df_upTrends$Duration)
  df_upTrends$`False trend` <- ifelse(df_upTrends$`Pct. Change` < thresh_pct_uptrend, 1, 0)
  #df_upTrends$`True trend` <- ifelse(df_upTrends$`Pct. Change` < thresh_pct_uptrend, 0, 1)
  df_upTrends$`Start date` <- as.Date(df_upTrends$`Start date`)
  df_upTrends$`Stop date` <- as.Date(df_upTrends$`Stop date`)
  #----------------------
  #Downtrends
  df_a <- in_df
  colnames(df_a)[2:ncol(df_a)] <- paste("Start", colnames(df_a)[2:ncol(df_a)])
  df_a$`Start date` <- NA
  df_a$`Start date`[ind_dnBeg] <- as.character(df_a$date[ind_dnBeg])
  df_a <- subset(df_a, is.na(`Start date`) != T)
  df_a$Date <- NULL
  df_b <- in_df
  colnames(df_b)[2:ncol(df_b)] <- paste("Stop", colnames(df_b)[2:ncol(df_b)])
  df_b$`Stop date` <- NA
  df_b$`Stop date`[ind_dnFin] <- as.character(df_b$date[ind_dnFin])
  df_b <- subset(df_b, is.na(`Stop date`) != T)
  df_b$Date <- NULL
  #--
  df_a$date <- NULL
  df_b$date <- NULL
  df_dnTrends <- cbind(df_a, df_b)
  #--
  df_dnTrends$Change <- df_dnTrends$`Stop p` - df_dnTrends$`Start p`
  df_dnTrends$`Pct. Change` <- 100 * df_dnTrends$Change / df_dnTrends$`Start p`
  df_dnTrends$Duration <- ind_dnFin - ind_dnBeg
  df_dnTrends$`Pct. Change/Time` <- df_dnTrends$`Pct. Change` / df_dnTrends$Duration
  df_dnTrends$`False trend` <- ifelse(df_dnTrends$`Pct. Change` > thresh_pct_dntrend, 1, 0)
  df_dnTrends$`Start date` <- as.Date(df_dnTrends$`Start date`)
  df_dnTrends$`Stop date` <- as.Date(df_dnTrends$`Stop date`)
  #=======================================================
  outlist <- list(df_upTrends, df_dnTrends)
  return(outlist)
}
#-----------------------------------------------------------------------------
# For consolidating more than two back-to-back trends
consolidateTrnds <- function(df, thresh_timeBetwn = 21, show_seqs = F){
  df$`TimeBetwn` <- c(NA, df$`Start date`[-1] - df$`Stop date`[-nrow(df)])
  ind <- which(df$TimeBetwn <= thresh_timeBetwn)
  #ind <- c(3, 5, 8, 15:17, 22, 26, 30:33, 36:37, 55, 59, 64, 73:78, 90, 93, 105)
  diffInd <- diff(ind)
  indDiffDiff1 <- which(diffInd == 1)
  seqUnbrok <- ind[indDiffDiff1]
  diffSeq <- diff(seqUnbrok)
  indBreaks <- which(diffSeq != 1)
  indBreaks <- c(indBreaks, length(seqUnbrok))
  seqEndPoints <- seqUnbrok[indBreaks]
  n_seq <- length(seqEndPoints)
  list_seq <- list()
  list_indrm <- list()
  this_break <- 0
  for(i in 1:n_seq){
    #this_endPoint <- seqEndPoints[i]
    indStart <- this_break + 1
    this_break <- indBreaks[i]
    this_seq <- seqUnbrok[indStart:this_break]
    this_seq <- c(this_seq, this_seq[length(this_seq)] + 1)
    list_seq[[i]] <- this_seq
    list_indrm[[i]] <- setdiff(this_seq, this_seq[length(this_seq)])
  }
  ind_rm <- do.call(c, list_indrm)
  
  df$Mark <- NA
  df$Mark[ind] <- 1
  df <- df[-ind_rm, ]
  ind <- which(df$Mark == 1)
  df$`Start date`[ind] <- df$`Start date`[ind - 1]
  df$`Start date_chr`[ind] <- df$`Start date_chr`[ind - 1]
  df$`Start p`[ind] <- df$`Start p`[ind - 1]
  df$`Pct. Change`[ind] <- 100 * (df$`Stop p`[ind] / df$`Start p`[ind] - 1)
  df$Duration[ind] <- df$`Stop date`[ind] - df$`Start date`[ind]
  df$`Pct. Change/Time`[ind] <- df$`Pct. Change`[ind] / df$Duration[ind]
  df <- df[-(ind - 1), ]
  
  if(show_seqs){
    print(list_seq)
  }
  
  return(df)
  
}
#-----------------------------------------------------------------------------
# For percentile oscillator
pctileFun <- function(x){
  out <- ecdf(x)(x[length(x)])
  return(out)
}
#-----------------------------------------------------------------------------
# Days (t-steps) above/below x percentile
daysAbovePctile <- function(df_pctile, thresh_pctiles = c(0.05, 0.95)){
  #pctileSeries <- df_plot$Pctile[which(!is.na(df_plot$Pctile))]
  #df_pctile <- df_plot
  ind_na <- which(is.na(df_pctile$Pctile))
  thresh_up <- thresh_pctiles[2]
  thresh_lo <- thresh_pctiles[1]
  #---
  indUp <- which(df_pctile$Pctile > thresh_up)
  breaks_indUp <- which(diff(indUp) != 1)
  breaks_indUp <- c(breaks_indUp, length(indUp))
  n_runs <- length(breaks_indUp)
  this_break <- 0
  #  df_pctile$`T-steps above thresh` <- NA
  df_pctile$`T-steps above thresh` <- 0
  df_pctile$`T-steps above thresh`[ind_na] <- NA
  for(i in 1:n_runs){
    indStart_indUp <- this_break + 1 
    this_break <- breaks_indUp[i]
    ind_thisRun <- indUp[indStart_indUp:this_break]
    numTstepsOverThresh <- 1:length(ind_thisRun)
    df_pctile$`T-steps above thresh`[ind_thisRun] <- numTstepsOverThresh
  }
  
  indLo <- which(df_pctile$Pctile < thresh_lo)
  breaks_indLo <- which(diff(indLo) != 1)
  breaks_indLo <- c(breaks_indLo, length(indLo))
  n_runs <- length(breaks_indLo)
  this_break <- 0
  #  df_pctile$`T-steps below thresh` <- NA
  df_pctile$`T-steps below thresh` <- 0
  df_pctile$`T-steps below thresh`[ind_na] <- NA
  for(i in 1:n_runs){
    indStart_indLo <- this_break + 1 
    this_break <- breaks_indLo[i]
    ind_thisRun <- indLo[indStart_indLo:this_break]
    numTstepsOverThresh <- 1:length(ind_thisRun)
    df_pctile$`T-steps below thresh`[ind_thisRun] <- numTstepsOverThresh
    
  }
  
  return(df_pctile)
  
}
#-----------------------------------------------------------------------------
# Visual inspection function
visuallyInspect <- function(df_plot, n_cols = 5){
  df_plot$date_chr <- as.character(df_plot$date)
  my_breaks <- df_plot$date_chr[seq.int(1, length(df_plot$date_chr), length.out = 5)]
  gg <- ggplot(df_plot, aes(x = date_chr, y = Value, group = 1))
  gg <- gg + geom_line()
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + facet_wrap(~Item, ncol = n_cols, scales = "free_y")
  gg <- gg + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  print(gg)
  
}
#-----------------------------------------------------------------------------
# Get principal components loadings function
get_S_and_corrXS <- function(mat_X_in){
  # mat_P = eigenvectors of the data correlation matrix
  # mat_G = corresponding eigenvalues
  mat_X_centered <- scale(mat_X_in, scale = F)
  # out_svd <- svd(mat_X_centered)
  # sing_values <- out_svd$d
  # n_obs <- nrow(mat_X_centered)
  # eig_values <- sing_values^2 / (n_obs - 1)
  # mat_P <- out_svd$v
  
  mat_X_in[which(is.na(mat_X_in[, 15])), 15]
  
  
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
# Define function to order data by group
group_fn <- function(group_info){
  list_groups <- group_info[[1]]
  group_names <- group_info[[2]]
  group_colors <- group_info[[3]]
  varNames_ordered <- do.call(c, list_groups)
  n_groups <- length(group_names)
  n_items <- length(varNames_ordered)
  if(is.na(group_colors)){
    bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * n_groups)
    group_colors <- sample(bag_of_colors, n_groups)
    #group_colors <- viridis::viridis_pal(option = "D")(length(group_names))
  }
  #if(reverse_order){group_colors <- rev(group_colors)}
  #varNames_ordered <- colnames(mat_pctDiff)
  group_vec <- rep(NA, n_items)
  group_color_vec <- rep(NA, n_items)
  for(i in 1:n_groups){
    this_group_vec <- list_groups[[i]]
    this_group_name <- group_names[i]
    this_group_color <- group_colors[i]
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
#=============================================================================
#=============================================================================
# End definition of functions
#=============================================================================
#=============================================================================
this_filename <- "fxFutData.csv"
this_filepath <- paste0(this_folder, this_filename)
df_raw <- read.csv(this_filepath, stringsAsFactors = F)
df_raw <- subset(df_raw, Item != "XBI")
df_raw$date <- as.POSIXct(df_raw$date, format = "%Y-%m-%d")
df_raw <- df_raw %>% group_by(Item, Element) %>%
  tq_transmute(select     = Value,
               mutate_fun = apply.weekly,
               FUN        = mean)
df_p <- subset(df_raw, Element == "p")
df_vol <- subset(df_raw, Element == "Volume")
df_dif <- subset(df_raw, Element == "diffHiLo")
df_p$Element <- NULL
df_vol$Element <- NULL
df_dif$Element <- NULL
#--------------------------------------------------------------------------------
# Detrend price
df_dt <- df_p %>% spread(Item, Value)
df_dt[, -1] <- as.data.frame(na.approx(df_dt[, -1]))
indNA <- which(is.na(df_dt$`CC=F`))
df_dt <- df_dt[-indNA, ]
per_ema_for_detrend <- 21
df_dt[, -1] <- as.data.frame(apply(df_dt[, -1], 2, EMA, per_ema_for_detrend))
df_dt <- df_dt[-c(1:(per_ema_for_detrend - 1)), ]
gathercols <- colnames(df_dt)[-1]
df_dt <- df_dt %>% gather_("Item", "ema", gathercols)
df_dt <- merge(df_dt, df_p, by = c("date", "Item"))
df_dt$Value <- df_dt$Value - df_dt$ema
df_dt$ema <- NULL
#--------------------------------------------------------------------------------
# Visually inspect
# Price
visuallyInspect(df_p, n_cols = 6)
# Price detrended
visuallyInspect(df_dt, n_cols = 6)
# High-low difference
visuallyInspect(df_dif, n_cols = 6)
# Volume
visuallyInspect(df_vol, n_cols = 4)
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Get trends
df <- subset(df_p, Item == "usd/eur")
df$date_chr <- as.character(df$date)
df <- df[, c("date", "date_chr", "Value")]
colnames(df)[ncol(df)] <- "p"
per_ema_for_detrend <- 21
df$ema <- EMA(df$p, per_ema_for_detrend)
df$dt <- df$p - df$ema

indList <- getTsTrendInds(df, thresh_pct_uptrend = 0.7, thresh_pct_dntrend = -0.7)
df_up <- indList[[1]]
df_down <- indList[[2]]
these_cols <- c("Start date_chr", "Stop date_chr", "Start date", "Stop date",
                "Start p", "Stop p", "Pct. Change", "Duration", "Pct. Change/Time", "False trend")
df_up <- subset(df_up[, these_cols], `False trend` == 0)
these_cols <- c("Start date_chr", "Stop date_chr", "Start date", "Stop date",
                "Start p", "Stop p", "Pct. Change", "Duration", "Pct. Change/Time", "False trend")
df_down <- subset(df_down[, these_cols], `False trend` == 0)
#-----------------------------------------------------------------------------
# Consolidate trends
thresh_timeBetwn <- 21
df_down <- consolidateTrnds(df_down, thresh_timeBetwn, show_seqs = F)
df_up <- consolidateTrnds(df_up, thresh_timeBetwn, show_seqs = F)
df_down <- subset(df_down, `Pct. Change` <= -2)
df_up <- subset(df_up, `Pct. Change` >= 2)
#-----------------------------------------------------------------------------
# Visually inspect trend capture
# hist(df_up$`Pct. Change`)
# hist(df_down$`Pct. Change`)
# hist(df_down$`Pct. Change/Time`)
# hist(df_up$`Pct. Change/Time`)
df_plot <- df[, c("date", "date_chr", "p", "ema")]
colnames(df_plot)[ncol(df_plot)] <- paste("ema", per_ema_for_detrend)
gathercols <- colnames(df_plot)[3:ncol(df_plot)]
df_plot <- df_plot %>% gather_("Type", "p", gathercols)
# Evenly spaced breaks
my_breaks <- df_plot$date_chr[seq.int(1, length(df_plot$date_chr), length.out = 30)]

gg <- ggplot()
# gg <- gg + geom_rect(data = df_upTrnds, aes(xmin = `Start date`, xmax = `Stop date`,
#                                             ymin = -Inf, ymax = Inf, fill = factor(`True uptrend`)), alpha = 0.7)
# gg <- gg + scale_fill_manual(values = c("magenta", "green"))
gg <- gg + geom_rect(data = df_up, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                       ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + geom_rect(data = df_down, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                         ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green")
gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = p, group = Type, color = Type))
gg <- gg + scale_x_discrete(breaks = my_breaks)
gg <- gg + theme(legend.position = "none",
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 legend.title = element_blank())
gg <- gg + guides(fill = guide_colorbar(title.position="top",
                                        title.hjust = 0.5),
                  color = "none")
gg_p <- gg
#-----------------------------------------------------------------------------
rollWind <- 144
df$Pctile <- rollapply(df$p, rollWind, pctileFun, fill = NA, align = "right")
df_pctile <- df[, c("date", "date_chr", "Pctile")]
thresh_pctiles <- c(0.05, 0.95)
df_pctileDays <- daysAbovePctile(df_pctile, thresh_pctiles)
df_pctileDays$Pctile <- NULL
#df_pctileDays <- df_pctileDays %>% gather(Item, Value, `T-steps above thresh`:`T-steps below thresh`)
df_pctileDays$`Pctile days above` <- rollapply(df_pctileDays$`T-steps above thresh`, rollWind, pctileFun, fill = NA, align = "right")
df_pctileDays$`Pctile days below` <- rollapply(df_pctileDays$`T-steps below thresh`, rollWind, pctileFun, fill = NA, align = "right")
df_pctileDays$`T-steps above thresh` <- NULL
df_pctileDays$`T-steps below thresh` <- NULL
#df_pctileDays <- df_pctileDays %>% gather(Item, Value, `Pctile days above`:`Pctile days below`)
df_pctileDays$pctileDiff <- df_pctileDays$`Pctile days below` - df_pctileDays$`Pctile days above`

df_plot <- df_pctile
gg <- ggplot()
gg <- gg + geom_rect(data = df_up, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                       ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + geom_rect(data = df_down, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                         ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green",
                                guide = "none")
gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = Pctile, group = 1))
gg <- gg + scale_x_discrete(breaks = my_breaks)
gg <- gg + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                 axis.title.x = element_blank(),
                 legend.position = "none")
gg_pctile <- gg

df_plot <- df_pctileDays
gg <- ggplot()
gg <- gg + geom_rect(data = df_up, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                       ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + geom_rect(data = df_down, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                         ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green",
                                guide = "none")
gg <- gg + geom_line(data = df_plot, aes(x = date_chr, y = pctileDiff, group = 1))#,
                                         #y = Value, group = Item, color = Item))
gg <- gg + geom_hline(yintercept = 0, color = "red")
gg <- gg + scale_color_manual(values = c("orange", "black"))
gg <- gg + scale_x_discrete(breaks = my_breaks)
gg <- gg + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.title = element_blank())
gg_pctileDays <- gg
#-----------------------------------------------------------------------------
gg_p + gg_pctileDays + gg_pctile + 
  plot_layout(ncol = 1, heights = c(1, 1 / 3, 1 / 3), guides = "collect") &
  theme(legend.position = "top")
#=============================================================================
rm(df_pctile, df_pctileDays)
#=============================================================================
#=============================================================================
# Calculate price percentile oscillator for all
rollWind <- 144
df_pctl <- df_p
df_pctl$Element <- NULL
df_pctl <- df_pctl %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
#-----------------------------------------------------------------------------
# Calculate days above/below threshold percentile oscillator for all
df_pctlDays <- df_pctl
df_pctlDays$Element <- NULL
colnames(df_pctlDays)[ncol(df_pctlDays)] <- "Pctile"
thresh_pctiles <- c(0.05, 0.95)
item_vec <- unique(df_pctlDays$Item); n_items <- length(item_vec)
list_df <- list()
for(i in 1:n_items){
  list_df[[i]] <- subset(df_pctlDays, Item == item_vec[i])
}
list_out <- lapply(list_df, daysAbovePctile, thresh_pctiles)
df_pctlDays <- as.data.frame(do.call(rbind, list_out))
df_pctlDays <- df_pctlDays %>% group_by(Item) %>%
  mutate(`Pctile days above` = rollapply(`T-steps above thresh`, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
df_pctlDays <- df_pctlDays %>% group_by(Item) %>%
  mutate(`Pctile days below` = rollapply(`T-steps below thresh`, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
df_pctlDays$Pctile <- NULL
df_pctlDays$`T-steps above thresh` <- NULL
df_pctlDays$`T-steps below thresh` <- NULL
df_pctlDays$pctlDaysDiff <- df_pctlDays$`Pctile days below` - df_pctlDays$`Pctile days above`
df_pctlDaysAbove <- df_pctlDays[, c("date", "Item", "Pctile days above")]
df_pctlDaysBelow <- df_pctlDays[, c("date", "Item", "Pctile days below")]
df_pctlDaysDif <- df_pctlDays[, c("date", "Item", "pctlDaysDiff")]
colnames(df_pctlDaysAbove)[3] <- "Value"
colnames(df_pctlDaysBelow)[3] <- "Value"
colnames(df_pctlDaysDif)[3] <- "Value"
#=============================================================================
# PCA
df_wide <- df_pctlDaysDiff %>% spread(Item, Value)
rows_rm <- 1:(rollWind - 1)
mat_X_in <- na.approx(df_wide[-rows_rm, -c(1)])
ind_rm <- which(is.na(mat_X_in[, 2]))
mat_X_in <- mat_X_in[-ind_rm, ]
# o <- apply(mat_X_in, 2, function(x) sum(is.na(x)))
# table(o)
list_out <- get_S_and_corrXS(mat_X_in)
mat_L <- list_out[[2]]
n_sigs <- 4
mat_L <- mat_L[, 1:n_sigs]

fx_vec <- item_vec[grep("/", unique(df_p$Item))]
us_vec <- c("ES=F", "ZB=F", "ZN=F", "ZF=F",
            "XLF", "XLI", "XLK", "IYR")
commod_vec <- c("GC=F", "NG=F", "CL=F", "CT=F", "KC=F", "CC=F", "SB=F")
list_groups <- list(fx_vec, us_vec, commod_vec)
group_names <- c("FX pairs", "US stocks\n& bonds", "Commodities")
n_groups <- length(list_groups)
bag_of_colors <- randomcoloR::distinctColorPalette(k = 5 * n_groups)
group_colors <- sample(bag_of_colors, n_groups)
group_info <- list(list_groups, group_names, group_colors)
plot_corrXS_barchart(mat_L, group_info = group_info)
# Varimax rotated loadings
mat_Lrot <- varimax(mat_L)[[1]]
mat_Lrot <- matrix(as.numeric(mat_Lrot),
                   attributes(mat_Lrot)$dim,
                   dimnames = attributes(mat_Lrot)$dimnames)
mat_R <- varimax(mat_L)[[2]]
mat_R <- matrix(as.numeric(mat_R),
                attributes(mat_R)$dim,
                dimnames = attributes(mat_R)$dimnames)

xAxis_title <- "Varimax Rotated Correlation"
plot_corrXS_barchart(mat_Lrot, group_info = group_info, xAxis_title, sigNames = NULL)
#=============================================================================
#=============================================================================
# Train ML model
#=============================================================================
#=============================================================================
df_upStart <- data.frame(date = df_up$`Start date`, y = "Buy")
df_upStop <- data.frame(date = df_up$`Stop date`, y = "Sell")
df_yUp <- as.data.frame(rbind(df_upStart, df_upStop))
df_yUp <- df_yUp[order(df_yUp$date), ]
df_yUp$y <- as.character(df_yUp$y)
df_dnStart <- data.frame(date = df_down$`Start date`, yDn = "Sell")
df_dnStop <- data.frame(date = df_down$`Stop date`, yDn = "Buy")
df_yDn <- as.data.frame(rbind(df_dnStart, df_dnStop))
df_yDn$yDn <- as.character(df_yDn$yDn)
df_yDn <- df_yDn[order(df_yDn$date), ]
#---
df_modP <- df_p %>% spread(Item, Value)
df_modP[, -1] <- na.approx(df_modP[, -1])
apply(df_modP[, -1], 2, function(x) sum(is.na(x)))
ind_rm <- which(is.na(df_modP[, which(colnames(df_modP) == "CC=F")]))
df_modP <- df_modP[-ind_rm, ]
colnames(df_modP)[-1] <- paste(colnames(df_modP)[-1], "p")
#-----------------------------------------------------------------------------
# Check to make sure the dependent var lines up with up and down trends
df_check <- plyr::join_all(list(df_modP, df_yUp, df_yDn), by = "date")
#df_check[which(!is.na(df_check$y) & !is.na(df_check$yDn)), ]
ind_dn <- which(!is.na(df_check$yDn))
df_check$y[ind_dn] <- df_check$yDn[ind_dn]
df_check$yDn <- NULL
df_plot <- subset(df_check[, c(1, which(colnames(df_check) == "usd/eur"))])
df_plot$date <- as.character(df_plot$date)
colnames(df_plot)[2] <- "p"
my_breaks <- df_plot$date[seq.int(1, length(df_plot$date), length.out = 20)]

gg <- ggplot()
gg <- gg + geom_rect(data = df_up, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                       ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + geom_rect(data = df_down, aes(xmin = `Start date_chr`, xmax = `Stop date_chr`,
                                         ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green")
gg <- gg + geom_line(data = df_plot, aes(x = date, y = p, group = 1))
gg <- gg + scale_x_discrete(breaks = my_breaks)
gg <- gg + theme(legend.position = "none",
                 axis.title.x = element_blank(),
                 axis.text.x = element_text(angle = 60, hjust = 1),
                 legend.title = element_blank())
gg <- gg + guides(fill = guide_colorbar(title.position="top",
                                        title.hjust = 0.5),
                  color = "none")
#---
df_plot <- subset(df_mod[, c(1, which(colnames(df_mod) %in% c("usd/eur", "y")))])
df_plot$date <- as.character(df_plot$date)
colnames(df_plot)[2] <- "p"
df_plot$colr <- NA
df_plot$colr[which(df_plot$y == "Buy")] <- "green"
df_plot$colr[which(df_plot$y == "Sell")] <- "red"
ind <- which(!is.na(df_plot$y))
df_plot$y[ind] <- df_plot$p[ind]
df_plot$y <- as.numeric(df_plot$y)
gg <- gg + geom_point(data = df_plot, aes(x = date, y = y, color = colr), size = 2)
gg <- gg + scale_color_manual(values = c("green", "red"))
gg
#-----------------------------------------------------------------------------
# Get feature datasets set up for input into model
df_pctlVol <- df_vol
df_pctlVol$Element <- NULL
df_pctlVol <- df_pctlVol %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
df_modVol <- df_pctlVol %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
df_modVol <- df_modVol[-ind_rm, ]
df_modVol[, -1] <- as.data.frame(na.approx(df_modVol[, -1], na.rm = F))
#df_modVol[, -1] <- log(df_modVol[, -1])
apply(df_modVol[, -1], 2, function(x) sum(is.na(x)))
ind_rm <- which(is.na(df_modVol[, which(colnames(df_modVol) == "CC=F")]))
df_modVol <- df_modVol[-ind_rm, ]
colnames(df_modVol)[-1] <- paste(colnames(df_modVol)[-1], "vol")
#---
df_pctlDif <- df_dif
df_pctlDif$Element <- NULL
df_pctlDif <- df_pctlDif %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
df_modDif <- df_pctlDif %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
df_modDif <- df_modDif[-ind_rm, ]
df_modDif[, -1] <- na.approx(df_modDif[, -1])
apply(df_modDif[, -1], 2, function(x) sum(is.na(x)))
ind_rm <- which(is.na(df_modDif[, which(colnames(df_modDif) == "CC=F")]))
df_modDif <- df_modDif[-ind_rm, ]
colnames(df_modDif)[-1] <- paste(colnames(df_modDif)[-1], "dif")
#---
df_modPctl <- df_pctl %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
df_modPctl <- df_modPctl[-ind_rm, ]
df_modPctl[, -1] <- na.approx(df_modPctl[, -1])
apply(df_modPctl[, -1], 2, function(x) sum(is.na(x)))
ind_rm <- which(is.na(df_modPctl[, which(colnames(df_modPctl) == "CC=F")]))
df_modPctl <- df_modPctl[-ind_rm, ]
colnames(df_modPctl)[-1] <- paste(colnames(df_modPctl)[-1], "pctl")
#---
df_modPctlDaysDif <- df_pctlDaysDif %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
df_modPctlDaysDif <- df_modPctlDaysDif[-ind_rm, ]
df_modPctlDaysDif[, -1] <- na.approx(df_modPctlDaysDif[, -1])
apply(df_modPctlDaysDif[, -1], 2, function(x) sum(is.na(x)))
ind_rm <- which(is.na(df_modPctlDaysDif[, which(colnames(df_modPctlDaysDif) == "CC=F")]))
df_modPctlDaysDif <- df_modPctlDaysDif[-ind_rm, ]
colnames(df_modPctlDaysDif)[-1] <- paste(colnames(df_modPctlDaysDif)[-1], "pctlDaysDif")
#---
df_modP[, -1] <- log(df_modP[, -1])
#---
#list_df <- list(df_modDif, df_modPctl, df_modPctlDaysDif, df_yUp, df_yDn)
#-----------------------------------------------------------------------------
# Which features to include
#list_df <- list(df_modP, df_modDif, df_modPctl, df_yUp, df_yDn)
list_df <- list(df_modPctl, df_modVol, df_modDif, df_modPctlDaysDif, df_yUp, df_yDn)
#-----------------------------------------------------------------------------
df_mod <- plyr::join_all(list_df, by = "date")
this_condition <- (length(grep(" p| vol| dif", colnames(df_mod))) != 0 &
                     length(grep(" pctl| pctlDaysDif| pctlDaysAbove| pctlDaysBelow", colnames(df_mod))))
if(this_condition){
  ind_rm <- 1:(rollWind - 1)
  df_mod <- df_mod[-ind_rm, ]
}
#df_mod[which(!is.na(df_mod$y) & !is.na(df_mod$yDn)), ]
ind_dn <- which(!is.na(df_mod$yDn))
df_mod$y[ind_dn] <- df_mod$yDn[ind_dn]
df_mod$yDn <- NULL
df_mod$y[which(is.na(df_mod$y))] <- "Hold"
#---
numNA_vec <- apply(df_mod[, -1], 2, function(x) sum(is.na(x)))
table(numNA_vec)
#df_mod$date[which(is.na(df_mod$`CL=F p`))]
df_mod[, -c(1, ncol(df_mod))] <- na.approx(df_mod[, -c(1, ncol(df_mod))])
#---
# How correlated are the features?
df_cor <- df_mod %>% select(-date, -y) %>% cor() %>% as.data.frame()
#corrplot::corrplot(method = 'circle',tl.pos='n')
#---
df_mod$y <- as.factor(df_mod$y)
df_mod$y <- relevel(df_mod$y, ref = "Hold")
contrasts(df_train$y)

trainDat_pctot <- .7
indtrain_beg <- 1
indtrain_end <- round(nrow(df_mod) * trainDat_pctot)
indtest_beg <- indtrain_end + 1
indtest_end <- nrow(df_mod)
indtrain <- indtrain_beg:indtrain_end
indtest <- indtest_beg:indtest_end
df_train <- df_mod[indtrain, ]
df_test <- df_mod[indtest, ]

# nrow(df_mod)
# nrow(df_train)
# nrow(df_test)
# Multinomial logistic regression using nnet package
#https://medium.com/@PAdhokshaja/using-anova-and-multinomial-logistic-regression-to-predict-human-activity-cd2101a5e8bf
df_train$date <- NULL
model <- nnet::multinom(y ~ ., data = df_train, maxit = 400)
#summary(model)
#library(AER)
#coeftest(model)
library(car)
Anova(model)
# df_broom <- as.data.frame(broom::tidy(model))
# df_broom[which(df_broom$p.value < 0.05), ]
#Prediction
mat_pred <- predict(model, df_test, type = "class")#type = "probs")
# ind_na <- which(is.na(mat_pred))
# ind_na
# mat_pred[ind_na] <- "Hold"
unique(mat_pred)
df_pred <- as.data.frame(mat_pred)
colnames(df_pred)[1] <- "yPred"
df_pred$yObs <- df_test$y
#misClasificError <- mean(df_pred$yPred != df_pred$yObs)
library(caret)
postResample(df_pred$yObs, df_pred$yPred)
# Confusion matrix
x <- confusionMatrix(df_pred$yPred, df_pred$yObs)
confmat <- x$table
confmat <- round(confmat %*% solve(diag(colSums(confmat))), 3)
confmat <- as.table(confmat)
colnames(confmat) <- rownames(confmat)
names(dimnames(confmat))[2] <- "Reference"
print(confmat)
print(x$table)
#class(confmat)
df_plot <- as.data.frame(confmat)
gg <- ggplot(df_plot) + geom_tile(aes(x = Prediction, y = Reference, fill = Freq))
gg <- gg + scale_fill_gradient(low = "orange", high = "cyan")
gg
#-----------------------------------------------------------------------------
# If binary model
# model <- glm(y ~., family = binomial(link = 'logit'), data = df_train)
# summary(model)
# anova(model, test="Chisq")
# fitted.prob <- predict(model, newdata = df_test, type = 'response')
# fitted.binry <- ifelse(fitted.prob > 0.6, 1, 0)
# df_compare <- data.frame(predicted_prob = fitted.prob, predicted_binry = fitted.binry, observed = df_test$y)
# df_compare$observed <- ifelse(df_test$y == "Start", 1, 0)
# misClasificError <- mean(df_compare$predicted_binry != df_compare$observed)
# print(paste('Accuracy', round(1 - misClasificError, 2)))
# misClasificError <- mean(abs(df_compare$observed - df_compare$predicted_prob))
# print(paste('Accuracy (prob)', round(misClasificError, 2)))
#-----------------------------------------------------------------------------
df_compare$Date <- df_ML_date[indtest]
df_plot <- fortify(xts_cp_mat[, this_ts_name])
colnames(df_plot) <- c("Date", "cp")
df_plot <- left_join(df_plot, df_compare)
indtest_cp <- c(which(df_plot$Date == df_compare$Date[1]):nrow(df_plot))
df_plot <- df_plot[indtest_cp,]
df_plot_up_true <- subset(df_upTrends, `False uptrend` == 0)
df_plot_up_false <- subset(df_upTrends, `False uptrend` == 1)
n_bins <- length(indtest_cp)
df_probGradient <- data.frame(xmin = df_plot$Date[-n_bins], xmax = df_plot$Date[-1])
df_probGradient$predProbs <- df_plot$predicted_prob[-1]
u <- df_probGradient$predProbs
df_probGradient$predProbs[which(is.na(u))] <- 0.5
gg <- ggplot()
gg <- gg + geom_rect(data = df_probGradient, aes(xmin = xmin, xmax = xmax,
                                                 ymin = -Inf, ymax = Inf, fill = predProbs), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "green", midpoint = 0.5)
gg <- gg + geom_line(data = df_plot, aes(x = Date, y = cp))
gg <- gg + theme_bw()
gg <- gg + ggtitle("Backtest: Predicted start in green, false start in red. Shade = confidence.")
gg <- gg + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_text(hjust = 0.5))
gg
gg_predicted <- gg


indtest_trends <- which(df_upTrends$UpStartDate %in% df_plot$Date)
df_upTrends_test <- df_upTrends[indtest_trends,]
gg <- ggplot()
gg <- gg + geom_rect(data = df_upTrends_test, aes(xmin = UpStartDate, xmax = UpStopDate,
                                                  ymin = -Inf, ymax = Inf, fill = `Pct. Change/Time`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green")
gg <- gg + geom_line(data = df_plot, aes(x = Date, y = cp))
gg <- gg + scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")
gg <- gg + geom_vline(data = df_plot_up_false, aes(xintercept = UpStartDate), color = "darkmagenta", alpha = 0.4)
gg <- gg + geom_vline(data = df_plot_up_false, aes(xintercept = UpStopDate), color = "darkmagenta", alpha = 0.4)
#gg <- gg + geom_vline(data = df_compare, aes(xintercept = Date), color = "green", size = 0.5, linetype = "dotted")
# gg <- gg + geom_hline(data = df_plot_mean_line, aes(yintercept = mu_line), color = "blue")
# gg <- gg + geom_hline(data = df_plot_sd_lines1, aes(yintercept = sd_line), color = "orange")
# gg <- gg + geom_hline(data = df_plot_sd_lines2, aes(yintercept = sd_line), color = "orange")
gg <- gg + theme_bw()
gg <- gg + ggtitle("Uptrends marked in green, false trends in red. Shade = return intenstity")
gg <- gg + theme(axis.title.x = element_blank(),
                 # axis.text.x = element_blank(),
                 # axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_text(hjust = 0.5))
gg
gg_observed <- gg
#library(gridExtra)
grid.arrange(gg_predicted, gg_observed)
