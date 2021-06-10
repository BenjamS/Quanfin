#https://online.stat.psu.edu/stat510/lesson/8/8.2
library(tidyverse)
library(tidyquant)
library(patchwork)
library(lubridate)
library(pracma)
library(caret)
library(pls)
library(car)
#--------------------------------------------------------------------------------
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
#=============================================================================
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
  
  if("Element" %in% colnames(df_plot)){
      gg <- ggplot(df_plot, aes(x = date_chr, y = Value,
                                group = Element, color = Element))
  }else{
     gg <- ggplot(df_plot, aes(x = date_chr, y = Value, group = 1))
  }
  gg <- gg + geom_line()
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + facet_wrap(~Item, ncol = n_cols, scales = "free_y")
  gg <- gg + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "top")
  print(gg)
  
}

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
  ind_neg <- which(eig_values < 0)
  if(length(ind_neg != 0)){
    for(i in 1:length(ind_neg)){
      if(eig_values[ind_neg[i]] < -10^-4){
        print("Problem: large negative eigenvalue")
      }
        
    }
    eig_values <- abs(eig_values)
  }
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

#-----------------------------------------------------------------------------
getTsTrendInds <- function(in_df,
                           pctThresh_uptrend = 1,
                           pctThresh_dntrend = -1){
  # Captures start and stop indices of up and down trends of a time series.
  # Columns of in_df must be "date", "targ series"
  #=======================================================
  #Get index of peaks and valleys
  ind_critPos <- pracma::findpeaks(in_df$`targ series`)[, 2]
  ind_critNeg <- pracma::findpeaks(-in_df$`targ series`)[, 2]
  #Get trend start/finish points by taking differences
  ind_upFin <- ind_critPos[which(diff(ind_critPos) != 1) + 1]
  ind_upBeg <- ind_critNeg[which(diff(ind_critNeg) != 1) + 1]
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
  df_upTrends$Change <- df_upTrends$`Stop targ series` - df_upTrends$`Start targ series`
  df_upTrends$`Pct. Change` <- 100 * df_upTrends$Change / df_upTrends$`Start targ series`
  df_upTrends$Duration <- ind_upFin - ind_upBeg
  df_upTrends$`Pct. Change/Time` <- df_upTrends$`Pct. Change` * exp(-0.05 * df_upTrends$Duration)
  df_upTrends$`False trend` <- ifelse(df_upTrends$`Pct. Change` < pctThresh_uptrend, 1, 0)
  # df_upTrends$`Start date` <- as.Date(df_upTrends$`Start date`)
  # df_upTrends$`Stop date` <- as.Date(df_upTrends$`Stop date`)
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
  df_dnTrends$Change <- df_dnTrends$`Stop targ series` - df_dnTrends$`Start targ series`
  df_dnTrends$`Pct. Change` <- 100 * df_dnTrends$Change / df_dnTrends$`Start targ series`
  df_dnTrends$Duration <- ind_dnFin - ind_dnBeg
  df_dnTrends$`Pct. Change/Time` <- df_dnTrends$`Pct. Change` / exp(-0.05 * df_dnTrends$Duration)
  df_dnTrends$`False trend` <- ifelse(df_dnTrends$`Pct. Change` > pctThresh_dntrend, 1, 0)
  # df_dnTrends$`Start date` <- as.Date(df_dnTrends$`Start date`)
  # df_dnTrends$`Stop date` <- as.Date(df_dnTrends$`Stop date`)
  #------------------------------------------------------------------------
  list_out <- list(df_upTrends, df_dnTrends)
  return(list_out)
  
}
#-----------------------------------------------------------------------------
# For consolidating more than two back-to-back trends
consolidateTrnds <- function(df, thresh_timeBetwn = 21, show_seqs = F){
  # If time between two trends is less than thresh_timeBetwn, they are
  # stitched together. The threshold thresh_timeBetwn must be defined in
  # days.
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

#---------------------------------------------------------------------------
# Get df of independent variables lagged at lags that are
# correlated with target series 
get_correlatedLaggedSeries <- function(df_ccf, max_lag = 25, max_n_lags = NULL){
  # df_ccf is the df of series that are to be compared with the target series
  # df_targMo is the df of the target series
  # Both dfs have a date vec. df_ccf is wide.
  list_bestLags <- list()
  list_bestCors <- list()
  list_dfLags <- list()
  targ_series <- df_ccf$targ_series
  date_vec <- df_ccf$date
  df_ccf$date <- NULL
  for(i in 1:ncol(df_ccf)){
    thisOther_guy <- colnames(df_ccf)[i]
    print(thisOther_guy)
    thisOther_vec <- df_ccf[, i]
    this_title <- paste("target", "&", thisOther_guy)
    x <- ccf(targ_series, thisOther_vec, main = this_title, max_lag)
    Sys.sleep(1)
    cor_vec <- as.numeric(x[[1]])
    lag_vec <- as.numeric(x[[4]])
    ind_critPos <- findpeaks(cor_vec)[, 2]
    ind_critNeg <- findpeaks(-cor_vec)[, 2]
    ind_critPeaks <- c(ind_critNeg, ind_critPos)
    N <- length(cor_vec)
    sigThresh95 <- exp(2*1.96/sqrt(N-3)-1)/exp(2*1.96/sqrt(N-3)+1)
    sigThresh98 <- sigThresh95 * 1.5
    
    ind_crit98 <- which(cor_vec > sigThresh98 | cor_vec < -sigThresh98)
    
    ind_indCritPeaks <- which(cor_vec[ind_critPeaks] > sigThresh95 | cor_vec[ind_critPeaks] < -sigThresh95)
    ind_crit <- unique(c(ind_critPeaks[ind_indCritPeaks], ind_crit98))
    
    if(length(ind_crit) != 0){
      ind_negLag <- which(lag_vec <= 0)
      cor_vec <- cor_vec[ind_negLag]
      lag_vec <- lag_vec[ind_negLag]
      #plot(cor_vec[ind_crit])
      plot(lag_vec, cor_vec,type='l', main = this_title)
      points(lag_vec[ind_crit], cor_vec[ind_crit])
      abline(h = sigThresh95, lty = "dashed", col = "blue")
      abline(h = -sigThresh95, lty = "dashed", col = "blue")
      abline(h = 0, lty = "dotted", col = "red")
      
      Sys.sleep(1)
      ind_crit <- ind_crit[order(ind_crit)]
      these_lags <- as.numeric(na.omit(lag_vec[ind_crit]))
      these_cors <- as.numeric(na.omit(cor_vec[ind_crit]))
      list_bestLags[[thisOther_guy]] <- these_lags
      list_bestCors[[thisOther_guy]] <- these_cors
      #---
      # Put lags in order of highest correlation to least
      these_lags <- these_lags[order(these_cors, decreasing = T)]
      these_cors <- these_cors[order(these_cors, decreasing = T)]
      #---
      if(length(these_lags) != 0){
        n_lags <- length(these_lags)
        if(!is.null(max_n_lags)){
          n_lags <- min(max_n_lags, n_lags)
        }
        list_lagged <- list()
        for(j in 1:n_lags){
          list_lagged[[j]] <- Hmisc::Lag(thisOther_vec, these_lags[j])
        }
        df_lags <- as.data.frame(do.call(cbind, list_lagged))
        colnames(df_lags) <- paste(thisOther_guy, these_lags[1:n_lags])
        df_lags$date <- date_vec
        df_lags <- df_lags[, c("date", colnames(df_lags)[-ncol(df_lags)])]
        list_dfLags[[thisOther_guy]] <- df_lags
      }else{
        print("No significant cross correlation.")
      }
    }else{
      print("No significant cross correlation.")
    }
  }
  # all_lags <- do.call(c, list_bestLags)
  # max_lag <- min(all_lags)
  # hist(all_lags)
  # unique(all_lags)
  df_lagged <- plyr::join_all(list_dfLags, by = "date")
  df_lagged$`targ series 0` <- NULL
  df_lagged$`targ series` <- targ_series
  return(df_lagged)
}

#=============================================================================
#=============================================================================
# End function definition
#=============================================================================
#=============================================================================
# Import the price/vol/Hi-Lo difference data
this_filename <- "fxFutData.csv"
this_filepath <- paste0(this_folder, this_filename)
df_raw <- read.csv(this_filepath, stringsAsFactors = F)
df_raw <- subset(df_raw, Item != "XBI")
# Change to weekly/monthly
# df_raw$date <- as.POSIXct(df_raw$date, format = "%Y-%m-%d")
# df_raw <- df_raw %>% group_by(Element, Item) %>%
#   tq_transmute(select     = Value,
#                mutate_fun = apply.weekly,
#                FUN        = mean,
#                na.rm = T) %>%
#   as.data.frame()
df_pvhl <- df_raw
#-----------------------------------------------------------------------------
# Import the Commitment of Traders (CoT) data
# Financials CoT data
this_filename <- "Finance CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_finCoT <- read.csv(this_filepath, stringsAsFactors = F)
# Commodities CoT data
this_filename <- "Commodity CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_comCoT <- read.csv(this_filepath, stringsAsFactors = F)
#-----------------------------------------------------------------------------
df_p <- subset(df_pvhl, Element == "p")
df_vol <- subset(df_pvhl, Element == "Volume")
df_dif <- subset(df_pvhl, Element == "diffHiLo")
df_p$Element <- NULL
df_vol$Element <- NULL
df_dif$Element <- NULL
# df_x <- subset(df_finCoT, Element == "Smart money net position (% of OI)")
# df_x$Element <- NULL
# df_x <- df_x %>% spread(Item, Value)
# df_x <- df_x[, c(1:2)]
# df_y <- df_p %>% spread(Item, Value)
# df_y <- df_y[, c(1:2)]
# start_date <- df_x$date[1]
# ind_start <- which(df_y$date == start_date)
# df_y <- df_y[ind_start:nrow(df_y), ]
# # Change to weekly/monthly
# df_y$date <- as.POSIXct(df_y$date, format = "%Y-%m-%d")
# colnames(df_y)[2] <- "Value"
# df_y <- df_y %>% tq_transmute(select     = Value,
#                mutate_fun = apply.weekly,
#                FUN        = mean,
#                na.rm = T) %>%
#   as.data.frame()
# df_x <- merge(df_x, df_y, by = "date")
#-----------------------------------------------------------------------------
# Get detrended price series
df_dt <- df_p %>% spread(Item, Value)
df_dt[, -1] <- as.data.frame(na.approx(df_dt[, -1]))
indNA <- which(is.na(df_dt$`CC=F`))
if(length(indNA) != 0){df_dt <- df_dt[-indNA, ]}
per_ema_for_detrend <- 21
df_dt[, -1] <- as.data.frame(apply(df_dt[, -1], 2, EMA, per_ema_for_detrend))
df_dt <- df_dt[-c(1:(per_ema_for_detrend - 1)), ]
gathercols <- colnames(df_dt)[-1]
df_dt <- df_dt %>% gather_("Item", "ema", gathercols)
df_dt <- merge(df_dt, df_p, by = c("date", "Item"))
df_dt$Value <- df_dt$Value - df_dt$ema
df_dt$ema <- NULL
#--------------------------------------------------------------------------------
# # Visually inspect
# # Price
# visuallyInspect(df_p, n_cols = 6)
# # High-low difference
# visuallyInspect(df_dif, n_cols = 6)
# # Volume
# visuallyInspect(df_vol, n_cols = 4)
# Detrended price
# visuallyInspect(df_dt, n_cols = 6)
# visuallyInspect(df_finCoT, n_cols = 4)
# visuallyInspect(df_comCoT, n_cols = 4)
#--------------------------------------------------------------------------------
# Set rolling percentile window size
rollWind <- 144 #144 for daily data
#--------------------------------------------------------------------------------
df_pctlP <- df_p %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()

df_pctlDif <- df_dif %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()

df_pctlVol <- df_vol %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()

df_difP <- df_p %>% group_by(Item) %>%
  mutate(Value = c(NA, diff(Value))) %>% as.data.frame()
df_pctlDifP <- df_difP %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, (rollWind - 1), pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
#---
# Set rolling percentile window size for weekly CoT data
rollWind_CoT <- 55
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
# Extract price seasonality series
daily_freq <- 261
weekly_freq <- 52
monthly_freq <- 12
this_freq <- daily_freq
#tsStart <- c(year(df_p$date)[1], week(df_p$date)[1])
#tsStart <- c(year(df_p$date)[1], month(df_p$date)[1])
item_vec <- unique(df_p$Item)
n_items <- length(item_vec)
ratio_vec <- c()
ratioCV_vec <- c()
list_df <- list()
pctDaysBswan_vec <- c()
pctDaysGswan_vec <- c()
for(i in 1:n_items){
  this_item <- item_vec[i]
  print(this_item)
  this_df <- subset(df_p, Item == this_item)
  this_df <- na.trim(this_df)
  this_df$Value <- na.approx(this_df$Value)
  this_dateVec <- this_df$date
  tsStart <- c(year(this_dateVec)[1], yday(this_dateVec)[1])
  this_ts <- ts(this_df$Value, start = tsStart, frequency = this_freq)
  #plot.ts(this_ts)
  ts_decomp <- this_ts %>%
    decompose(type = "multiplicative") #"additive" or "multiplicative"
  #plot(ts_decomp)
  df_out <- data.frame(date = this_dateVec, Item = this_item, Value = as.numeric(ts_decomp$seasonal))
  list_df[[i]] <- df_out
  # mu_ratio <- mean(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T)
  # ratio_vec[i] <- mu_ratio
  # ratioCV_vec[i] <- sd(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T) / mu_ratio
  randFluc <- as.numeric(ts_decomp$random)
  sd_rand <- sd(randFluc, na.rm = T)
  mu_rand <- mean(randFluc, na.rm = T)
  # Num black swan days
  thresh_up <- mu_rand + 2 * sd_rand
  thresh_lo <- mu_rand - 2 * sd_rand
  n_bswan <- length(which(randFluc > thresh_up | randFluc < thresh_lo))
  pct_bswan <- 100 * n_bswan / length(randFluc)
  # Num gray swan days
  thresh_up <- mu_rand + 1 * sd_rand
  thresh_lo <- mu_rand - 1 * sd_rand
  n_gswan <- length(which(randFluc > thresh_up | randFluc < thresh_lo))
  pct_gswan <- 100 * n_gswan / length(randFluc)
  
  pctDaysBswan_vec[i] <- pct_bswan
  pctDaysGswan_vec[i] <- pct_gswan
  #print(autoplot(ts_decomp))
  #Sys.sleep(2)
}
df_randFluc <- data.frame(Item = item_vec, pctDaysBswan = pctDaysBswan_vec, pctDaysGswan = pctDaysGswan_vec)
# hist(pctDaysBswan_vec)
# hist(pctDaysGswan_vec)
# ind_look <- which(pctDaysGswan_vec < 20)
# item_vec[ind_look]
# pctDaysBswan_vec[ind_look]
df_s <- as.data.frame(do.call(rbind, list_df))
#----------------------------------------------------------------------------
# this_series <- subset(df_p, Item == this_guy)$Value
# this_ts <- ts(this_series, start = c(2002, 9), frequency = 52)
# #this_ts <- ts(this_series, start = c(2002, 2), frequency = 12)
# #plot.ts(this_ts)
# ts_decomp <-  this_ts %>% 
#   decompose(type = "multiplicative") #"additive" or "multiplicative"
# #plot(ts_decomp)
# autoplot(ts_decomp)
#-----------------------------------------------------------------------------
# Visually inspect percentile oscillators, seasonality
#visuallyInspect(df_pctlP, n_cols = 6)
#visuallyInspect(df_pctlDif, n_cols = 6)
#visuallyInspect(df_pctlVol, n_cols = 6)
#visuallyInspect(df_pctlDifP, n_cols = 6)
#visuallyInspect(df_s, n_cols = 6)
# visuallyInspect(df_comCoT, n_cols = 5)
# visuallyInspect(df_finCoT, n_cols = 5)
#--------------------------------------------------------------------------------
# Long to wide format
dfWide_p <- df_p %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
dfWide_p <- dfWide_p[-ind_rm, ]
dfWide_p[, -1] <- as.data.frame(na.approx(dfWide_p[, -1]))
numNA_vec <- apply(dfWide_p[, -1], 2, function(x) sum(is.na(x)))
indNA <- which(is.na(dfWide_p$`CC=F`))
if(length(indNA) != 0){dfWide_p <- dfWide_p[-indNA, ]}
#---
dfWide_pctlP <- df_pctlP %>% spread(Item, Value)
dfWide_p <- df_p %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
dfWide_pctlP <- dfWide_pctlP[-ind_rm, ]
#dfWide_pctlP[, -1] <- as.data.frame(na.approx(dfWide_pctlP[, -1]))
# numNA_vec <- apply(dfWide_pctlP[, -1], 2, function(x) sum(is.na(x)))
# indNA <- which(is.na(dfWide_pctlP$`CC=F`))
# if(length(indNA) != 0){dfWide_pctlP <- dfWide_pctlP[-indNA, ]}
#thisPctlP_vec <- dfWide_pctlP[, this_guy]
#
dfWide_pctlDif <- df_pctlDif %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
dfWide_pctlDif <- dfWide_pctlDif[-ind_rm, ]
#dfWide_pctlDif[, -1] <- as.data.frame(na.approx(dfWide_pctlDif[, -1]))
# numNA_vec <- apply(dfWide_pctlDif[, -1], 2, function(x) sum(is.na(x)))
# indNA <- which(is.na(dfWide_pctlDif$`CC=F`))
# if(length(indNA) != 0){dfWide_pctlDif <- dfWide_pctlDif[-indNA, ]}
#
dfWide_pctlVol <- df_pctlVol %>% spread(Item, Value)
ind_rm <- 1:(rollWind - 1)
dfWide_pctlVol <- dfWide_pctlVol[-ind_rm, ]
#dfWide_pctlVol[, -1] <- as.data.frame(na.approx(dfWide_pctlVol[, -1]))
# numNA_vec <- apply(dfWide_pctlVol[, -1], 2, function(x) sum(is.na(x)))
# indNA <- which(is.na(dfWide_pctlVol$`CC=F`))
# if(length(indNA) != 0){dfWide_pctlVol <- dfWide_pctlVol[-indNA, ]}
#
dfWide_difP <- df_difP %>% spread(Item, Value)
#dfWide_pctlDifP[, -1] <- as.data.frame(na.approx(dfWide_pctlDifP[, -1]))
# numNA_vec <- apply(dfWide_pctlDifP[, -1], 2, function(x) sum(is.na(x)))
# indNA <- which(is.na(dfWide_pctlDifP$`CC=F`))
# if(length(indNA) != 0){dfWide_pctlDifP <- dfWide_pctlDifP[-indNA, ]}
dfWide_s <- df_s %>% spread(Item, Value)
# numNA_vec <- apply(dfWide_s[, -1], 2, function(x) sum(is.na(x)))
dfWide_dt <- df_dt %>% spread(Item, Value)
#
dfWide_pctlFinCoT <- df_pctlFinCoT %>% spread(Item, Value)
dfWide_pctlComCoT <- df_pctlComCoT %>% spread(Item, Value)
ind_rm <- 1:(rollWind_CoT - 1)
dfWide_pctlFinCoT <- dfWide_pctlFinCoT[-ind_rm, ]
dfWide_pctlComCoT <- dfWide_pctlComCoT[-ind_rm, ]
#--------------------------------------------------------------------------------
# Select target series
this_guy <- "jpy/aud"#"CT=F"
df_targ <- subset(df_p, Item == this_guy)
df_targ$Item <- NULL
#df_targ <- dfWide_pctlP[, c("date", this_guy)]
#df_targ <- dfWide_dt[, c("date", this_guy)]
colnames(df_targ)[2] <- "targ series"
# Remove leading/trailing NAs if any
df_targ <- na.trim(df_targ)
sum(is.na(df_targ$`targ series`))
#df_targ[which(is.na(df_targ$`targ series`)), ]
# Interpolate remaining NAs
df_targ$`targ series` <- na.approx(df_targ$`targ series`)
# Aggregate to monthly
df_targ$date <- as.POSIXct(df_targ$date, format = "%Y-%m-%d")
colnames(df_targ)[2] <- "targ_series"
df_targ <- df_targ %>% tq_transmute(select = targ_series,
             mutate_fun = apply.monthly,
             FUN        = mean,
             na.rm = T) %>%
  as.data.frame()
df_targ$date <- as.yearmon(df_targ$date)
#--------------------------------------------------------------------------------
# Cross correlation analysis
#https://stats.stackexchange.com/questions/405204/time-series-confused-about-identification-of-possibly-an-armap-q-model
#https://otexts.com/fpp2/seasonal-plots.html
#https://stats.stackexchange.com/questions/79312/why-sinusoid-pattern-in-correlogram
#https://rstudio-pubs-static.s3.amazonaws.com/468354_117d044751f54c389d2f7a741f13a31e.html
#https://rpubs.com/davoodastaraky/TSA1
# myts2 <- window(myts, start=1983)
# gglagplot(myts2)
# Select type of variable in which to look for correlation with target
#df_ccf <- dfWide_dt %>% select(!c(this_guy))
#df_ccf <- dfWide_pctlFinCoT
#df_ccf <- dfWide_pctlComCoT
#df_ccf <- dfWide_pctlP %>% select(!c(this_guy))
#df_ccf <- dfWide_pctlDif
#df_ccf <- dfWide_s
#df_ccf <- dfWide_pctlVol
colnames(dfWide_pctlDif)[-1] <- paste("pctlDif", colnames(dfWide_pctlDif)[-1])
colnames(dfWide_difP)[-1] <- paste("difP", colnames(dfWide_difP)[-1])
colnames(dfWide_pctlP)[-1] <- paste("pctlP", colnames(dfWide_pctlP)[-1])
colnames(dfWide_pctlVol)[-1] <- paste("pctlVol", colnames(dfWide_pctlVol)[-1])
colnames(dfWide_dt)[-1] <- paste("dt", colnames(dfWide_dt)[-1])
list_df <- list(dfWide_pctlDif, dfWide_pctlP, dfWide_pctlVol)
df_ccf <- plyr::join_all(list_df, by = "date")
# Aggregate to monthly
df_ccf$date <- as.POSIXct(df_ccf$date, format = "%Y-%m-%d")
df_ccf <- df_ccf %>% tq_transmute(select = colnames(df_ccf)[-1],
               mutate_fun = apply.monthly,
               FUN        = mean,
               na.rm = T) %>%
  as.data.frame()
df_ccf$date <- as.yearmon(df_ccf$date)
# Merge target series with covariates
df_ccf <- merge(df_ccf, df_targ, by = "date")
#------
# Quick check for NAs and NaNs
numNA_vec <- apply(df_ccf[, -1], 2, function(x) sum(is.na(x)))
table(numNA_vec)
numNaN_vec <- apply(df_ccf[, -1], 2, function(x) sum(is.nan(x)))
table(numNaN_vec)
#----------------------------------------------------------------------------
# Remove highly correlated items
mat_ccf <- as.matrix(df_ccf[, -1])
#cormat <- cor(mat_ccf)
cormat <- sandwich::lrvar(mat_ccf)
cormat <- diag(1 / sqrt(diag(cormat))) %*% cormat %*% diag(1 / sqrt(diag(cormat)))
colnames(cormat) <- colnames(mat_ccf)
row.names(cormat) <- colnames(mat_ccf)
#corrplot::corrplot(cormat)
cormat[lower.tri(cormat)] <- 0
diag(cormat) <- 0
hiCor_vec <- c()#; t <- 0
for(i in 1:ncol(cormat)){
  this_col <- cormat[, i]
  ind_hiCor <- which(this_col > 0.8 | this_col < -0.8)
  if(length(ind_hiCor) != 0){
    hiCor_vec[colnames(cormat)[i]] <- length(ind_hiCor)
  }
}
hiCor_vec
rm_cols <- which(colnames(df_ccf) %in% names(hiCor_vec))
#rm_cols <- which(colnames(df_ccf) %in% c("XLI", "XLF", "ZN=F", "usd/eur"))
#rm_cols <- which(colnames(df_ccf) %in% c("XLI", "XLK", "XLF", "ZN=F"))
#rm_cols <- which(colnames(df_ccf) %in% c("XLI", "XLK", "XLF", "ZN=F", "usd/eur"))
#rm_cols <- which(colnames(df_ccf) %in% c("ES=F", "ZN=F", "usd/eur"))
# rm_cols <- which(colnames(df_ccf) %in% c("S&P 500 Consolidated - CME",
#                                           "E-MINI S&P 500 STOCK INDEX - CME",
#                                           "U.S. DOLLAR INDEX - ICE"))
df_ccf <- df_ccf[, -rm_cols]
mat_ccf <- as.matrix(df_ccf[, -1])
#cormat <- cor(mat_ccf)
cormat <- sandwich::lrvar(mat_ccf)
cormat <- diag(1 / sqrt(diag(cormat))) %*% cormat %*% diag(1 / sqrt(diag(cormat)))
colnames(cormat) <- colnames(mat_ccf)
row.names(cormat) <- colnames(mat_ccf)
corrplot::corrplot(cormat)
#----------------------------------------------------------------------------
# Get correlated lagged series
max_lag <- 25
max_n_lags <- 4
df_feat <- get_correlatedLaggedSeries(df_ccf, max_lag, max_n_lags)
max_lag <- max(as.numeric(gsub(".*-([0-9]+)", "\\1", colnames(df_feat))), na.rm = T)
df_feat$`targ_series 0` <- NULL
targSeries_mod <- df_feat$`targ series`
indRow_rm <- (nrow(df_feat) - max_lag + 1):nrow(df_feat)
df_feat <- df_feat[-indRow_rm, ]
# If doing CoT
# df_feat <- as.data.frame(cbind(date_vec, df_ccf))
# colnames(df_feat)[1] <- "date"
# df_feat <- merge(df_feat, df_targMo, by = "date")
#----------------------------------------------------------------------------
#df_feat <- merge(df_feat, dfWide_s)
#----------------------------------------------------------------------------
df_mod <- df_feat
df_mod$date <- NULL
df_mod <- df_mod[, -grep("chf/usd", colnames(df_mod))]
#df_mod <- log(df_mod)
#colnames(df_mod)[1:3] <- c("x1", "x2", "x3")
#mod <- dynlm::dynlm(`targ series` ~., df_mod)
mod <- lm(`targ series` ~., df_mod)
Anova(mod)
plot(mod$fitted.values, mod$residuals)
# vif(mod)
# vifMod <- vif(mod)
# ind_rm <- which(vifMod == max(vifMod))
# df_mod <- df_mod[, -ind]
# mod <- lm(`targ series` ~., df_mod)
# Anova(mod)
# plot(mod$fitted.values, mod$residuals)

# df_train = df_mod %>%
#   sample_frac(0.6)
# df_test = df_mod %>%
#   setdiff(df_train)
# mod <- lm(`targ series` ~., df_train)
# Anova(mod)
# plot(mod$fitted.values, mod$residuals)
# pred = predict(mod, df_test[, -1])
# resid <- pred - df_test$`targ series`
# mean(resid^2)
# plot(pred, resid)
# plot(df_test$`targ series`, pred)

x <- Anova(mod)$`Pr(>F)`
pvals <- x[-length(x)]
item_vec <- colnames(df_feat)[-c(1, ncol(df_feat))]
signif_items <- item_vec[(which(pvals <= 0.1))]
signif_items <- gsub("`", "", signif_items)
ind_keep <- which(colnames(df_mod) %in% c("targ series", signif_items))
df_mod <- df_mod[, ind_keep]
mod <- lm(`targ series` ~., df_mod)
Anova(mod)
plot(mod$fitted.values, mod$residuals)
#vif(mod)
#car::Anova(mod)
#coeftest(mod, vcov. = sandwich)
# mat_mod <- as.matrix(df_mod[-indRow_rm, -ncol(df_mod)])
# cormat <- cor(mat_mod)
# cormat_inv <- MASS::ginv(cormat)
# colnames(cormat_inv) <- colnames(df_mod)[-ncol(df_mod)]
# row.names(cormat_inv) <- colnames(df_mod)[-ncol(df_mod)]
# corrplot::corrplot(cormat_inv, method = "number", is.corr = F)
# car::linearHypothesis(mod, colnames(df_mod)[1:3])
# PCR 1
df_mod <- df_feat
df_mod$date <- NULL
mod_pcr <- pls::pcr(`targ series` ~., data = df_mod, scale = TRUE, validation = "CV")#, ncomp = 27)
#summary(mod_pcr)
validationplot(mod_pcr, val.type="MSEP")
predplot(mod_pcr)

df_train = df_mod %>%
  sample_frac(0.6)

df_test = df_mod %>%
  setdiff(df_train)

this_ncomp <- 10
mod_pcr = pcr(`targ series` ~., data = df_train, scale = F, validation = "CV", ncomp = this_ncomp)
validationplot(mod_pcr, val.type = "MSEP")


pcr_pred = predict(mod_pcr, df_test[, -ncol(df_test)], ncomp = this_ncomp)
resid <- pcr_pred - df_test$`targ series`
mean(resid^2)
plot(pcr_pred, resid)
plot(df_test$`targ series`, pcr_pred)
#============================================================================
# PCA
mat_X_in <- as.matrix(df_mod)
# o <- apply(mat_X_in, 2, function(x) sum(is.na(x)))
# table(o)
list_out <- get_S_and_corrXS(mat_X_in)
eigVals <- list_out[[3]]
mat_L <- list_out[[2]]
n_sigs <- 25
mat_L <- mat_L[, 1:n_sigs]
#eigVals[1:n_sigs]

fx_vec <- colnames(df_mod)[grep("/", colnames(df_mod))]
pattern_us <- "ES=F|ZB=F|ZN=F|ZF=F|XLF|XLI|XLK|IYR"
us_vec <- colnames(df_mod)[grep(pattern_us, colnames(df_mod))]
pattern_commod <- "GC=F|NG=F|CL=F|CT=F|KC=F|CC=F|SB=F"
commod_vec <- colnames(df_mod)[grep(pattern_commod, colnames(df_mod))]
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
#============================================================================
# PCR 2
df_mod <- df_feat
df_mod$date <- NULL
df_mod$`targ series` <- NULL
mat_X_in <- as.matrix(df_mod)
list_out <- get_S_and_corrXS(mat_X_in)
mat_S <- list_out[[1]]
mat_S <- mat_S[, 1:n_sigs]
df_mod <- as.data.frame(cbind(df_feat$`targ series`, mat_S))
colnames(df_mod) <- c("targ series", paste("Signal", 1:n_sigs))
#mod <- dynlm::dynlm(`targ series` ~., df_mod)
mod <- lm(`targ series` ~., df_mod)
Anova(mod)
plot(mod$fitted.values, mod$residuals)
#vif(mod)

df_mod <- df_feat
df_mod$date <- NULL
df_mod$`targ series` <- NULL
# df_train = df_mod %>%
#   sample_frac(0.6)
# df_test = df_mod %>%
#   setdiff(df_train)

trainPctTot <- 0.6
ind_train <- 1:round(trainPctTot * nrow(df_mod))
ind_test <- (ind_train[length(ind_train)] + 1):(nrow(df_mod) - length(ind_train))
mat_train <- as.matrix(df_mod[ind_train, ])
mat_test <- as.matrix(df_mod[ind_test, ])
list_out <- get_S_and_corrXS(mat_train)
mat_S_train <- list_out[[1]]; mat_S_train <- mat_S_train[, 1:n_sigs]
mat_P <- list_out[[4]]
mat_S_test <- scale(mat_test, scale = F) %*% mat_P
mat_S_test <- mat_S_test[, 1:n_sigs]

df_train <- as.data.frame(cbind(df_feat$`targ series`[ind_train], mat_S_train))
colnames(df_train) <- c("targ series", paste("Signal", 1:n_sigs))
mod <- lm(`targ series` ~., df_train)
Anova(mod)
plot(mod$fitted.values, mod$residuals)
df_test <- as.data.frame(cbind(df_feat$`targ series`[ind_test], mat_S_test))
colnames(df_test) <- c("targ series", paste("Signal", 1:n_sigs))
pcr_pred = predict(mod, df_test[, -1])
resid <- pcr_pred - df_test$`targ series`
mean(resid^2)
plot(pcr_pred, resid)
plot(df_test$`targ series`, pcr_pred)

df_featPCA <- df_mod
df_featPCA$date <- df_feat$date
rearrange_cols <- c("date", colnames(df_featPCA)[-ncol(df_featPCA)])
df_featPCA <- df_featPCA[, rearrange_cols]
#============================================================================
#============================================================================
#============================================================================
# Get trends
#df <- subset(df_p, Item == this_guy)
#df$date <- as.Date(df$date)
# df <- na.trim(df)
# df$Value <- na.approx(df$Value)
# #df$date_chr <- as.character(df$date)
# #df <- df[, c("date", "date_chr", "Value")]
# df$Item <- NULL
# df <- df[, c("date", "Value")]
#colnames(df)[ncol(df)] <- "p"
# df <- merge(df, df_feat, by = "date")
# df <- df[, c("date", "p")]

#df <- df_feat[, c("date", "targ series")]
df <- df_featPCA[, c("date", "targ series")]
ind_critPos <- findpeaks(df$`targ series`)[, 2]
ind_critNeg <- findpeaks(-df$`targ series`)[, 2]
#diff(ind_critPos)
ind_crit <- c(ind_critPos, ind_critNeg)
df_plot <- df
df_plot$Crit <- NA
df_plot$Crit[ind_crit] <- df_plot$`targ series`[ind_crit]
gg <- ggplot(df_plot, aes(x = date, y = `targ series`))
gg <- gg + geom_line()
gg <- gg + geom_point(aes(x = date, y = Crit))
gg


list_out <- getTsTrendInds(df,
                           pctThresh_uptrend = 3,
                           pctThresh_dntrend = -3)
df_up <- list_out[[1]]
df_down <- list_out[[2]]
df_up <- subset(df_up, `False trend` == 0)
df_down <- subset(df_down, `False trend` == 0)
# Consolidate trends
# thresh_timeBetwn <- 31 # Has to be expressed in days, eg. 3 months = 92 days.
# df_down <- consolidateTrnds(df_down, thresh_timeBetwn, show_seqs = F)
# df_up <- consolidateTrnds(df_up, thresh_timeBetwn, show_seqs = F)
# df_down <- subset(df_down, `Pct. Change` <= -2)
# df_up <- subset(df_up, `Pct. Change` >= 2)
#-----------------------------------------------------------------------------
# Visually inspect trend capture
# hist(df_up$`Pct. Change`)
# hist(df_down$`Pct. Change`)
# hist(df_down$`Pct. Change/Time`)
# hist(df_up$`Pct. Change/Time`)
df_plot <- df
df_plot$date <- as.Date(df_plot$date)
#df_plot$date_chr <- as.character(df$date)
df_plotUp <- df_up
df_plotDown <- df_down
df_plotUp$`Start date` <- as.Date(as.yearmon(df_plotUp$`Start date`))
df_plotUp$`Stop date` <- as.Date(as.yearmon(df_plotUp$`Stop date`))
df_plotDown$`Stop date` <- as.Date(as.yearmon(df_plotDown$`Stop date`))
df_plotDown$`Start date` <- as.Date(as.yearmon(df_plotDown$`Start date`))
# df_plotUp$`Start date_chr` <- as.character(df_plotUp$`Start date`)
# df_plotUp$`Stop date_chr` <- as.character(df_plotUp$`Stop date`)
# df_plotDown$`Stop date_chr` <- as.character(df_plotDown$`Stop date`)
# df_plotDown$`Start date_chr` <- as.character(df_plotDown$`Start date`)

# Evenly spaced breaks
my_breaks <- df_plot$date[seq.int(1, length(df_plot$date), length.out = 30)]

gg <- ggplot()
gg <- gg + geom_rect(data = df_plotUp, aes(xmin = `Start date`, xmax = `Stop date`,
                                       ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + geom_rect(data = df_plotDown, aes(xmin = `Start date`, xmax = `Stop date`,
                                         ymin = -Inf, ymax = Inf, fill = `Pct. Change`), alpha = 0.7)
gg <- gg + scale_fill_gradient2(low = "darkmagenta", mid = "khaki2", high = "green")
gg <- gg + geom_line(data = df_plot, aes(x = date, y = `targ series`, group = 1))
gg <- gg + scale_x_date(breaks = my_breaks)
gg <- gg + theme(legend.position = "top",
                 # axis.title.x = element_blank(),
                 axis.text.x = element_text(angle = 60, hjust = 1))#,
                 #legend.title = element_blank())
gg <- gg + guides(fill = guide_colorbar(title.position="top",
                                        title.hjust = 0.5),
                  color = "none")
#gg_p <- gg
gg
#-----------------------------------------------------------------------------
# Get the dependent variable sorted out
df_upStart <- data.frame(date = df_up$`Start date`, yUp = "Buy")
df_upStop <- data.frame(date = df_up$`Stop date`, yUp = "Sell")
df_yUp <- as.data.frame(rbind(df_upStart, df_upStop))
df_yUp$yUp <- as.character(df_yUp$yUp)
#---
df_dnStart <- data.frame(date = df_down$`Start date`, yDn = "Sell")
df_dnStop <- data.frame(date = df_down$`Stop date`, yDn = "Buy")
df_yDn <- as.data.frame(rbind(df_dnStart, df_dnStop))
df_yDn$yDn <- as.character(df_yDn$yDn)
#df_yDn <- df_yDn[order(df_yDn$date), ]
df_y <- merge(df_yUp, df_yDn, by = "date", all = T)
df_y$y <- df_y$yUp
ind <- which(!is.na(df_y$yDn))
df_y$y[ind] <- df_y$yDn[ind]
df_y <- df_y[, c("date", "y")]
df_y$date <- as.yearmon(df_y$date)
df_y <- df_y[order(df_y$date), ]
#-----------------------------------------------------------------------------
# Insert intervening months between action (buy/sell) months
# df_yUp$date <- as.Date(as.yearmon(df_yUp$date))
# df_yUp <- df_yUp[order(df_yUp$date), ]
# intDate_vec <- as.yearmon(seq(df_yUp$date[1], df_yUp$date[nrow(df_yUp)], by = "month"))
# df_intDate <- data.frame(date = intDate_vec)
# df_yUp$date <- as.yearmon(df_yUp$date)
# df_yUp <- merge(df_intDate, df_yUp, by = "date", all = T)
# ind_buy <- which(df_yUp$yUp == "Buy")
# ind_sell <- which(df_yUp$yUp == "Sell")
# for(i in 1:length(ind_buy)){
#   df_yUp$yUp[ind_buy[i]:ind_sell[i]] <- "Up"
# }
# #---
# df_yDn$date <- as.Date(as.yearmon(df_yDn$date))
# df_yDn <- df_yDn[order(df_yDn$date), ]
# intDate_vec <- as.yearmon(seq(df_yDn$date[1], df_yDn$date[nrow(df_yDn)], by = "month"))
# df_intDate <- data.frame(date = intDate_vec)
# df_yDn$date <- as.yearmon(df_yDn$date)
# df_yDn <- merge(df_intDate, df_yDn, by = "date", all = T)
# ind_buy <- which(df_yDn$yDn == "Buy")
# ind_sell <- which(df_yDn$yDn == "Sell")
# for(i in 1:length(ind_sell)){
#   df_yDn$yDn[ind_sell[i]:ind_buy[i]] <- "Dn"
# }
# #---
# df_y <- merge(df_yUp, df_yDn, by = "date", all = T)
# df_y$y <- df_y$yUp
# ind <- which(!is.na(df_y$yDn) & is.na(df_y$yUp))
# df_y$y[ind] <- df_y$yDn[ind]
# df_y$y[which(is.na(df_y$y))] <- "Hold"
# df_y <- df_y[, c("date", "y")]
#-----------------------------------------------------------------------------
#df_mod <- merge(df_y, df_feat, by = "date", all = T)
df_mod <- merge(df_y, df_featPCA, by = "date", all = T)
#df_mod$y[grep("Buy|Sell", df_mod$y)] <- "Act"
#ind_keep <- which(colnames(df_mod) %in% c("date", "y", "targ series", signif_items))
#df_mod <- df_mod[, ind_keep]
# df_mod$`targ series` <- NULL
# rm_col <- paste(this_guy, "0")
# rm_col <- which(colnames(df_mod) == rm_col)
# df_mod <- df_mod[, -rm_col]
# rm_col <- grep("targ_series", colnames(df_mod))
# rm_col <- grep(this_guy, colnames(df_mod))
# rm_col <- grep("/", colnames(df_mod))
# df_mod <- df_mod[, -rm_col]
#---
ind_exclud <- which(colnames(df_mod) %in% c("date", "y"))
numNA_vec <- apply(df_mod[, -ind_exclud], 2, function(x) sum(is.na(x)))
table(numNA_vec)
# df_mod <- na.trim(df_mod)
# df_mod[, -ind_exclud] <- na.approx(df_mod[, -ind_exclud])
df_mod$y[which(is.na(df_mod$y))] <- "Hold"
df_mod$y <- as.factor(df_mod$y)
df_mod$y <- relevel(df_mod$y, ref = "Hold")
contrasts(df_mod$y)
df_mod$date <- NULL

trainDat_pctot <- .7

indtrain_beg <- 1
indtrain_end <- round(nrow(df_mod) * trainDat_pctot)
indtest_beg <- indtrain_end + 1
indtest_end <- nrow(df_mod)
indtrain <- indtrain_beg:indtrain_end
indtest <- indtest_beg:indtest_end
df_train <- df_mod[indtrain, ]
df_test <- df_mod[indtest, ]

df_train = df_mod %>%
  sample_frac(trainDat_pctot)
df_test = df_mod %>%
  setdiff(df_train)



# Multinomial logistic regression using nnet package
#https://medium.com/@PAdhokshaja/using-anova-and-multinomial-logistic-regression-to-predict-human-activity-cd2101a5e8bf

# model <- nnet::multinom(y ~ ., data = df_mod, maxit = 2500)
# Anova(model)

model <- nnet::multinom(y ~ ., data = df_train, maxit = 2500)
Anova(model)
#summary(model)
#library(AER)
#coeftest(model)
pvals <- Anova(model)$`Pr(>Chisq)`
ind_signif <- which(pvals <= 0.1) + 1

cols_reduced <- c(colnames(df_train)[ind_signif], "y")
df_train <- df_train[, cols_reduced]
df_test <- df_test[, cols_reduced]
# # df_broom <- as.data.frame(broom::tidy(model))
# # df_broom[which(df_broom$p.value < 0.05), ]
model <- nnet::multinom(y ~ ., data = df_train, maxit = 2500)
Anova(model)
#Prediction
mat_pred <- predict(model, df_test, type = "class")
unique(mat_pred)
df_pred <- as.data.frame(mat_pred)
colnames(df_pred)[1] <- "yPred"
df_pred$yObs <- df_test$y
#misClasificError <- mean(df_pred$yPred != df_pred$yObs)
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
#===========================================================================
#===========================================================================
#===========================================================================
# normalized gini function taked from:
# https://www.kaggle.com/c/ClaimPredictionChallenge/discussion/703
normalizedGini <- function(aa, pp) {
  Gini <- function(a, p) {
    if (length(a) !=  length(p)) stop("Actual and Predicted need to be equal lengths!")
    temp.df <- data.frame(actual = a, pred = p, range=c(1:length(a)))
    temp.df <- temp.df[order(-temp.df$pred, temp.df$range),]
    population.delta <- 1 / length(a)
    total.losses <- sum(a)
    null.losses <- rep(population.delta, length(a)) # Hopefully is similar to accumulatedPopulationPercentageSum
    accum.losses <- temp.df$actual / total.losses # Hopefully is similar to accumulatedLossPercentageSum
    gini.sum <- cumsum(accum.losses - null.losses) # Not sure if this is having the same effect or not
    sum(gini.sum) / length(a)
  }
  Gini(aa,pp) / Gini(aa,aa)
}

# create the normalized gini summary function to pass into caret
giniSummary <- function (data, lev = "Hold", model = NULL) {
  #levels(data$y) <- c("Hold", "Buy", "Sell")
  out <- normalizedGini(as.numeric(levels(data$y))[data$y], data[, lev[1]])  
  names(out) <- "NormalizedGini"
  out
}
#===========================================================================


trainDat_pctot <- .7
indtrain_beg <- 1
indtrain_end <- round(nrow(df_mod) * trainDat_pctot)
indtest_beg <- indtrain_end + 1
indtest_end <- nrow(df_mod)
indtrain <- indtrain_beg:indtrain_end
indtest <- indtest_beg:indtest_end
df_train <- df_mod[indtrain, ]
df_test <- df_mod[indtest, ]

#this_method <- "glmnet"
#this_method <- "glm"
#this_method <- "nb"
#this_method <- "gbm" #good
#this_method <- "xgbLinear"
this_method <- "xgbTree" #better
#this_method <- "naive_bayes"
#this_method <- "svmRadial"
#this_method <- "cforest"
#this_method <- "rf"




# create the training control object. Two-fold CV to keep the execution time under the kaggle
# limit. You can up this as your compute resources allow. 
trControl = trainControl(method = "repeatedcv", 
                         number = 5, 
                         repeats = 3, 
                         verboseIter = T)

# trControl = trainControl(
#   method = 'cv',
#   number = 2,
#   summaryFunction = giniSummary,
#   classProbs = TRUE,
#   verboseIter = TRUE,
#   allowParallel = TRUE)

# create the tuning grid. Again keeping this small to avoid exceeding kernel memory limits.
# You can expand as your compute resources allow. 
tuneGridXGB <- expand.grid(
  nrounds=c(350),
  max_depth = c(4, 6),
  eta = c(0.05, 0.1),
  gamma = c(0.01),
  colsample_bytree = c(0.75),
  subsample = c(0.50),
  min_child_weight = c(0))

#start <- Sys.time()

# train the xgboost learner
xgbmod <- train(
  x = df_train[, -ncol(df_train)],
  y = df_train$y,
  method = 'xgbTree',
  #metric = 'NormalizedGini',
  trControl = trControl,
  tuneGrid = tuneGridXGB)


#print(Sys.time() - start)
confusionMatrix(data = predict(xgbmod, df_test),
                reference = df_test$y)



# make predictions
preds <- predict(xgbmod, newdata = df_test, type = "prob")
#preds_final <- predict(xgbmod, newdata = dtest, type = "prob")


# convert test target values back to numeric for gini and roc.plot functions
# levels(y_test) <- c("0", "1")
# y_test_raw <- as.numeric(levels(y_test))[y_test]

# Diagnostics
print(xgbmod$results)
print(xgbmod$resample)

# plot results (useful for larger tuning grids)
#plot(xgbmod)

# score the predictions against test data
normalizedGini(df_test$y, preds$Yes)

# plot the ROC curve
#roc.plot(df_test$y, preds$Hold, plot.thres = c(0.02, 0.03, 0.04, 0.05))






# prep the predictions for submissions
# sub <- data.frame(id = as.integer(dtest$id), target = preds_final$Yes)
# write to csv
# write.csv(sub, 'xgb_submission.csv', row.names = FALSE)








































objControl <- trainControl(method = 'repeatedcv', number = 5)#, repeats = 3)
objModel <- train(df_train[, -ncol(df_train)], df_train$y,
                  method = this_method,
                  trControl = objControl)
var_importance <- varImp(objModel, scale = T)
print(var_importance)
predictions <- predict(object = objModel, df_test[, featNames], type = "raw")#type = 'raw') #type='prob')
print(postResample(pred = predictions, obs = as.factor(df_test[, yName])))












# This is a minimal framework for training xgboost in R using caret to do the cross-validation/grid tuning
# and using the normalized gini metric for scoring. The # of CV folds and size of the tuning grid
# are limited to remain under kaggle kernel limits. To improve the score up the nrounds and expand
# the tuning grid.

library(data.table)
library(caret)
library(xgboost)
library(verification)

# Read train and test data
dtrain <- fread('../input/train.csv')
dtest <- fread('../input/test.csv')

# Check data size in memory
print("Training data size in RAM:");
print(object.size(dtrain), units = 'Mb')

# print training data dimensions
print(dim(dtrain))

# collect names of all categorical variables
cat_vars <- names(dtrain)[grepl('_cat$', names(dtrain))]

# turn categorical features into factors
dtrain[, (cat_vars) := lapply(.SD, factor), .SDcols = cat_vars]
dtest[, (cat_vars) := lapply(.SD, factor), .SDcols = cat_vars]

# one hot encode the factor levels
dtrain <- as.data.frame(model.matrix(~. - 1, data = dtrain))
dtest <- as.data.frame(model.matrix(~ . - 1, data = dtest))

# create index for train/test split
train_index <- sample(c(TRUE, FALSE), size = nrow(dtrain), replace = TRUE, prob = c(0.8, 0.2))

# perform x/y ,train/test split.
x_train <- dtrain[train_index, 3:ncol(dtrain)]
y_train <- as.factor(dtrain$target[train_index])

x_test <- dtrain[!train_index, 3:ncol(dtrain)]
y_test <- as.factor(dtrain$target[!train_index])

# Convert target factor levels to 0 = "No" and 1 = "Yes" to avoid this error when predicting class probs:
# https://stackoverflow.com/questions/18402016/error-when-i-try-to-predict-class-probabilities-in-r-caret
levels(y_train) <- c("No", "Yes")
levels(y_test) <- c("No", "Yes")

