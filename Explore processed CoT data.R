library(tidyverse)
library(tidyquant)
#===============================================================================
# Define functions
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
#-----------------------------------------------------------------------------
# Visual inspection function
visuallyInspect <- function(df_plot, n_cols = 5){
  df_plot$date_chr <- as.character(df_plot$date)
  my_breaks <- df_plot$date_chr[seq.int(1, length(df_plot$date_chr), length.out = 5)]
  gg <- ggplot(df_plot, aes(x = date_chr, y = Value, group = Element, color = Element))
  gg <- gg + geom_line()
  gg <- gg + scale_x_discrete(breaks = my_breaks)
  gg <- gg + facet_wrap(~Item, ncol = n_cols, scales = "free_y")
  gg <- gg + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                   legend.position = "top",
                   legend.title = element_blank())
  print(gg)
  
}
#-----------------------------------------------------------------------------
# For percentile oscillator
pctileFun <- function(x){
  out <- ecdf(x)(x[length(x)])
  return(out)
}
#=============================================================================
#=============================================================================
# End definition of functions
#=============================================================================
#=============================================================================
# Import the processed Commitment of Traders (CoT) data
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
# Financials CoT data
this_filename <- "Finance CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_finCoT <- read.csv(this_filepath, stringsAsFactors = F)
# Commodities CoT data
this_filename <- "Commodity CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
df_comCoT <- read.csv(this_filepath, stringsAsFactors = F)
#=============================================================================
# Explore the commodity CoT data
# ICE Futures, Commodity Exchange Inc., Chicago Board of Trade,
# New York Mercantile Exchange, Chicago Mercantile Exchange
df_check <- df_comCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_check$series_length)
max(u)
length(which(u < max(u)))
df_check <- subset(df_check, series_length == max(u))
visuallyInspect(df_check, n_cols = 5)
# Explore the finance CoT data
# Chicago mercantile exchange, CBOE Futures, ICE Futures
df_check <- df_finCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_check$series_length)
max(u)
length(which(u < max(u)))
df_check <- subset(df_check, series_length == max(u))
visuallyInspect(df_check, n_cols = 5)
#-----------------------------------------------------------------------------
# Get percentile oscillator series
# Set rolling percentile window size
rollWind <- 52
# Commodity oscillator
df_pctlComCoT <- df_comCoT %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
# Financial oscillator
df_pctlFinCoT <- df_finCoT %>% group_by(Item) %>%
  mutate(Value = rollapply(Value, rollWind, pctileFun, fill = NA, align = "right")) %>%
  as.data.frame()
#-----------------------------------------------------------------------------
# PCA Commodity CoT
df_pca <- df_pctlComCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_pca$series_length)
df_pca <- subset(df_pca, series_length == max(u) &
                   Element == "Smart money net position (% of OI)")
df_pca$series_length <- NULL
df_pca$Element <- NULL
df_pca <- df_pca %>% spread(Item, Value)
# Remove leading NAs if using percentile oscillator
ind_rm <- 1:(rollWind - 1)
df_pca <- df_pca[-ind_rm, ]
# Shorten names
colnames(df_pca) <- gsub("ICE FUTURES U.S.", "ICE", colnames(df_pca))
colnames(df_pca) <- gsub("CHICAGO BOARD OF TRADE", "CBoT", colnames(df_pca))
colnames(df_pca) <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", colnames(df_pca))
colnames(df_pca) <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", colnames(df_pca))
colnames(df_pca) <- gsub("COMMODITY EXCHANGE INC.", "CE", colnames(df_pca))
mat_X_in <- as.matrix(df_pca[, -1])
row.names(mat_X_in) <- df_pca$date
out <- get_S_and_corrXS(mat_X_in)
mat_L <- out[[2]]
mat_L <- mat_L[, 1:4]
plot_corrXS_barchart(mat_L, group_info = NULL, xAxis_title = NULL, sigNames = NULL)
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
mat_Lrot <- mat_Lrot[, 1:8]
plot_corrXS_barchart(mat_Lrot, group_info = NULL, xAxis_title, sigNames = NULL)
#-----------------------------------------------------------------------------
# PCA Finance CoT
df_pca <- df_pctlFinCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_pca$series_length)
df_pca <- subset(df_pca, series_length == max(u) &
                   Element == "Smart money net position (% of OI)")
df_pca$series_length <- NULL
df_pca$Element <- NULL
df_pca <- df_pca %>% spread(Item, Value)
# Remove leading NAs if using percentile oscillator
ind_rm <- 1:(rollWind - 1)
df_pca <- df_pca[-ind_rm, ]
# Shorten names
colnames(df_pca) <- gsub("ICE FUTURES U.S.", "ICE", colnames(df_pca))
colnames(df_pca) <- gsub("CHICAGO BOARD OF TRADE", "CBoT", colnames(df_pca))
colnames(df_pca) <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", colnames(df_pca))
colnames(df_pca) <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", colnames(df_pca))
colnames(df_pca) <- gsub("COMMODITY EXCHANGE INC.", "CE", colnames(df_pca))
colnames(df_pca) <- gsub("CBOE FUTURES EXCHANGE", "CBOE", colnames(df_pca))
mat_X_in <- as.matrix(df_pca[, -1])
row.names(mat_X_in) <- df_pca$date
out <- get_S_and_corrXS(mat_X_in)
mat_L <- out[[2]]
mat_L <- mat_L[, 1:4]
plot_corrXS_barchart(mat_L, group_info = NULL, xAxis_title = NULL, sigNames = NULL)
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
mat_Lrot <- mat_Lrot[, 1:8]
plot_corrXS_barchart(mat_Lrot, group_info = NULL, xAxis_title, sigNames = NULL)
