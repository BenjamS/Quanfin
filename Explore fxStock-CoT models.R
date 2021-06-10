#https://online.stat.psu.edu/stat510/lesson/8/8.2
library(tidyverse)
library(tidyquant)
library(patchwork)
library(lubridate)
#--------------------------------------------------------------------------------
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
#=============================================================================
# Define functions
# For percentile oscillator
pctileFun <- function(x){
  out <- ecdf(x)(x[length(x)])
  return(out)
}

#----------------------------------------------------------------------------
# Extract seasonality
getSeasons <- function(df_in,
                       freq = 52,
                       mod_type = "multiplicative",
                       show_graphs = T){
  # daily_freq <- 104
  # weekly_freq <- 52
  # monthly_freq <- 12
  # freq <- daily_freq
  tsStart <- c(year(df_in$date)[1], week(df_in$date)[1])
  #tsStart <- c(year(df_in$date)[1], month(df_in$date)[1])
  item_vec <- unique(df_in$Item)
  n_items <- length(item_vec)
  ratio_vec <- c()
  ratioCV_vec <- c()
  list_df <- list()
  for(i in 1:n_items){
    this_item <- item_vec[i]
    print(this_item)
    this_dateVec <- subset(df_in, Item == this_item)$date
    this_series <- subset(df_in, Item == this_item)$Value
    names(this_series) <- this_dateVec
    this_series <- na.approx(this_series)
    this_series <- na.trim(this_series)
    this_dateVec <- names(this_series)
    this_ts <- ts(this_series, start = tsStart, frequency = freq)
    #plot.ts(this_ts)
    ts_decomp <- this_ts %>%
      decompose(type = mod_type) #"additive" or "multiplicative"
    df_out <- data.frame(date = this_dateVec, Item = this_item, Value = as.numeric(ts_decomp$seasonal))
    list_df[[i]] <- df_out
    mu_ratio <- mean(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T)
    ratio_vec[i] <- mu_ratio
    ratioCV_vec[i] <- sd(abs(ts_decomp$seasonal / ts_decomp$random), na.rm = T) / mu_ratio
    if(show_graphs){
      plot(ts_decomp)
      Sys.sleep(2)
    }

  }
  # names(ratio_vec) <- item_vec
  # names(ratioCV_vec) <- item_vec
  # hist(ratio_vec)
  # hist(ratioCV_vec)
  df_s <- as.data.frame(do.call(rbind, list_df))
  list_out <- list(df_s, ratio_vec, ratioCV_vec)
  return(list_out)
}
#=============================================================================
#=============================================================================
# End function definition
#=============================================================================
#=============================================================================
# Read in the price data
this_filename <- "fxFutData.csv"
this_filepath <- paste0(this_folder, this_filename)
df_raw <- read.csv(this_filepath, stringsAsFactors = F)
df_raw <- subset(df_raw, Item != "XBI")
df_pvhl <- df_raw
df_p <- subset(df_pvhl, Element == "p")
df_vol <- subset(df_pvhl, Element == "Volume")
df_dif <- subset(df_pvhl, Element == "diffHiLo")
df_p$Element <- NULL
df_vol$Element <- NULL
df_dif$Element <- NULL
#-----------------------------------------------------------------------------
# Read in the Commitment of Traders (CoT) data
# Financials CoT data
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
df_in$Value <- na.trim(df_in$Value)
out_list <- getSeasons(df_in,
                       freq = 52,
                       mod_type = "additive",
                       show_graphs = T)
df_s <- out_list[[1]]
ratio1 <- out_list[[2]]
ratio2 <- out_list[[3]]
