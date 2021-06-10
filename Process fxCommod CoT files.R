# Import raw Committment of Traders (CoT) data files downloaded from the Commodity
# Futures Tradicing Commission and process into features for model training.
# Raw files downloaded from:
# https://www.cftc.gov/MarketReports/CommitmentsofTraders/HistoricalCompressed/index.htm
# Advice on what the different long/short positions are:
# https://freecotdata.com/how-to-use/
# "In the commodity markets, speculators are composed of CTAs (commodity
# trading advisors), CPOs (commodity pool operators), hedge funds, other
# reportables, and non-reportables. Other reportables are traders other
# than CTAs, CPOs, and HFs that carry positions above the CFTC's reporting
# limits. Non-reportables are typically called small speculators. These are
# retail traders who own small positions under the reporting limits."

# "In financial futures, speculators are composed of two main trader types:
# asset managers and leveraged money. The asset manager category includes
# mutual funds, endowments, and pension funds. The leveraged money category
# includes the previously mentioned CTAs, CPOs, and hedge funds. In
# addition to these main two categories, other reportables and
# non-reportables are included to form my composite "speculator" category."

#Also: https://www.danielstrading.com/2014/05/15/how-to-read-the-commitment-of-traders-report
#==============================================================================
library(tidyverse)
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
#==============================================================================
yr_vec <- 2010:2021
n_yrs <- length(yr_vec)
#------------------------------------------------------------------------------
# First get financial CoT
this_subfolder <- "financial CoT/"
this_rawDataFolder <- paste0(this_folder, this_subfolder)
list_df <- list()
for(i in 1:n_yrs){
  this_yr <- as.character(yr_vec[i])
  print(this_yr)
  this_file <- paste0("FinFutYY_", this_yr, ".xls")
  this_filepath <- paste0(this_rawDataFolder, this_file)
  this_df <- as.data.frame(readxl::read_xls(this_filepath))
  these_cols <- c("Market_and_Exchange_Names",
                  "Report_Date_as_MM_DD_YYYY",
                  colnames(this_df)[grep("Pct_", colnames(this_df))])
  this_df <- this_df[, these_cols]
  colnames(this_df)[1:2] <- c("Item", "date")
  smart_cols <- c("Pct_of_OI_Dealer_Long_All", "Pct_of_OI_Dealer_Short_All")
  smart_df <- this_df[, c("date", "Item", smart_cols)]
  smart_df$Pct_of_OI_Dealer_Short_All <- -smart_df$Pct_of_OI_Dealer_Short_All
  smart_df <- smart_df %>% gather_("Element", "Value", smart_cols)
  smart_df <- smart_df %>% group_by(date, Item) %>%
    summarize(Value = sum(Value, na.rm = T)) %>% as.data.frame()
  smart_df$Element <- "Smart money net position (% of OI)"
  # this_df$`Smart money net position (% of OI)` <- this_df$Pct_of_OI_Dealer_Long_All -
  #   this_df$Pct_of_OI_Dealer_Short_All
  spec_cols <- c("Pct_of_OI_Asset_Mgr_Long_All",
                 "Pct_of_OI_Lev_Money_Long_All",
                 "Pct_of_OI_Other_Rept_Long_All",
                 "Pct_of_OI_NonRept_Long_All",
                 "Pct_of_OI_Asset_Mgr_Short_All",
                 "Pct_of_OI_Lev_Money_Short_All",
                 "Pct_of_OI_Other_Rept_Short_All",
                 "Pct_of_OI_NonRept_Short_All")
  spec_df <- this_df[, c("date", "Item", spec_cols)]
  short_cols <- colnames(spec_df)[grep("Short", colnames(spec_df))]
  spec_df[, short_cols] <- -spec_df[, short_cols]
  spec_df <- spec_df %>% gather_("Element", "Value", spec_cols)
  spec_df <- spec_df %>% group_by(date, Item) %>%
    summarize(Value = sum(Value, na.rm = T)) %>% as.data.frame()
  spec_df$Element <- "Speculators net position (% of OI)"
  # this_df$`Speculators net position (% of OI)` <- this_df$Pct_of_OI_Asset_Mgr_Long_All +
  #   this_df$Pct_of_OI_Lev_Money_Long_All +
  #   this_df$Pct_of_OI_Other_Rept_Long_All +
  #   this_df$Pct_of_OI_NonRept_Long_All - 
  #   this_df$Pct_of_OI_Other_Rept_Short_All -
  #   this_df$Pct_of_OI_NonRept_Short_All -
  #   this_df$Pct_of_OI_Asset_Mgr_Short_All -
  #   this_df$Pct_of_OI_Asset_Mgr_Short_All
  out_df <- as.data.frame(rbind(smart_df, spec_df))
  out_df <- out_df[, c("date", "Item", "Element", "Value")]
  list_df[[i]] <- out_df
}

df_finCoT <- as.data.frame(do.call(rbind, list_df))
df_finCoT <- df_finCoT[order(df_finCoT$date), ]
# Keep only complete series
df_finCoT <- df_finCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_finCoT$series_length)
df_finCoT <- subset(df_finCoT, series_length == max(u))
df_finCoT$series_length <- NULL
# Shorten names
df_finCoT$Item <- gsub("ICE FUTURES U.S.", "ICE", df_finCoT$Item)
df_finCoT$Item <- gsub("CHICAGO BOARD OF TRADE", "CBoT", df_finCoT$Item)
df_finCoT$Item <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", df_finCoT$Item)
df_finCoT$Item <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", df_finCoT$Item)
df_finCoT$Item <- gsub("COMMODITY EXCHANGE INC.", "CE", df_finCoT$Item)
df_finCoT$Item <- gsub("CBOE FUTURES EXCHANGE", "CBOE", df_finCoT$Item)
# Write file
this_filename <- "Finance CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
write.csv(df_finCoT, this_filepath, row.names = F)
#------------------------------------------------------------------------------
# Now get commodities CoT
this_subfolder <- "commodity CoT/"
this_rawDataFolder <- paste0(this_folder, this_subfolder)
list_df <- list()
for(i in 1:n_yrs){
  this_yr <- as.character(yr_vec[i])
  print(this_yr)
  this_file <- paste0("f_year_", this_yr, ".xls")
  this_filepath <- paste0(this_rawDataFolder, this_file)
  this_df <- as.data.frame(readxl::read_xls(this_filepath))
  these_cols <- c("Market_and_Exchange_Names",
                  "Report_Date_as_MM_DD_YYYY",
                  colnames(this_df)[grep("Pct_", colnames(this_df))])
  this_df <- this_df[, these_cols]
  colnames(this_df)[1:2] <- c("Item", "date")
  smart_cols <- c("Pct_of_OI_Prod_Merc_Long_All", "Pct_of_OI_Prod_Merc_Short_All")
  smart_df <- this_df[, c("date", "Item", smart_cols)]
  smart_df$Pct_of_OI_Prod_Merc_Short_All <- -smart_df$Pct_of_OI_Prod_Merc_Short_All
  smart_df <- smart_df %>% gather_("Element", "Value", smart_cols)
  smart_df <- smart_df %>% group_by(date, Item) %>%
    summarize(Value = sum(Value, na.rm = T)) %>% as.data.frame()
  smart_df$Element <- "Smart money net position (% of OI)"
  # this_df$`Smart money net position (% of OI)` <- this_df$Pct_of_OI_Dealer_Long_All -
  #   this_df$Pct_of_OI_Dealer_Short_All
  spec_cols <- colnames(this_df)[grep("_All", colnames(this_df))]
  spec_cols <- setdiff(spec_cols, smart_cols)
  pattern_rm <- "_Open_Interest_All|_Spread_|_Tot_Rept"
  spec_cols <- spec_cols[-grep(pattern_rm, spec_cols)]
  spec_df <- this_df[, c("date", "Item", spec_cols)]
  short_cols <- colnames(spec_df)[grep("Short", colnames(spec_df))]
  spec_df[, short_cols] <- -spec_df[, short_cols]
  spec_df <- spec_df %>% gather_("Element", "Value", spec_cols)
  spec_df <- spec_df %>% group_by(date, Item) %>%
    summarize(Value = sum(Value, na.rm = T)) %>% as.data.frame()
  spec_df$Element <- "Speculators net position (% of OI)"
  out_df <- as.data.frame(rbind(smart_df, spec_df))
  out_df <- out_df[, c("date", "Item", "Element", "Value")]
  list_df[[i]] <- out_df
}

df_comCoT <- as.data.frame(do.call(rbind, list_df))
df_comCoT <- df_comCoT[order(df_comCoT$date), ]
# Keep only complete series
df_comCoT <- df_comCoT %>% group_by(Item, Element) %>%
  mutate(series_length = length(Value)) %>% as.data.frame()
u <- unique(df_comCoT$series_length)
df_comCoT <- subset(df_comCoT, series_length == max(u))
df_comCoT$series_length <- NULL
# Shorten names
df_comCoT$Item <- gsub("ICE FUTURES U.S.", "ICE", df_comCoT$Item)
df_comCoT$Item <- gsub("CHICAGO BOARD OF TRADE", "CBoT", df_comCoT$Item)
df_comCoT$Item <- gsub("NEW YORK MERCANTILE EXCHANGE", "NYME", df_comCoT$Item)
df_comCoT$Item <- gsub("CHICAGO MERCANTILE EXCHANGE", "CME", df_comCoT$Item)
df_comCoT$Item <- gsub("COMMODITY EXCHANGE INC.", "CE", df_comCoT$Item)
# Write file
this_filename <- "Commodity CoT_processed.csv"
this_filepath <- paste0(this_folder, this_filename)
write.csv(df_comCoT, this_filepath, row.names = F)
