library(tidyverse)
library(tidyquant)
library(lubridate)
#=============================================================================
av_api_key("HQBAWJK4Y3YW81VG")
#=============================================================================
fx_vec <- c("eur_gbp", "eur_usd", "aud_usd", "usd_jpy", "aud_jpy", "eur_jpy",
            "usd_cad", "usd_chf")
stock_vec <- c("ES=F", "GC=F", "NG=F", "CL=F", "CT=F", "KC=F", "CC=F",
               "SB=F", "ZB=F", "ZN=F", "ZF=F", "XLF", "XLI", "XLK",
               "IYR", "XLV", "XLY", "XLP")
get_these <- c(fx_vec, stock_vec)
fromdate <- "2002-03-01"
#-----------------------------------------------------------------------------
fx_emaPer_vec <- rep(34, length(fx_vec))
stock_emaPer_vec <- rep(55, length(stock_vec))
emaPer_vec <- c(fx_emaPer_vec, stock_emaPer_vec)
#-----------------------------------------------------------------------------
time_step <- "weekly" #1min, 5min, 15min, 30min, 60min, daily, weekly
#=============================================================================
list_df <- list()
for(i in 1:length(get_these)){
  #---------------------------------------------------------------------
  # AlphaVantage allows only 5 queries per minute.
  # So, have to wait a minute after every 5 queries.
  if(i > 1 & (i - 1) %% 5 == 0){print("waiting a minute..."); Sys.sleep(61)}
  #---------------------------------------------------------------------
  get_this <- get_these[i]
  print(get_this)
  if(length(grep("_", get_this)) != 0){
    fx_pair <- strsplit(get_this, "_")[[1]]
    symb_currency_from <- fx_pair[1]
    symb_currency_to <- fx_pair[2]
    stock_symbol <- NULL
  }else{
    stock_symbol <- get_this
  }
  
  #----------------------------------------------------------------------------
  # Set parameters
  this_emaPer <- emaPer_vec[i]
  per_ema_for_detrend <- this_emaPer
  #----------------------------------------------------------------------------
  time_step_unit <- as.character(stringr::str_extract_all(time_step, "[a-z]+")[[1]])
  #----------------------------------------------------------------------------
  # What is this?
  if(is.null(stock_symbol)){
    if(is.null(symb_currency_from)){
      symbol_domain <- "cryptocurrency"
      
    }else{
      symbol_domain <- "forex"
    }
  }else{
    symbol_domain <- "stocks"
  }
  #----------------------------------------------------------------------------
  if(symbol_domain == "forex"){
    #--------------------------------------------------
    # If forex
    #--------------------------------------------------
    if(time_step_unit == "min"){
      # If intraday
      tbl_ohlcv <- tq_get("", get = "alphavantager", av_fun = "FX_INTRADAY", interval = time_step, from_symbol = symb_currency_from, to_symbol = symb_currency_to, outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
      
    }
    if(time_step_unit == "daily"){
      # If daily
      tbl_ohlcv <- tq_get("", get = "alphavantager", av_fun = "FX_DAILY", from_symbol = symb_currency_from, to_symbol = symb_currency_to, outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
    }
    if(time_step_unit == "weekly"){
      # If weekly
      tbl_ohlcv <- tq_get("", get = "alphavantager", av_fun = "FX_WEEKLY", from_symbol = symb_currency_from, to_symbol = symb_currency_to, outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
    }
    df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
    colnames(df_ohlcv)[1] <- "date"
    this_df <- df_ohlcv[, c("date", "p", "diffHiLo")]
    fx_symbol <- paste0(symb_currency_to, "/", symb_currency_from)
    colnames(this_df)[-1] <- paste(fx_symbol, colnames(this_df)[-1])
    
    
  }
  if(symbol_domain == "stocks"){
    #--------------------------------------------------
    # If stock
    #--------------------------------------------------
    if(time_step_unit == "min"){
      # If intraday
      tbl_ohlcv <- stock_symbol %>%
        tq_get(get = "alphavantager", av_fun = "TIME_SERIES_INTRADAY", interval = time_step, outputsize = "full") 
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, c(3:5)]) / 3
    }
    if(time_step_unit == "daily"){
      # If daily
      #fromdate <- Sys.Date() - starting_how_many_days_ago
      tbl_ohlcv <- tq_get(stock_symbol, get = "stock.prices", from = fromdate)
      # tbl_ohlcv <- stock_symbol %>%
      #   tq_get(get = "alphavantager", av_fun = "TIME_SERIES_DAILY_ADJUSTED", outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- df_ohlcv$adjusted
      
    }
    if(time_step_unit == "weekly"){
      # If weekly
      tbl_ohlcv <- tq_get(stock_symbol, get = "stock.prices", from = fromdate)
      tbl_ohlcv <- stock_symbol %>%
        tq_get(get = "alphavantager", av_fun = "TIME_SERIES_WEEKLY_ADJUSTED", outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- df_ohlcv$adjusted.close
    }
    df_ohlcv$diffHiLo <- df_ohlcv$high - df_ohlcv$low
    colnames(df_ohlcv)[1] <- "date"
    this_df <- df_ohlcv[, c("date", "p", "volume", "diffHiLo")]
    colnames(this_df)[-1] <- paste(stock_symbol, colnames(this_df)[-1])
  }
  
  if(symbol_domain == "cryptocurrency"){
    #--------------------------------------------------
    # If cryptocurrency
    #--------------------------------------------------
    if(time_step_unit == "min"){
      print("Intraday data not available for cryptocurrencies. Try daily or weekly.")
    }
    # If daily
    if(time_step_unit == "daily"){
      tbl_ohlcv <- crypto_symbol_from %>% tq_get(get = "alphavantager", av_fun = "DIGITAL_CURRENCY_DAILY", market = crypto_symbol_to)
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, 2:4]) / 3
    }
    if(time_step_unit == "weekly"){
      # If weekly
      tbl_ohlcv <- crypto_symbol_from %>% tq_get(get = "alphavantager", av_fun = "DIGITAL_CURRENCY_WEEKLY", market = crypto_symbol_to, outputsize = "full")
      df_ohlcv <- as.data.frame(tbl_ohlcv)
      df_ohlcv$p <- rowSums(df_ohlcv[, 2:4]) / 3
    }
    
  }
  
  
  list_df[[i]] <- this_df
  
}

df <- plyr::join_all(list_df, by = "date")

cols_vol <- grep("volume", colnames(df))
cols_dif <- grep("diffHiLo", colnames(df))
cols_p <- setdiff(colnames(df), colnames(df)[c(cols_vol, cols_dif)])

df_p <- df[, cols_p]
df_vol <- df[, c(1, cols_vol)]
df_dif <- df[, c(1, cols_dif)]

gathercols <- colnames(df_vol)[-1]
df_vol <- df_vol %>% gather_("Item", "Value", gathercols)
df_vol$Element <- "Volume"
df_vol$Item <- gsub(" volume", "", df_vol$Item)

gathercols <- colnames(df_dif)[-1]
df_dif <- df_dif %>% gather_("Item", "Value", gathercols)
df_dif$Element <- "diffHiLo"
df_dif$Item <- gsub(" diffHiLo", "", df_dif$Item)

gathercols <- colnames(df_p)[-1]
df_p <- df_p %>% gather_("Item", "Value", gathercols)
df_p$Element <- "p"
df_p$Item <- gsub(" p", "", df_p$Item)

df_out <- as.data.frame(do.call(rbind, list(df_p, df_vol, df_dif)))
#============================================================================
this_folder <- "C:/Users/bensc/OneDrive/Documents/Data/Trading/"
this_filename <- "fxFutData_weekly.csv"
this_filepath <- paste0(this_folder, this_filename)
write.csv(df_out, this_filepath, row.names = F)
