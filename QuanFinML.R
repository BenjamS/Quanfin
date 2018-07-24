#setwd("D:/OneDrive - CGIAR/Documents")
source('./getTsTrends.R', echo=TRUE)
source('./getGlobTrnds.R', echo=TRUE)
source('./tradeSim.R', echo=TRUE)
source('./collectiveModes.R', echo=TRUE)
library(plyr)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(GGally)
library(cluster)
library(mclust)
library(randomForest)
library(corrplot)
library(quantmod)
library(lubridate)
library(xgboost)
library(caret)
library(scales)

fromdate<-"2012-01-01"; todate <- "2018-07-10"
getGlobTrnds(fromdate, todate)
rmcol <- which(colnames(cpgtetfmat) %in% c("EXS1.DE", "^IRX", "EWH", "AIA", "XTN", "NLR", "VXX", "PGD", "MIDD"))
xts_cp_mat <- cpgtetfmat[, -rmcol]
xts_vol_mat <- volgtetfmat[, -rmcol]
namevec <- colnames(xts_cp_mat)
n_names <- length(namevec)

o <- apply(xts_cp_mat, 2, function(x) length(which(is.na(x))))
table(o)
rmcols <- which(o > 60)
colnames(xts_cp_mat)[rmcols]
xts_cp_mat <- xts_cp_mat[, -rmcols]
xts_vol_mat <- xts_vol_mat[, -rmcols]
o <- apply(xts_cp_mat, 1, function(x) length(which(is.na(x))))
table(o)
xts_cp_mat <- na.spline(xts_cp_mat)
xts_vol_mat <- na.spline(xts_vol_mat)
o <- apply(xts_cp_mat, 1, function(x) length(which(is.na(x))))
table(o)
#---------------
n_ts <- ncol(xts_cp_mat)
#---------------

df_cp_mat <- fortify(xts_cp_mat)
colnames(df_cp_mat)[1] <- "Date"
df_zcp_mat <- as.data.frame(scale(df_cp_mat[, 2:ncol(df_cp_mat)]))
df_zcp_mat$Date <- as.character(df_cp_mat$Date)
df_zcp_mat <- df_zcp_mat[, c(ncol(df_zcp_mat), 1:(ncol(df_zcp_mat) - 1))]
df_ts <- df_zcp_mat[-c(1:500),]
#--
GroupInfo_raw <- read.csv("GlobTrndsID.csv", stringsAsFactors = F)
GroupInfo_raw$X <- NULL
colnames(GroupInfo_raw)[1] <- "name"
kept_items <- colnames(df_ts)[2:ncol(df_ts)]
ind_keep <- which(GroupInfo_raw$name %in% kept_items)
df_group <- GroupInfo_raw[ind_keep, c("name", "Sub.type")]
#--
date_vec <- df_ts$Date
df_ts$Date <- NULL
mat_diff <- diff(as.matrix(df_ts))
date_vec <- date_vec[-1]
nrow(mat_diff) / ncol(mat_diff)
out_cm <- collectiveModes(mat_diff, date_vec, df_group,
                          Contrib_as_ModeSq = F,
                          AggregateContributions = F,
                          plot_eigenportfolios = F)





in_ts <- xts_cp_mat[, "NIB"]
out <- tradeSim(in_ts, 
                per_ema = 3, 
                per_slope = 3,
                thresh_pct_uptrend = 1.5,
                thresh_pct_dntrend = -1.5,
                Commission = 7,
                invest_t0 = 500,
                quietly = F)

cat("This ts: ", colnames(in_ts), "\n",
    "NetGain_up_naive: ", out[1], "\n",
    "NetGain_up_shrewd: ", out[2]
)








# #Grid search for best per_ema + per_slope combo
# in_ts <- xts_cp_mat[, "WTI"]
# per_ema_vec <- c(3, 5, 8, 13, 21, 34, 55, 89, 144)
# per_slope_vec <- c(3, 5, 8, 13, 21, 34, 55, 89, 144)
# grid_naive <- matrix(NA, nrow = length(per_ema_vec), ncol = length(per_slope_vec))
# grid_shrewd <- matrix(NA, nrow = length(per_ema_vec), ncol = length(per_slope_vec))
# for(i in 1:length(per_ema_vec)){
#   print(i)
#   for(j in 1:length(per_slope_vec)){
#     print(j)
#     out <- tradeSim(in_ts, 
#                        per_ema = per_ema_vec[i], 
#                        per_slope = per_slope_vec[j],
#                        thresh_pct_uptrend = 1.5,
#                        thresh_pct_dntrend = -1.5,
#                        Commission = 7,
#                        invest_t0 = 500,
#                        quietly = T)
#     grid_naive[i, j] <- out[1]
#     grid_shrewd[i, j] <- out[2]
#     
#   }
# }
# df_naive <- as.data.frame(grid_naive)
# colnames(df_naive) <- paste("slope", per_slope_vec)
#==============================
# Search for ts with highest return
NetGain_up_naive_vec <- c()
NetGain_up_shrewd_vec <- c()
for(i in 1:ncol(xts_cp_mat)){
  in_ts <- xts_cp_mat[, i]
  out <- tradeSim(in_ts, 
                  per_ema = 3, 
                  per_slope = 3,
                  thresh_pct_uptrend = 1.5,
                  thresh_pct_dntrend = -1.5,
                  Commission = 7,
                  invest_t0 = 500,
                  quietly = F)
  
  cat("This ts: ", colnames(in_ts), "\n",
    "NetGain_up_naive: ", out[1], "\n",
    "NetGain_up_shrewd: ", out[2]
    )
  NetGain_up_naive_vec[i] <- out[1] 
  NetGain_up_shrewd_vec[i] <- out[2]
  
}

df_x <- data.frame(name = colnames(xts_cp_mat), 
                   Naive.Net.Gain = NetGain_up_naive_vec,
                   Shrewd.Net.Gain = NetGain_up_shrewd_vec)
gg <- ggplot(df_x, aes(x = Naive.Net.Gain, y = Shrewd.Net.Gain)) + geom_point()
gg <- gg + scale_x_log10() + scale_y_log10()
gg <- gg + geom_text(aes(label = name))
gg


# Take a closer look
in_ts <- xts_cp_mat[, "WTI"]
out <- tradeSim(in_ts, 
                per_ema = 3, 
                per_slope = 3,
                thresh_pct_uptrend = 1.5,
                thresh_pct_dntrend = -1.5,
                Commission = 7,
                invest_t0 = 500,
                quietly = F)

cat("This ts: ", colnames(in_ts), "\n",
    "Naive net gain: ", out[1], "\n",
    "Shrewd net gain: ", out[2], "\n",
    "Number naive trades: ", out[3] - out[4], "\n",
    "Number shrewd trades: ", out[4], "\n",
    "Total number of trades: ", out[3]
)



























datevec <- index(xts_cp_mat)
xts_cp_matEMA <- xts(apply(xts_cp_mat, MARGIN = 2, FUN = "EMA", n = per_ema), datevec)
#---------------
list_df <- list()
for(i in 1:n_ts)
{
  this_sec <- colnames(xts_cp_mat)[i]
  out_df <- fortify(VWMA(xts_cp_mat[, i], xts_vol_mat[, i]))
  colnames(out_df)[2] <- this_sec
  list_df[[i]] <- out_df
}
df_vwma <- join_all(list_df)
datevec <- index(xts_cp_mat)
xts_cp_matVWMA <- xts(df_vwma[, -1], datevec)
# rm_rows <- which(is.na(df_vwma[, 2]))
# df_vwma <- df_vwma[-rm_rows, ]
# datevec <- df_vwma$Index
#====================================
#Get mean, cv of major groups
df_CV <- as.data.frame(t(df_vwma[, 2:ncol(df_vwma)]))
colnames(df_CV) <- as.character(df_vwma$Index)
df_CV$name <- rownames(df_CV)
GroupInfo_raw <- read.csv("GlobTrndsID.csv", stringsAsFactors = F)
GroupInfo_raw$X <- NULL
colnames(GroupInfo_raw)[1] <- "name"
#Excluding "groups" with 2 or less ts
GroupInfo_this <- GroupInfo_raw
rm_these <- names(table(GroupInfo_this$Sub.type)[which(table(GroupInfo_this$Sub.type) <= 2)])
rm_rows <- which(GroupInfo_this$Sub.type %in% rm_these)
GroupInfo_this <- GroupInfo_this[-rm_rows, ]
df_CV <- merge(df_CV, GroupInfo_this[, c("name", "Sub.type")], by = "name")
df_CV$name <- NULL
gathercols <- colnames(df_CV)
df_CVmu <- df_CV %>% group_by(Sub.type) %>% summarise_all(mean)
df_CVsd <- df_CV %>% group_by(Sub.type) %>% summarise_all(sd)
df_groupCV <- as.data.frame(t(df_CVsd[, -1] / df_CVmu[, -1]))
df_groupMU <- as.data.frame(t(df_CVmu[, -1]))
colnames(df_groupCV) <- paste(df_CVmu$Sub.type, "CV")
colnames(df_groupMU) <- paste(df_CVmu$Sub.type, "MU")
o <- apply(df_groupCV, 2, function(x) length(which(is.nan(x))))
table(o)
rmcols <- which(o > 10)
colnames(df_groupCV)[rmcols]
df_groupCV <- df_groupCV[, -rmcols]
df_groupMU <- df_groupMU[, -rmcols]
df_group <- cbind(df_groupCV, df_groupMU)
xts_group <- xts(df_group, as.Date(rownames(df_group)))
#====================================
#Covariance matrix, collective modes
xts_collectiveModes <- xts_cpVWMA
o <- apply(xts_collectiveModes, 1, function(x) length(which(is.na(x))))
table(o)
rmrows <- as.integer(which(o > 0))
xts_collectiveModes_diff <- diff(log(xts_collectiveModes[-rmrows, ]))
xts_collectiveModes_diff <- xts_collectiveModes_diff[-1, ]
o <- apply(xts_collectiveModes_diff, 2, function(x) length(which(is.nan(x))))
table(o)
rmcols <- as.integer(which(o > 0))
xts_collectiveModes_diff <- xts_collectiveModes_diff[, -rmcols]
# which(is.nan(xts_collectiveModes_diff[, "INR"]))
# which(is.nan(xts_collectiveModes_diff[, "BAL"]))
#nrow(xts_collectiveModes)
cormat <- cor(xts_collectiveModes_diff[-1, ])
image(cormat)
e <- eigen(cormat)
#class(e$values)
df_e <- data.frame(values = e$values)
N_t <- nrow(xts_collectiveModes_diff)
N_c <- ncol(xts_collectiveModes_diff)
Q <- N_t / N_c
#s_sq <- 1 - max(df_e$values) / N_c
s_sq <- 1
lam_rand_max <- s_sq * (1 + 1 / Q + 2 / sqrt(Q))
lam_rand_min <- s_sq * (1 + 1 / Q - 2 / sqrt(Q))
lam <- seq(0, lam_rand_max * 1.1, 0.01)
dens_rand <- Q / (2 * pi * s_sq) * sqrt((lam_rand_max - lam) * (lam - lam_rand_min)) / lam
gg <- ggplot() + geom_density(data = df_e, aes(x = values)) + coord_cartesian(xlim = c(0, 4))
gg <- gg + geom_line(data = data.frame(x = lam, y = dens_rand), aes(x = x, y = y))
gg
ind_deviating_from_noise <- which(df_e$values > lam_rand_max)
CollectiveModes <- as.matrix(e$vectors[, ind_deviating_from_noise])
df_collectiveModes <- as.data.frame(CollectiveModes)
n_collectiveModes <- ncol(CollectiveModes)
#colnames(xts_collectiveModes)
names_all_ts <- colnames(xts_collectiveModes_diff)
n_all_ts <- length(names_all_ts)
GroupInfo_clean <- subset(GroupInfo_raw, name %in% names_all_ts)
nameSubtypes <- GroupInfo_clean$Sub.type
ArbGroupTypes <- unique(nameSubtypes)
n_groups <- length(ArbGroupTypes)
n_in_group <- c()
Pvec_list <- list()
for(i in 1:n_groups){
  ind <- which(nameSubtypes == ArbGroupTypes[i])
  ts_in_group <- GroupInfo_clean$name[ind]
  n_in_group <- length(ts_in_group)
  Pvec <- rep(0, n_all_ts)
  Pvec[ind] <- 1 / n_in_group
  Pvec_list[[i]] <- Pvec
}
Pmat <- as.matrix(do.call(cbind, Pvec_list))
nrow(Pmat)
nrow(df_collectiveModes)
Xkl <- t(Pmat) %*% CollectiveModes^2
#barplot(Xkl[, 9])
df_contrib <- as.data.frame(Xkl)
colnames(df_contrib) <- paste("lambda", c(1:n_collectiveModes))
df_contrib$ArbGroup <- ArbGroupTypes
gathercols <- colnames(df_contrib)[1:n_collectiveModes]
df_contrib <- gather_(df_contrib, "Lambda", "Value", gathercols)
gg <- ggplot(df_contrib, aes(x = ArbGroup, y = Value)) + geom_bar(stat="identity")
gg <- gg + facet_wrap(~ Lambda, ncol = 3) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
gg
#--
PorterTom <- rnorm(n_all_ts, 0, 1)
colnames(df_collectiveModes) <- paste("lambda", c(1:n_collectiveModes))
df_collectiveModes$PorterTom <- PorterTom
gathercols <- colnames(df_collectiveModes)
df_collectiveModes <- gather_(df_collectiveModes, "Lambda", "Value", gathercols)
gg <- ggplot(df_collectiveModes, aes(x = Value)) + geom_density()
gg <- gg + facet_wrap(~ Lambda, ncol = 3)
gg

#====================================
#Experiment with tstrends()
xts_inMat <- cbind(xts_cpVWMA, xts_group)
this_one <- "NIB"
in_ts <- xts_inMat[, this_one]
print(colnames(in_ts))
out <- tsTrends(in_ts, quietly = F)
#====================================
#xts_inMat <- cbind(xts_cpVWMA, xts_group)
xts_inMat <- xts_cpVWMA
list_df_dydxmu <- list()
list_df_ldydxcv <- list()
list_df_sigs <- list()
n_ts <- ncol(xts_inMat)
for(i in 1:n_ts){
  in_ts <- xts_inMat[, i]
  this_ts <- colnames(in_ts)
  print(this_ts)
  out <- tsTrends(in_ts, before_window = 0, aft_window = 0, quietly = T)
  df_dydxmu_this <- out[[1]][c("Index", "dydxmu")]
  df_ldydxcv_this <- out[[1]][c("Index", "ldydxcv")]
  df_sigs_this <- out[[1]][c("Index", "SigUpWin")]
  df_sigs_this[, 2] <- as.factor(df_sigs_this[, 2])
  df_sigs_this <- as.data.frame(model.matrix(~.-1, data = df_sigs_this))
  df_sigs_this$Index <- df_dydxmu_this$Index
  colnames(df_dydxmu_this)[2] <- paste(this_ts, colnames(df_dydxmu_this)[2])
  colnames(df_ldydxcv_this)[2] <- paste(this_ts, colnames(df_ldydxcv_this)[2])
  colnames(df_sigs_this) <-   gsub("SigUpWin", "", colnames(df_sigs_this))
  n_col_this <- ncol(df_sigs_this)
  colnames(df_sigs_this)[2:n_col_this] <- paste(this_ts, colnames(df_sigs_this)[2:n_col_this])
  list_df_dydxmu[[i]] <- df_dydxmu_this
  list_df_ldydxcv[[i]] <- df_ldydxcv_this
  list_df_sigs[[i]] <- df_sigs_this
}
df_dydxmu <- join_all(list_df_dydxmu, by = "Index")
df_ldydxcv <- join_all(list_df_ldydxcv, by = "Index")
df_sigs <- join_all(list_df_sigs, by = "Index")
# colnames(df_dydxmu) <- c("Index", paste(colnames(xts_inMat), "dydxmu"))
# colnames(df_ldydxcv) <- c("Index", paste(colnames(xts_inMat), "ldydxcv"))
# colnames(df_sigs) <- c("Index", paste(colnames(xts_inMat), "sig"))
#====================================
#====================================
#====================================
#====================================
#Which time series are you trying to model
name_y <- "SPY"
#====================================
#Get features
# df_feat1 <- df_ldydxcv
# df_feat2 <- df_dydxmu
# n <- ncol(df_feat1)
# colnames(df_feat1)[1:(n - 1)] <- paste(colnames(df_feat1)[1:(n - 1)], "A")
#First remove signals of predicted var
rm_cols <- grep(name_y, colnames(df_sigs))
df_sigs <- df_sigs[, -rm_cols]
#Now assemble features matrix
df_groupCV$Index <- as.Date(rownames(df_groupCV))
df_groupMU$Index <- as.Date(rownames(df_groupMU))
# df_feat <- join_all(list(df_ldydxcv, df_dydxmu, df_groupCV, df_groupMU, df_sigs), by = "Index")
# df_feat <- join_all(list(df_ldydxcv, df_dydxmu, df_groupCV, df_groupMU), by = "Index")
#df_feat <- join_all(list(df_ldydxcv, df_dydxmu, df_groupCV), by = "Index")
df_feat <- join_all(list(df_ldydxcv, df_dydxmu, df_groupCV, df_sigs), by = "Index")
#Get predicted var
in_ts <- xts_cpEMA[, name_y]
print(colnames(in_ts))
before_window <- 5
aft_window <- 2
out_y <- tsTrends(in_ts, quietly = F, before_window = before_window, aft_window = aft_window)
#Trim ts of predicted var to first uptrend start point and last downtrend stop point
#so that ML training not messed up
ind_tsStart <- as.numeric(out_y[[4]][1])
ind_tsFinish <- as.numeric(out_y[[4]][2])
df_y <- out_y[[1]][ind_tsStart:ind_tsFinish, c("Index", "SigUpWin")]
colnames(df_y)[2] <- "y"
n_ts_trimd <- nrow(df_y)
print(n_ts_trimd)
df_mod <- join_all(list(df_y, df_feat), by = "Index")
rownames(df_mod) <- df_mod$Index
df_mod$Index <- NULL
#df_mod$y <- as.factor(df_mod$y)
#Check for NAs and NaNs
length(which(is.na(df_mod$y) == T))
length(which(is.nan(df_mod$y) == T))
o <- apply(df_mod, 2, function(x) length(which(is.na(x))))
table(o)
o <- apply(df_mod, 2, function(x) length(which(is.nan(x))))
table(o)
#df_mod <- as.data.frame(PCA(df_feat, ncp = 14)$ind$coord)

predictorsNames <- colnames(df_mod)[2:ncol(df_mod)]
outcomeName <- colnames(df_mod)[1]


#====================================
# rm_rows <- which(as.character(df_mod$y) == "Hold")
# df_mod <- df_mod[-rm_rows, ]
# df_mod$y <- as.factor(as.character(df_mod$y))
# unique(df_mod$y)
ind <- which(df_mod$y == "Hold")
df_mod <- df_mod[-ind, ]
ind <- which(df_mod$y %in% c("Uptrend False Start", "Uptrend False Stop"))
df_mod$y[ind] <- "Hold"
df_mod$y <- as.factor(df_mod$y)
unique(df_mod$y)
fitDat_pctot <- .95
indfit_beg <- 1
indfit_end <- round(nrow(df_mod) * fitDat_pctot)
indvalid_beg <- indfit_end + 1
indvalid_end <- nrow(df_mod)
fit_rows <- indfit_beg:indfit_end
valid_rows <- indvalid_beg:indvalid_end
df_fit <- df_mod[fit_rows, ]
df_valid <- df_mod[valid_rows, ]
#====================================


# set.seed(1234)
# splitIndex <- createDataPartition(df_fit[, outcomeName], p = .5, list = FALSE, times = 1)
# df_train <- df_fit[ splitIndex,]
# df_test <- df_fit[-splitIndex,]

trainDat_pctot <- .5
indtrain_beg <- 1
indtrain_end <- round(nrow(df_mod) * trainDat_pctot)
indtest_beg <- indtrain_end + 1
indtest_end <- nrow(df_mod)
train_rows <- indtrain_beg:indtrain_end
test_rows <- indtest_beg:indtest_end
df_train <- df_mod[train_rows, ]
df_test <- df_mod[test_rows, ]


nrow(df_mod)
nrow(df_train)
nrow(df_test)
nrow(df_valid)
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

# gbmGrid <-  expand.grid(n.trees = 50, interaction.depth =  c(1, 5, 9),
#                         shrinkage = 0.01, n.minobsinnode = )
# run model
# head(trainDF[,predictorsNames])
# head(trainDF[,outcomeName])
# head(trainDF)
# params <- list(booster = "gbtree", objective = "multi:softmax", eta = 0.1,
#                gamma = 0, max_depth = 25, min_child_weight = 1, subsample = 0.5,
#                colsample_bytree = 0.5, num_class = 5)
#=============================
#modelLookup(this_method)
#=============================
#number = 5 is better (and slower)
objControl <- trainControl(method = 'repeatedcv', number = 3)#, repeats = 3)
# objControl <- trainControl(
#   method = "repeatedcv",
#   number = 5,
#   repeats = 2,
#   returnData = FALSE,
#   classProbs = TRUE,
#   summaryFunction = multiClassSummary
# )

# df_train$y <- as.factor(make.names(as.character(df_train$y)))
# class(df_train$y)
# df_test$y <- as.factor(make.names(as.character(df_test$y)))

# tune_grid <- expand.grid(nrounds=c(25, 40, 60),
#                         max_depth = c(25),
#                         eta = c(0.05, 0.1, 0.3),
#                         gamma = c(0),
#                         colsample_bytree = c(0.5, 0.75),
#                         subsample = c(0.50),
#                         min_child_weight = c(0))

objModel <- train(df_train[, predictorsNames], df_train[, outcomeName],
                  method = this_method,
                  trControl = objControl)
#-----------------------------
var_importance <- varImp(objModel, scale = T)
print(var_importance)
#plot(var_importance)
#-----------------------------
#Predict using fitted model
# probabilities ("prob") or integer ("raw")
predictions <- predict(object = objModel, df_test[, predictorsNames], type = 'raw') #type='prob')
print(postResample(pred = predictions, obs = as.factor(df_test[,outcomeName])))
#head(predictions)
#df_plot <- rownames(df_test)
#----------------------------
#Confusion matrix
x <- confusionMatrix(predictions, df_test[, outcomeName])
confmat <- x$table
confmat <- round(confmat %*% solve(diag(colSums(confmat))), 3)
confmat <- as.table(confmat)
colnames(confmat) <- rownames(confmat)
names(dimnames(confmat))[2] <- "Reference"
print(confmat)
print(x$table)
class(confmat)
df_plot <- as.data.frame(confmat)
gg <- ggplot(df_plot) + geom_tile(aes(x = Prediction, y = Reference, fill = Freq))
gg <- gg + scale_fill_gradient2(low = muted("green"), high = muted("blue"))
gg

#Again on validate set
levels(df_test$y)
predictions <- predict(object = objModel, df_valid[, predictorsNames], type = 'raw') #type='prob')
print(postResample(pred = predictions, obs = as.factor(df_valid[,outcomeName])))
x <- confusionMatrix(predictions, df_valid[, outcomeName])
confmat <- x$table
confmat <- round(confmat %*% solve(diag(colSums(confmat))), 3)
confmat <- as.table(confmat)
colnames(confmat) <- rownames(confmat)
names(dimnames(confmat))[2] <- "Reference"
print(confmat)
print(x$table)
class(confmat)
df_plot <- as.data.frame(confmat)
gg <- ggplot(df_plot) + geom_tile(aes(x = Prediction, y = Reference, fill = Freq))
gg <- gg + scale_fill_gradient2(low = muted("green"), high = muted("blue"))
gg
#

df_plot <- out_y[[1]][ind_tsStart:ind_tsFinish, c("Index", "ts", "SigUpWin")]
ind_up <- which(df_plot$SigUpWin == "Uptrend Start")
ind_dn <- which(df_plot$SigUpWin == "Uptrend Stop")
df_plot <- df_plot[valid_rows, ]
ind <- which(df_plot$SigUpWin %in% c(""))
df_plot$SigUpWin
gg <- ggplot(df_plot, aes(x = Index, y = ts)) + geom_line()
gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Act")], color = "blue", alpha = 0.3)
gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Uptrend Stop")], color = "red", alpha = 0.3)

# gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Uptrend Start")], color = "green", alpha = 0.3)
# gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Uptrend Stop")], color = "red", alpha = 0.3)
# gg <- gg + geom_vline(xintercept = df_plot$Index[which(df_plot$SigUpWin == "Uptrend Start")], color = "blue", linetype = "dotted")
# gg <- gg + geom_vline(xintercept = df_plot$Index[which(df_plot$SigUpWin == "Uptrend Stop")], color = "violet", linetype = "dotted")
gg

#
predictions <- predict(object = objModel, df_testRecent[, predictorsNames], type = 'prob') #type='prob')
df_plot$pred <- predictions$Act
mark <- which(predictions$Act > 0.75)
df_plot <- df_plot %>% gather(Type, Value, ts:pred)
gg <- ggplot(df_plot, aes(x = Index, y = Value)) + geom_line()
gg <- gg + facet_wrap(~Type, ncol = 1, scales = "free")
gg <- gg + geom_vline(xintercept = df_plot$Index[mark])
gg <- gg + geom_vline(xintercept = df_plot$Index[mark])
gg









predictions <- predict(object = objModel, df_valid[, predictorsNames], type = 'raw') #type='prob')
print(postResample(pred = predictions, obs = as.factor(df_valid[,outcomeName])))

df_plot <- out_y[[1]][ind_tsStart:ind_tsFinish, c("Index", "ts")]
df_plot <- df_plot[valid_rows,]
gg <- ggplot(df_plot, aes(x = Index, y = ts)) + geom_line()
gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Act")])
gg



predictions <- predict(object = objModel, df_testRecent[, predictorsNames], type = 'prob') #type='prob')

df_plot <- out_y[[1]][ind_tsStart:ind_tsFinish, c("Index", "ts")]
df_plot <- df_plot[valid_rows,]
df_plot$pred <- predictions$Act
df_plot <- df_plot %>% gather(Type, Value, ts:pred)
gg <- ggplot(df_plot, aes(x = Index, y = Value)) + geom_line()
#gg <- gg + geom_vline(xintercept = df_plot$Index[which(predictions == "Act")])
gg <- gg + facet_wrap(~Type, ncol = 1, scales = "free")
gg

















#----------------------------
#redux with only important variables
q <- as.numeric(quantile(var_importance$importance[, 1], probs = 0.85))
u_name <- rownames(var_importance$importance)
u_val <- var_importance$importance$Overall
impvar_names <- u_name[which(u_val > q)]
rm(u_val, u_name)
y_col <- which(colnames(df_mod) == "y")
keep_cols <- c(which(colnames(df_mod) %in% impvar_names), y_col)
df_mod2 <- df_mod[, keep_cols]
#-----------
predictorsNames <- colnames(df_mod2)[1:(ncol(df_mod2) - 1)]
outcomeName <- colnames(df_mod2)[ncol(df_mod2)]
set.seed(1234)
splitIndex <- createDataPartition(df_mod2[, outcomeName], p = .5, list = FALSE, times = 1)
#splitIndex <- createTimeSlices(df_mod[, outcomeName], 30)
df_train <- df_mod2[ splitIndex,]
df_test <- df_mod2[-splitIndex,]  
objControl <- trainControl(method = 'repeatedcv', number = 3)#, repeats = 3)
objModel <- train(df_train[, predictorsNames], df_train[, outcomeName],
                  method = this_method, trControl = objControl)
#-----------
var_importance <- varImp(objModel, scale = T)
print(var_importance)
#plot(var_importance)
#-----------------------------
# probabilities ("prob") or integer ("raw")
predictions <- predict(object = objModel, df_test[, predictorsNames], type = 'raw') #type='prob')
print(postResample(pred = predictions, obs = as.factor(df_test[,outcomeName])))
#head(predictions)
x <- confusionMatrix(predictions, df_test[, outcomeName])
confmat <- x$table
confmat <- round(confmat %*% solve(diag(colSums(confmat))), 3)
confmat <- as.table(confmat)
colnames(confmat) <- rownames(confmat)
names(dimnames(confmat))[2] <- "Reference"
print(confmat)
print(x$table)
class(confmat)
df_plot <- as.data.frame(confmat)
gg <- ggplot(df_plot) + geom_tile(aes(x = Prediction, y = Reference, fill = Freq))
gg <- gg + scale_fill_gradient2(low = muted("green"), high = muted("blue"))
gg
#----------------------------










df_pca <- df
GroupInfo <- read.csv("GlobTrndsID.csv", stringsAsFactors = F)
GroupInfo$X <- NULL
colnames(GroupInfo)[1] <- "name"
df_pca <- merge(df_pca, GroupInfo, by = "name")
rownames(df_pca) <- df_pca$name
df_pca$name <- NULL
df_pca$General.Type <- as.factor(df_pca$General.Type)
ind_num <- which(!(colnames(df_pca) %in% c("General.Type", "Sub.type", "Specific.Track")))
#---------------
res <- PCA(df_pca[, ind_num])
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res, habillage = df_pca$General.Type)
#====================================
df_dydx_pca <- df_dydxmu
#df_dydx_pca <- df_ldydxcv
df_dydx_pca$Month <- month(df_dydx_pca$Index)
gathercols <- colnames(df_dydx_pca)[1:(ncol(df_dydx_pca) - 2)]
df_dydx_pca <- gather_(df_dydx_pca, "name", "ts", gathercols)
df_dydx_pca <- df_dydx_pca[-which(is.na(df_dydx_pca$ts)), ]
df_volat <- df_dydx_pca
df_volat <- df_volat %>% group_by(name, Month) %>% summarize(ts_mu = mean(ts, na.rm = T), ts_sd = sd(ts, na.rm = T))
df_volat$cv <- df_volat$ts_sd / df_volat$ts_mu
#df_volat$date <- as.yearmon(paste(df_volat$Year, df_volat$Month, sep = "-"))
df_volat <- df_volat[which(is.nan(df_volat$cv) == F),]
df_volat <- df_volat[, c("name", "Month", "cv")]
df_volat <- df_volat %>% spread(Month, cv)
colnames(df_volat)[2:13] <- month.abb[as.numeric(colnames(df_volat)[2:13])]
#----------------------
df_pca <- merge(df_volat, GroupInfo, by = "name")
#df_pca <- df_pca[-which(df_pca$name == "VXX"), ]
#df_pca <- df_pca[-which(df_pca$name == "PGD"), ]
rownames(df_pca) <- df_pca$name
df_pca$name <- NULL
df_pca$General.Type <- as.factor(df_pca$General.Type)
ind_num <- which(!(colnames(df_pca) %in% c("General.Type", "Sub.type", "Specific.Track")))
#---------------
res <- PCA(df_pca[, ind_num])
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res, habillage = df_pca$General.Type)
#---------------
# mc <- Mclust(df_pca[, ind_num])
# summary(mc)
# fviz_cluster(mc, frame.type = "norm", geom = "point")




































y_train <- recode(df_train$y, Hold = "0", `Uptrend Start` = "1",
                  `Uptrend Stop` = "2", `Uptrend False Start` = "3",
                  `Uptrend False Stop` = "4")
y_test <- recode(df_test$y, Hold = "0", `Uptrend Start` = "1",
                 `Uptrend Stop` = "2", `Uptrend False Start` = "3",
                 `Uptrend False Stop` = "4")
y_train <- as.numeric(as.character(y_train))
y_test <- as.numeric(as.character(y_test))

mat_train <- as.matrix(df_train[, -1])
mat_test <- as.matrix(df_test[, -1])
l_train <- list(feats = mat_train, y = y_train)
l_test <- list(feats = mat_test, y = y_test)

dtrain <- xgb.DMatrix(data = l_train$feats, label = l_train$y)
dtest <- xgb.DMatrix(data = l_test$feats, label = l_test$y)

# xgbcv <- xgb.cv(params = params, data = dtrain, nrounds = 100, 
#                 nfold = 5, showsd = T, stratified = T, 
#                 print_every_n = 10, early_stop_round = 20,
#                 num_class = 12,
#                 maximize = F)
# 
# min(xgbcv$test.error.mean)

# k-fold cross validation, with timing
params <- list(booster = "gbtree", objective = "multi:softmax", eta = 0.1,
               gamma = 0, max_depth = 25, min_child_weight = 1, subsample = 0.5,
               colsample_bytree = 0.5, num_class = 5)
nround.cv = 500
xgboost.cv <- xgb.cv(params = params,
                     data = dtrain,
                     nfold = 10,
                     nrounds = nround.cv,
                     prediction = TRUE,
                     verbose=1)

# index of maximum auc:
max.auc.idx = which.max(xgboost.cv$dt[, test.auc.mean]) 
max.auc.idx 
## [1] 493
# minimum merror
xgboost.cv$dt[max.auc.idx,]

# real model fit training, with full data
xgb.bst <- xgboost(param=param, data=train.matrix, label=train$Survived, 
                   nrounds=max.auc.idx, verbose=1)
pred <- predict(xgb.bst,test.matrix)

















params <- list(booster = "gbtree", objective = "multi:softmax", eta = 0.1,
               gamma = 0, max_depth = 25, min_child_weight = 1, subsample = 0.5,
               colsample_bytree = 0.5)


watchlist <- list(train = dtrain, test = dtest)

#xgb <- xgboost(
xgb <- xgb.train(params = params,
                 data = dtrain,
                 nround = 35,
                 watchlist = watchlist,
                 seed = 1,
                 eval_metric = "merror",
                 num_class = 5,
                 verbose = 1
)

#xgb1 <- xgb.train (params = params, data = dtrain, nrounds = 79, watchlist = list(val=dtest,train=dtrain), print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")
#model prediction
xgbpred <- predict(xgb, dtest)
str(xgbpred)

x <- confusionMatrix (xgbpred, y_test)
print(x)
confmat <- x$table
confmat <- round(confmat %*% solve(diag(colSums(confmat))), 3)
confmat <- as.table(confmat)
colnames(confmat) <- rownames(confmat)
names(dimnames(confmat))[2] <- "Reference"
print(confmat)
print(x$table)
class(confmat)
df_plot <- as.data.frame(confmat)
gg <- ggplot(df_plot) + geom_tile(aes(x = Prediction, y = Reference, fill = Freq))
gg <- gg + scale_fill_gradient2(low = muted("green"), high = muted("blue"))
gg




mat <- xgb.importance (feature_names = colnames(mat_train), model = xgb)
xgb.plot.importance (importance_matrix = mat[1:20]) 


# predict values in test set
y_pred <- predict(xgb, data.matrix(df_test[,-1]))
head(y_pred)
# Get the feature real names
names <- dimnames(data.matrix(df_train[, -1]))[[2]]

# Compute feature importance matrix
importance_matrix <- xgb.importance(names, model = xgb)
# Nice graph
xgb.plot.importance(importance_matrix[1:10,])









#names(getModelInfo())
# boosted tree model (gbm) adjust learning rate and and trees
objControl <- trainControl(method = 'cv', number = 3, returnResamp = 'none', classProbs = TRUE)
# gbmGrid <-  expand.grid(interaction.depth =  c(1, 5, 9),
#                         n.trees = 50,
#                         shrinkage = 0.01)
# run model
#objModel <- train(Xmat_train[,predictorsNames], Xmat_train[,outcomeName], method='gbm', trControl=objControl, tuneGrid = gbmGrid, verbose=F)
#objControl <- trainControl(method='cv', classProbs = TRUE)
#names(getModelInfo())
#this_method <- "glmnet"
#this_method <- "glm"
#this_method <- "nb"
this_method <- "gbm" #good
#this_method <- "naive_bayes"
#this_method <- "svmRadial"
#this_method <- "cforest"
objModel <- train(Modmat_train[, predictorsNames], Modmat_train[, outcomeName],
                  method = this_method,
                  trControl = objControl,
                  preProc = c("center", "scale"))

predictions <- predict(object=objModel, Modmat_test[,predictorsNames], type='raw')
head(predictions)
print(postResample(pred = predictions, obs = as.factor(Modmat_test[, outcomeName])))

# probabilities 
predictions <- predict(object = objModel, Modmat_test[, predictorsNames], type = 'prob')
head(predictions)
postResample(pred = predictions_buy[[2]], obs = ifelse(Modmat_test[,outcomeName] == 'buy', 1 , 0))


varImp(objModel,scale=T)
var_importance <- varImp(objModel,scale=T)$importance
#class(var_importance$importance)
#plot(var_importance)
q <- as.numeric(quantile(var_importance[, 1], probs = 0.75))
varnames <- rownames(var_importance)
#relinf <- summary(objModel)[,"rel.inf"]
#print(relinf)
#vars_screened <- summary(objModel)[which(relinf > q),"var"]
vars_screened <- varnames[which(var_importance[, 1] > q)]
vars_screened <- as.character(vars_screened)
objModel <- train(Modmat_train[, vars_screened], Modmat_train[, outcomeName],
                  method=this_method, 
                  trControl=objControl,  
                  metric = "ROC",
                  preProc = c("center", "scale"))





















df_outPCA <- as.data.frame(res$ind$coord)
get_clust_tendency(df_cpO[, ind_num], n = 50, gradient = list(low = "steelblue",  high = "white"))
fviz_nbclust(df_cpO[, ind_num], kmeans, method = "gap_stat")
pam.res <- pam(df_cpO[, ind_num], 5)
fviz_cluster(pam.res)
mc <- Mclust(df_outPCA)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")

#-----------
# calculate correlation matrix
#for sectors
zdf_cp_sectr <- zdf_cp_w %>% group_by(Sector) %>% summarize_all(mean)
zdf_cp_sectr$Industry <- NULL
rownames(zdf_cp_sectr) <- zdf_cp_sectr$Sector
zdf_cp_sectr$Sector <- NULL
zdf_cp_sectr <- as.data.frame(t(zdf_cp_sectr))
corMatMy <- cor(zdf_cp_sectr)
corrplot(corMatMy, order = "hclust")
#Apply correlation filter at 0.70,
highlyCor <- colnames(df_cp_sectr)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
print(highlyCor)

# #for industries
# df_cp_indus <- df_cp_w %>% group_by(Industry) %>% summarize_all(mean)
# df_cp_indus$Sector <- NULL
# #rownames(df_cp_indus) <- df_cp_indus$Industry
# df_cp_indus$Industry <- NULL
# df_cp_indus <- as.data.frame(t(df_cp_indus))
# corMatMy <- cor(df_cp_indus)
# corrplot(corMatMy, order = "hclust")
# #Apply correlation filter at 0.70,
# highlyCor <- colnames(df_prop)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
# print(highlyCor)
#-----------


out_rf <- randomForest(t(df_cp_w))
# str(out_rf)
k = 2
MDSplot(out_rf, df_pam$pamClust, k)

mc <- Mclust(df)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")




























#fviz_pca_ind(res, axes = c(1,2), label= "none", habillage = df_$Class)#, addEllipses = T, ellipse.level=0.95)
#df_outPCA <- data.frame(res$ind$coord, res$ind$contrib, res$ind$cos2)
df_outPCA <- data.frame(res$ind$contrib)
get_clust_tendency(df_outPCA, n = 50, gradient = list(low = "steelblue",  high = "white"))
get_clust_tendency(zdf, n = 50, gradient = list(low = "steelblue",  high = "white"))
get_clust_tendency(df, n = 50, gradient = list(low = "steelblue",  high = "white"))
#Optimum number of clusters
fviz_nbclust(df_num, kmeans, method = "gap_stat")
fviz_nbclust(df, kmeans, method = "gap_stat")
# fviz_nbclust(zdf, kmeans, method = "gap_stat")
# fviz_nbclust(df, kmeans, method = "gap_stat")


pam.res <- pam(df_num, 10)
fviz_cluster(pam.res)
df_pam <- data.frame(Clust = pam.res$clustering)
df_pam <- cbind(df, df_pam)
df_pam$pamClust <- as.factor(df_pam$pamClust)
out_rf <- randomForest(df)
# str(out_rf)
k = 2
MDSplot(out_rf, df_pam$pamClust, k)

mc <- Mclust(df)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")
df_mc <- data.frame(Clust = mc$classification)

df_mc$mcClust <- as.factor(df_mc$mcClust)

#df_outPCA2 <- cbind(df_outPCA, df_mc)
df_outPCA2$mcClust <- as.factor(df_outPCA2$mcClust)
out_rf <- randomForest(df_outPCA)
# str(out_rf)
k = 3
MDSplot(out_rf, df_mc$mcClust, k)
out_rf.mds <- stats::cmdscale(1 - out_rf$proximity, eig = TRUE, k)
plot(out_rf.mds$points)
str(out_rf.mds)
class(out_rf.mds$points)
df_rf <- as.data.frame(out_rf.mds$points)
head(df_rf)
ggplot(df_rf, aes(x = V1, y = V2)) + geom_point()

mc <- Mclust(df_rf)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")
df_mc <- as.data.frame(mc$classification)
colnames(df_mc)[1] <- "mcClust"
head(df_mc)
rownames(df_mc)[which(df_mc$mcClust == 1)]
rownames(df_mc)[which(df_mc$mcClust == 3)]
rownames(df_mc)[which(df_mc$mcClust == 4)]
#
dfred_rf <- df[which(df_mc$mcClust == 1),]
fviz_nbclust(dfred_rf, pam, method = "gap_stat")
pam.res <- pam(dfred_rf, 5)
fviz_cluster(pam.res)
# fviz_nbclust(dfred_rf, kmeans, method = "gap_stat")
# km.res <- kmeans(dfred_rf, 3, nstart = 25)
# fviz_cluster(km.res, data = dfred_rf, frame.type = "convex") + theme_minimal()

#
fviz_nbclust(df_rf, kmeans, method = "gap_stat")
km.res <- kmeans(df_rf, 8, nstart = 25)
# Visualize
fviz_cluster(km.res, data = df_rf, frame.type = "convex") + theme_minimal()
str(km.res)
df_km <- as.data.frame(km.res$cluster)
colnames(df_km)[1] <- "kmClust"
head(df_km)
rownames(df_km)[which(df_km$kmClust == 1)]

# library(sp)
# dens <- kde2d(df_rf[, 1], df_rf[, 2], n = 100, lims=c(-0.3, 0.3, -0.3, 0.3))
# levels <- 0.9
# ls <- contourLines(dens, level=levels)
# inner <- point.in.polygon(df_rf[, 1], df_rf[, 2], ls[[2]]$x, ls[[2]]$y)
# out <- point.in.polygon(df_rf[, 1], df_rf[, 2], ls[[1]]$x, ls[[1]]$y)
# df_rf$region <- factor(inner + out)
# plot(V1 ~ V2, col=region, data=df_rf, pch=15)
# contour(dens, levels=levels, labels=prob, add=T)

# filled.contour(f1)
#z <- round(f1$z,3)
# dens <- interp.surface(f1, df_rf[,c(1,2)])
# head(dens)
