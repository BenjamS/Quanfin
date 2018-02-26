setwd("D:/OneDrive - CGIAR/Documents")
source('./tsSlope.R', echo=TRUE)
#source('./captureEvents.R', echo=TRUE)
source('./getGlobTrnds.R', echo=TRUE)
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

fromdate<-"2012-01-01"; todate <- "2018-02-26"
getGlobTrnds(fromdate, todate)
rmcol <- which(colnames(cpgtetfmat) %in% c("EXS1.DE", "^IRX", "EWH", "AIA", "XTN", "VXX"))
xts_cp <- cpgtetfmat[, -rmcol]
xts_vol <- volgtetfmat[, -rmcol]
namevec <- colnames(xts_cp)
n_names <- length(namevec)


o <- apply(xts_cp, 2, function(x) length(which(is.na(x))))
table(o)
rmcols <- which(o > 57)
colnames(xts_cp)[rmcols]
xts_cp <- xts_cp[, -rmcols]
xts_vol <- xts_vol[, -rmcols]
o <- apply(xts_cp, 1, function(x) length(which(is.na(x))))
table(o)
xts_cp <- na.spline(xts_cp)
xts_vol <- na.spline(xts_vol)
o <- apply(xts_cp, 1, function(x) length(which(is.na(x))))
table(o)
#---------------
n_ts <- ncol(xts_cp)
#---------------
list_df <- list()
for(i in 1:n_ts)
{
  this_sec <- colnames(xts_cp)[i]
  out_df <- fortify(VWMA(xts_cp[, i], xts_vol[, i]))
  colnames(out_df)[2] <- this_sec
  list_df[[i]] <- out_df
}
df_vwma <- join_all(list_df)
rm_rows <- which(is.na(df_vwma[, 2]))
df_vwma <- df_vwma[-rm_rows, ]
datevec <- df_vwma$Index
df_cp <- fortify(xts_cp[-rm_rows, ])
df_dt_vwma <- df_cp[, -1] - df_vwma[, -1]
#--------------
df_ema <- as.data.frame(apply(xts_cp, 2, function(x) EMA(x, n = 21)))
rm_rows <- which(is.na(df_ema[, 2]))
df_ema <- df_ema[-rm_rows, ]
df_cp <- fortify(xts_cp[-rm_rows, ])
df_dt_ema <- df_cp[, -1] - df_ema
#--------------
ts <- df_cp$FXI
ts_diff <- diff(ts, differences = 3)
df_plot <- data.frame(ts, ts_diff, ts_dt_vwma = df_dt_vwma$FXI, ts_dt_ema = df_dt_ema$FXI)



sd_diff <- sd(df_out$tsDiff, na.rm = T)
ind_diffs <- which(abs(df_out$tsDiff) > 3 * sd_diff)
df_plot <- df_out %>% gather(Series, Value, ts:dydx_sd)
gg <- ggplot(df_plot, aes(x = Date, y = Value))
gg <- gg + geom_line()
gg <- gg + geom_vline(xintercept = df_out$Date[ind_diffs], color = "green")
gg <- gg + facet_wrap(~Series, ncol = 1, scales = "free")
gg



gg1 <- ggplot(df_vwma_plot, aes(x = Index, y = SPY)) + geom_line()
gg2 <- ggplot(df_dt_plot, aes(x = Index, y = SPY)) + geom_line()
grid.arrange(gg1, gg2, ncol = 1)



#--------
selectnames <- c("FM", "UUP", "USO", "TTEK", "GRN", "SPY")
df_mnth <- df_evwma
df_mnth <- df_mnth[, which(colnames(df_mnth) %in% selectnames)]
df_mnth$date <- datevec
gathercols <- colnames(df_mnth)[1:(ncol(df_mnth) - 1)]
df_mnth <- gather_(df_mnth, "name", "cp", gathercols)
df_mnth$Year <- year(df_mnth$date)
df_mnth$Month <- month(df_mnth$date)
df_mnth$Day <- yday(df_mnth$date)
df_mnth <- df_mnth[-which(is.na(df_mnth$cp)), ]
df_volat <- df_mnth
df_volat <- df_volat %>% group_by(name, Year, Month) %>% summarize(cp_mu = mean(cp, na.rm = T), cp_sd = sd(cp, na.rm = T))
df_volat$cv <- df_volat$cp_sd / df_volat$cp_mu
df_volat$date <- as.yearmon(paste(df_volat$Year, df_volat$Month, sep = "-"))
#--------
df_volat$date <- as.Date(df_volat$date)
gg <- ggplot(df_volat, aes(x = date, y = cv, group = name, color = name)) + geom_line()
gg <- gg + facet_wrap(~Year, ncol = 1, scales = "free")
gg
#df_mnth <- df_mnth %>% group_by(Year, name) %>% mutate(scale(cp))




#---























gg <- ggplot(df_mnth, aes(x = Day, y = `scale(cp)`, group = name, color = name)) + geom_line()
gg <- gg + facet_wrap(~Year, ncol = 1)
gg
df_mnth$Year <- as.factor(df_mnth$Year)
gg <- ggplot(df_mnth, aes(x = Day, y = `scale(cp)`, group = Year, color = Year)) + geom_line()
gg <- gg + facet_wrap(~name, ncol = 2)
gg
#--------







































df <- df_dt[, which(!(colnames(df_dt) %in% c("PGD", "VXX", "EWH")))]
df <- scale(df)
df$date <- datevec
df$Month <- month(df$date, label = T)
df$date <- NULL
df <- df %>% group_by(Month) %>% summarise_each(funs(mean(., na.rm = TRUE))) #summarize_all(mean)
class(df)
df <- as.data.frame(df)
rownames(df) <- df$Month
df$Month <- NULL
df <- as.data.frame(t(df))
#df$Symb <- rownames(df)

res <- PCA(df)
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res)
#---------------
mc <- Mclust(df)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")
#---------------
df_cor <- df
corMatMy <- cor(df_cor)
corrplot(corMatMy, order = "hclust")
#--------------------
df_cor <- df_dt[-c(1:10), which(colnames(df_dt) != "PGD")]
corMatMy <- cor(df_cor)
corrplot(corMatMy, order = "hclust")
#--------------------
df














#---------------
df_cp <- fortify(apply.weekly(xts_cp, mean))
#df_cp_w <- fortify(apply.monthly(df_cp, mean))
#df_cp_w <- fortify(apply.yearly(df_cp, mean))
df_cp$Month <- months(as.Date(df_cp$Index))
df_cp$Index <- NULL
zdf_cp <- scale(df_cp[, -ncol(df_cp)])
zdf_cp <- data.frame(Month = df_cp$Month, zdf_cp)
zdf_cp <- zdf_cp %>% group_by(Month) %>% summarize_all(mean)
class(zdf_cp)
zdf_cp <- as.data.frame(zdf_cp)
rownames(zdf_cp) <- zdf_cp$Month
zdf_cp$Month <- NULL
zdf_cp <- as.data.frame(t(zdf_cp))
zdf_cp$Symb <- rownames(zdf_cp)
#--
GroupInfo <- read.csv("GlobTrndsID.csv", stringsAsFactors = F)
GroupInfo$X <- NULL
colnames(GroupInfo)[1] <- "Symb"
zdf_cp <- merge(zdf_cp, GroupInfo, by = "Symb")
rownames(zdf_cp) <- zdf_cp$Symb
#--
zdf_cp$Symb <- NULL
ind_num <- which(!(colnames(zdf_cp) %in% c("General.Type", "Sub.type", "Specific.Track")))
#--
corMatMy <- cor(zdf_cp[, ind_num])
corrplot(corMatMy, order = "hclust")
#Apply correlation filter at 0.70,
# highlyCor <- colnames(zdf_cp)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
# print(highlyCor)
#--
res <- PCA(zdf_cp[, ind_num])
fviz_screeplot(res, ncp=5)
zdf_cp$General.Type <- as.factor(zdf_cp$General.Type)
#fviz_pca_ind(res, habillage = zdf_cp$General.Type)
fviz_pca_biplot(res, habillage = zdf_cp$General.Type)#,
#-------------------
ema_per <- 21
datevec <- index(xts_cp)
xts_cpEMA <- xts(apply(xts_cp, MARGIN = 2, FUN = "EMA", n = ema_per), datevec)
class(xts_cpEMA)


rmcolsEVWMA <- which(colnames(xts_cp) %in% c("PGD", "^IRX"))
xts_cp_inEVWMA <- xts_cp[, -rmcolsEVWMA]
xts_vol_inEVWMA <- xts_vol[, -rmcolsEVWMA]
outlist <- list()
for(i in 1:ncol(xts_cp_inEVWMA)){this_symb <- colnames(xts_cp_inEVWMA)[i] 
this_series <- fortify(EVWMA(xts_cp_inEVWMA[, i],xts_vol_inEVWMA[, i]))
colnames(this_series)[2] <- this_symb;print(head(this_series))
outlist[[i]] <- this_series}

for(i in 1:length(outlist)){print(tail(outlist[[i]]))}

xts_cpEVWMA <- join_all(outlist)

#xts_cpEMA <- xts(apply(xts_cpW, MARGIN = 2, FUN = "EVWMA", volume = ), index(xts_cpW))
#num_na <- ema_per - 1
xts_cpEMA <- xts_cpEMA[c(ema_per:nrow(xts_cpEMA)),]
class(xts_cpEMA)
#volmat <- volmat[c(ema_per:nrow(volmat)),]; volmat <- round(volmat, 4)
# library(pracma)
# sOsc_per <- 13
# xts_cpsOsc <- xts(apply(xts_cpEMA, MARGIN = 2, FUN = "slopeOsc", wind = sOsc_per), index(xts_cpEMA))
# xts_cpsOsc <- xts_cpsOsc[c((sOsc_per + 1):nrow(xts_cpsOsc)),]
# class(xts_cpsOsc)
# df_cpsOsc <- fortify(xts_cpsOsc)
# rownames(df_cpsOsc) <- df_cpsOsc$Index
# df_cpsOsc$Index <- NULL
# df_cpO <- as.data.frame(t(df_cpsOsc))
# df_cpO$Symb <- rownames(df_cpO)
# df_cpO <- merge(df_cpO, df_raw_sectInd, by = "Symb")
# rownames(df_cpO) <- df_cpO$Symb
# df_cpO$Symb <- NULL
# ind_num <- which(!(colnames(df_cpO) %in% c("Sector", "Industry")))
# res <- PCA(df_cpO[, ind_num])
# df_cpO$Sector <- as.factor(df_cpO$Sector)
# fviz_pca_ind(res, habillage = df_cpO[, "Sector"])
# fviz_pca_biplot(res, habillage = df_cpO[, "Sector"])
# ggpairs(df_outPCA[, 1:5], aes(alpha = 0.4))
#--
# df_cpO_sectr <- df_cpO %>% group_by(General.Type) %>% summarize_all(mean)
# df_cpO_sectr$Industry <- NULL
# rownames(df_cpO_sectr) <- df_cpO_sectr$Sector
# df_cpO_sectr$Sector <- NULL
# df_cpO_sectr <- as.data.frame(t(df_cpO_sectr))
# corMatMy <- cor(df_cpO_sectr)
# corrplot(corMatMy, order = "hclust")
# #Apply correlation filter at 0.70,
# highlyCor <- colnames(df_cpO_sectr)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
# print(highlyCor)
#------------------
# class(xts_cpsOsc)
# monthvec <- months(index(xts_cpsOsc))
# df_cpsOsc <- data.frame(Month = monthvec, fortify(xts_cpsOsc))
# df_cpsOsc$Index <- NULL
# df_cpsOsc_month <- df_cpsOsc %>% group_by(Month) %>% summarize_all(mean)
# rownames(df_cpsOsc_month) <- df_cpsOsc_month$Month
# df_cpsOsc_month$Month <- NULL
# df_cpsOsc_month <- as.data.frame(t(df_cpsOsc_month))
# corMatMy <- cor(df_cpsOsc_month)
# corrplot(corMatMy, order = "hclust")
# #Apply correlation filter at 0.70,
# highlyCor <- colnames(df_cpsOsc_month)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
# print(highlyCor)
#-----------------

#=============================================
datevec <<- index(xts_cpEMA)
n_stocks <- ncol(xts_cpEMA)
outlist <- list()
outlist_buy <- list()
outlist_sell <- list()
for(i in 1:n_stocks)
{
  indata <- xts_cpEMA[, i]
  out <- getBuySellbinryVar(indata)
  # buy_sig <- as.integer(as.character(out[[1]]))
  # sell_sig <- as.integer(as.character(out[[2]]))
  buy_sig <- out[[3]]
  sell_sig <- out[[4]]
  df_buysell <- data.frame(Date = index(xts_cpEMA), buySig = buy_sig, sellSig = sell_sig)
  df_buysell$Month <- months(df_buysell$Date)
  df_buysell_m <- df_buysell %>% group_by(Month) %>% summarize(buy = sum(buySig), sell = sum(sellSig))
  df_buysell_m <- df_buysell_m %>% gather(Sig, Value, buy:sell)
  #df_buysell_m <- df_buysell_m %>% gather(Sig, Value, buy:sell) %>% unite(MonthSig, Month, Sig)
  #df_buysell_m <- df_buysell_m %>% spread(MonthSig, Value)
  df_buysell_m <- df_buysell_m %>% spread(Month, Value)
  #rownames(df_buysell_m) <- colnames(indata)
  rownames(df_buysell_m) <- paste(df_buysell_m$Sig, colnames(indata))
  df_buy_m <- df_buysell_m[1,]
  df_sell_m <- df_buysell_m[2,]
  rownames(df_buy_m) <- colnames(indata)
  rownames(df_sell_m) <- colnames(indata)
  outlist_buy[[i]] <- df_buy_m
  outlist_sell[[i]] <- df_sell_m
  
}

#df_SigCA <- do.call(rbind, outlist)
df_SigCA_buy <- do.call(rbind, outlist_buy)
df_SigCA_sell <- do.call(rbind, outlist_sell)
df_SigCA <- rbind(df_SigCA_buy, df_SigCA_sell)
#df_SigCA <- df_SigCA_sell
df_SigCA$Symb <- rownames(df_SigCA)
df_SigCA$Symb <- gsub("1", "", df_SigCA$Symb)
df_SigCA <- merge(df_SigCA, GroupInfo, by = "Symb")
df_SigCA$Symb <- NULL
df_SigCA$x <- paste(df_SigCA$Specific.Track, df_SigCA$Sig)
rownames(df_SigCA) <- df_SigCA$x
df_SigCA$x <- NULL
df_SigCA$Sig <- NULL
df_SigCA$Specific.Track <- NULL
#-
# cormat <- as.matrix(cor(df_SigCA))
# corrplot(cormat, is.corr=T)
#-
ind_num <- which(!(colnames(df_SigCA) %in% c("General.Type", "Sub.type")))
df_SigCA$General.Type <- as.factor(df_SigCA$General.Type)
df_SigCA$Sub.type <- as.factor(df_SigCA$Sub.type)
res <- CA(df_SigCA[, ind_num], ncp = 10)
#summary(res)
nrow(res$row$coord)
nrow(df_SigCA)
fviz_screeplot(res, addlabels = TRUE)
fviz_ca(res)
fviz_ca_row(res, habillage = df_SigCA[,"General.Type"])
col <- get_ca_col(res)
corrplot(col$cos2, is.corr=FALSE)
#-
df_outCA <- as.data.frame(res$row$coord)
get_clust_tendency(df_outCA, n = 50, gradient = list(low = "steelblue",  high = "white"))
fviz_nbclust(df_outCA, kmeans, method = "gap_stat")
pam.res <- pam(df_outCA, 5)
fviz_cluster(pam.res)
mc <- Mclust(df_outCA)
summary(mc)
fviz_cluster(mc, frame.type = "norm", geom = "point")

















#class(df_buysell$Date)
Modmat <- df_cpsOsc
Modmat$Month <- NULL
Modmat$Decision <- 0
Modmat$Decision[which(buy_sig == 1)] <- 1
Modmat$Decision[which(sell_sig == 1)] <- 1
#Modmat <- merge(XTmat, df_buysell, by = "Date")
# Modmat <- merge(XTmat, df_buysell, by = "Date")
# Modmat <- merge(Tmat, df_buysell, by = "Date")
Modmat$buySig <- ifelse(Modmat$buySig == 1,'buy', 'hold')
Modmat$buySig <- as.factor(Modmat$buySig)
outcomeName_buy <- "buySig"
Modmat$sellSig <- ifelse(Modmat$sellSig == 1,'sell', 'hold')
Modmat$sellSig <- as.factor(Modmat$sellSig)
outcomeName_sell <- "sellSig"
Modmat$anySig <- "hold"
ind_buy <- which(Modmat$buySig == "buy"); ind_sell <- which(Modmat$sellSig == "sell")
ind_event <- c(ind_buy, ind_sell)
Modmat$anySig[ind_event] <- "event"
Modmat$anySig <- as.factor(Modmat$anySig)
outcomeName_Any <- "anySig"


predictorsNames <- colnames(Modmat)[2:(ncol(Modmat) - 3)]
traindat_pctot <- .42
indtrain_beg <- 1
indtrain_end <- round(length(datevec)*traindat_pctot)
indtest_beg <- indtrain_end + 1
indtest_end <- nrow(Modmat)
Modmat_train <- Modmat[indtrain_beg:indtrain_end,]
Modmat_test <- Modmat[indtest_beg:indtest_end,]

#names(getModelInfo())
# boosted tree model (gbm) adjust learning rate and and trees
objControl <- trainControl(method='cv', number=3, returnResamp='none', summaryFunction = twoClassSummary, classProbs = TRUE)
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
#this_method <- "gbm" #good
this_method <- "naive_bayes"
#this_method <- "svmRadial"
#this_method <- "cforest"
objModel <- train(Modmat_train[, predictorsNames], Modmat_train[,outcomeName_buy],
                  method=this_method,
                  trControl=objControl,
                  preProc = c("center", "scale"))

# predictions <- predict(object=objModel, Modmat_test[,predictorsNames], type='raw')
# head(predictions)
# print(postResample(pred=predictions, obs=as.factor(Modmat_test[,outcomeName_buy])))

# probabilities 
predictions_buy <- predict(object=objModel, Modmat_test[,predictorsNames], type='prob')
head(predictions_buy)
postResample(pred=predictions_buy[[2]], obs=ifelse(Modmat_test[,outcomeName_buy]=='buy',1,0))


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
objModel <- train(Modmat_train[, vars_screened], Modmat_train[,outcomeName_buy],
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
