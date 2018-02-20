library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(GGally)
library(cluster)
library(mclust)
library(randomForest)
library(corrplot)

df_raw <- read.csv("fundamentals.csv", stringsAsFactors = F)
df_raw$X <- NULL
colnames(df_raw)[1] <- "Symb"
#rownames(df_raw)
#length((df_raw$Period.Ending[grep("2016-09", df_raw$Period.Ending)]))

df_raw_sectInd <- read.csv("securities.csv", stringsAsFactors = F)
df_raw_sectInd <- df_raw_sectInd[, c("Ticker.symbol", "GICS.Sector", "GICS.Sub.Industry")]
head(df_raw_sectInd)
colnames(df_raw_sectInd) <- c("Symb", "Sector", "Industry")

df <- df_raw %>% group_by(Symb) %>% summarize_all(mean)
df$Period.Ending <- NULL
df$For.Year <- NULL
df <- merge(df, df_raw_sectInd, by = "Symb")
rm(df_raw); gc()
o <- c()
for(i in 1:ncol(df)){x <- df[,i]; o[i] <- length(which(is.na(x)))} 
o
rmcols <- which(o > 0)
colnames(df)[rmcols]
df <- df[, -rmcols]
rownames(df) <- df$Symb
df$Symb <- NULL
df$Sector <- as.factor(df$Sector)
df$Industry <- as.factor(df$Industry)
#colnames(df)


df_num <- df[, -(which(colnames(df) %in% c("Sector", "Industry")))]
zdf <- scale(df_num, center=TRUE, scale=TRUE)
#Is data clusterable? Hopkins statistic must be clost to 0 (much < 0.5).
get_clust_tendency(df_num, n = 50, gradient = list(low = "steelblue",  high = "white"))
#get_clust_tendency(zdf, n = 50, gradient = list(low = "steelblue",  high = "white"))
res <- PCA(df_num, ncp = 5)
fviz_screeplot(res, ncp=5)
fviz_pca_ind(res, axes = c(1,2), habillage = df[, "Sector"])#, addEllipses = T, ellipse.level=0.95)
df_outPCA <- data.frame(coord = res$ind$coord)
df_outPCApair <- df_outPCA[, c(1:4)]
ggpairs(df_outPCApair, aes(alpha = 0.4))

# fviz_pca_biplot(res,
#                 habillage = df[, "Sector"],
#                 addEllipses = TRUE,
#                 col.var = "red", alpha.var ="cos2",
#                 label = "var") +
#   scale_color_brewer(palette="Dark2")
#=====================

df_AdjP_Vol <- read.csv("prices-split-adjusted.csv", stringsAsFactors = F)
head(df_AdjP_Vol)
class(df_AdjP_Vol)
colnames(df_AdjP_Vol) <- c("Date", "Symb", "op", "cp", "low", "high", "vol")
df_AdjP_Vol$Date <- as.Date(df_AdjP_Vol$Date)
#df_AdjP_Vol <- merge(df_AdjP_Vol, df_raw_sectInd, by = "Symb")
df_cp <- df_AdjP_Vol[, c("Date", "Symb", "cp")]
df_vol <- df_AdjP_Vol[, c("Date", "Symb", "vol")]
#rm(df_AdjP_Vol); gc()
df_cp <- df_cp %>% spread(Symb, cp)
df_vol <- df_vol %>% spread(Symb, vol)
o <- c()
for(i in 1:ncol(df_cp)){x <- df_cp[, i]; o[i] <- length(which(is.na(x)))}
o
rmcols <- which(o >= 754)
df_cp <- df_cp[, -rmcols]
rmrows <- c(1:754)
df_cp <- df_cp[-rmrows, ]
# o <- c()
# for(i in 1:ncol(df_cp)){x <- df_cp[, i]; o[i] <- length(which(is.na(x)))}
# o
# rmcols <- which(o >= 80)
# df_cp <- df_cp[, -rmcols]
# df_cp <- na.spline(df_cp)
datevec <- df_cp$Date
rownames(df_cp) <- df_cp$Date
df_cp$Date <- NULL
df_cp <- as.xts(df_cp, order.by = datevec)
#---------------
#df_cp_w <- fortify(apply.weekly(df_cp, mean))
df_cp_w <- fortify(apply.monthly(df_cp, mean))
#df_cp_w <- fortify(apply.yearly(df_cp, mean))
library(lubridate)
df_cp_w$Month <- month(as.Date(df_cp_w$Index))
df_cp_w$Index <- NULL
zdf_cp_w <- scale(df_cp_w[, -ncol(df_cp_w)])
zdf_cp_w <- data.frame(Month = df_cp_w$Month, zdf_cp_w)
#df_cp_w[1:5,469:472]
#rownames(df_cp_w) <- df_cp_w$Index
zdf_cp_w <- zdf_cp_w %>% group_by(Month) %>% summarize_all(mean)
class(zdf_cp_w)
zdf_cp_w <- as.data.frame(zdf_cp_w)
rownames(zdf_cp_w) <- zdf_cp_w$Month
zdf_cp_w$Month <- NULL
zdf_cp_w <- as.data.frame(t(zdf_cp_w))
zdf_cp_w$Symb <- rownames(zdf_cp_w)
zdf_cp_w <- merge(zdf_cp_w, df_raw_sectInd, by = "Symb")
rownames(zdf_cp_w) <- zdf_cp_w$Symb
zdf_cp_w$Symb <- NULL
ind_num <- which(!(colnames(zdf_cp_w) %in% c("Sector", "Industry")))
#--
corMatMy <- cor(zdf_cp_w[, ind_num])
corrplot(corMatMy, order = "hclust")
#Apply correlation filter at 0.70,
highlyCor <- colnames(zdf_cp_w)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
print(highlyCor)
#--
res <- PCA(zdf_cp_w[, ind_num])
fviz_screeplot(res, ncp=5)
zdf_cp_w$Sector <- as.factor(zdf_cp_w$Sector)
fviz_pca_ind(res, habillage = zdf_cp_w[, "Sector"])
fviz_pca_biplot(res, habillage = zdf_cp_w[, "Sector"])#,







df_cp_w <- as.data.frame(t(df_cp_w))
df_cp_w$Symb <- rownames(df_cp_w)
df_cp_w <- merge(df_cp_w, df_raw_sectInd, by = "Symb")
rownames(df_cp_w) <- df_cp_w$Symb
df_cp_w$Symb <- NULL
ind_num <- which(!(colnames(df_cp_w) %in% c("Sector", "Industry")))
zdf_cp_w <- as.data.frame(t(zdf_cp_w))
zdf_cp_w$Symb <- rownames(zdf_cp_w)
zdf_cp_w <- merge(zdf_cp_w, df_raw_sectInd, by = "Symb")
rownames(zdf_cp_w) <- zdf_cp_w$Symb
zdf_cp_w$Symb <- NULL
res <- PCA(zdf_cp_w[, ind_num])
df_cp_w$Sector <- as.factor(df_cp_w$Sector)
zdf_cp_w$Sector <- as.factor(zdf_cp_w$Sector)
fviz_pca_ind(res, habillage = zdf_cp_w[, "Sector"])
fviz_pca_biplot(res, habillage = zdf_cp_w[, "Sector"])#,
#                addEllipses = TRUE,
  #               col.var = "red", alpha.var ="cos2",
  #               label = "var") +
  # scale_color_brewer(palette="Dark2")
fviz_pca_ind(res, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")#,
             #repel = TRUE # Avoid text overlapping (slow if many points)
)

# # Contributions of variables to PC1
# fviz_contrib(res, choice = "ind", axes = 1, top = 30)
# # Contributions of variables to PC2
# fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

xts_cpW <- apply.weekly(df_cp, mean)
class(xts_cpW)
#xts_cpW <- apply.monthly(df_cp, mean)
ema_per <- 21
xts_cpEMA <- xts(apply(xts_cpW, MARGIN = 2, FUN = "EMA", n = ema_per), index(xts_cpW))
class(xts_cpEMA)
#xts_cpEMA <- xts(apply(xts_cpW, MARGIN = 2, FUN = "EVWMA", volume = ), index(xts_cpW))
#num_na <- ema_per - 1
xts_cpEMA <- xts_cpEMA[c(ema_per:nrow(xts_cpEMA)),]
class(xts_cpEMA)
#volmat <- volmat[c(ema_per:nrow(volmat)),]; volmat <- round(volmat, 4)
sOsc_per <- 13
xts_cpsOsc <- xts(apply(xts_cpEMA, MARGIN = 2, FUN = "slopeOsc", wind = sOsc_per), index(xts_cpEMA))
xts_cpsOsc <- xts_cpsOsc[c((sOsc_per + 1):nrow(xts_cpsOsc)),]
class(xts_cpsOsc)
df_cpsOsc <- fortify(xts_cpsOsc)
rownames(df_cpsOsc) <- df_cpsOsc$Index
df_cpsOsc$Index <- NULL
df_cpO <- as.data.frame(t(df_cpsOsc))
df_cpO$Symb <- rownames(df_cpO)
df_cpO <- merge(df_cpO, df_raw_sectInd, by = "Symb")
rownames(df_cpO) <- df_cpO$Symb
df_cpO$Symb <- NULL
ind_num <- which(!(colnames(df_cpO) %in% c("Sector", "Industry")))
res <- PCA(df_cpO[, ind_num])
df_cpO$Sector <- as.factor(df_cpO$Sector)
fviz_pca_ind(res, habillage = df_cpO[, "Sector"])
fviz_pca_biplot(res, habillage = df_cpO[, "Sector"])
ggpairs(df_outPCA[, 1:5], aes(alpha = 0.4))
#--
df_cpO_sectr <- df_cpO %>% group_by(Sector) %>% summarize_all(mean)
df_cpO_sectr$Industry <- NULL
rownames(df_cpO_sectr) <- df_cpO_sectr$Sector
df_cpO_sectr$Sector <- NULL
df_cpO_sectr <- as.data.frame(t(df_cpO_sectr))
corMatMy <- cor(df_cpO_sectr)
corrplot(corMatMy, order = "hclust")
#Apply correlation filter at 0.70,
highlyCor <- colnames(df_cpO_sectr)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
print(highlyCor)
#------------------
class(xts_cpsOsc)
monthvec <- months(index(xts_cpsOsc))
df_cpsOsc <- data.frame(Month = monthvec, fortify(xts_cpsOsc))
df_cpsOsc$Index <- NULL
df_cpsOsc_month <- df_cpsOsc %>% group_by(Month) %>% summarize_all(mean)
rownames(df_cpsOsc_month) <- df_cpsOsc_month$Month
df_cpsOsc_month$Month <- NULL
df_cpsOsc_month <- as.data.frame(t(df_cpsOsc_month))
corMatMy <- cor(df_cpsOsc_month)
corrplot(corMatMy, order = "hclust")
#Apply correlation filter at 0.70,
highlyCor <- colnames(df_cpsOsc_month)[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
print(highlyCor)
#-----------------

#=============================================
datevec <<- index(xts_cpEMA)
indata <- xts_cpEMA[, "AAL"]
out <- getBuySellbinryVar(indata)
buy_sig <- as.integer(as.character(out[[1]]))
sell_sig <- as.integer(as.character(out[[2]]))
df_buysell <- data.frame(Date = index(xts_cpEMA), buySig = buy_sig, sellSig = sell_sig)
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


