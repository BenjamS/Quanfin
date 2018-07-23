setwd('D:/OneDrive - CGIAR/Documents')
#devtools::install_github("jessevent/crypto")
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
library(forecast)
library(stats)
library(crypto)
source('./collectiveModes.R')
#--------------------------------
df_raw <- crypto_history() #getCoins()
class(df_raw$date)
symbvec <- unique(df_raw$symbol)
namevec <- unique(df_raw$name)
n_names <- length(namevec)
# ind_top <- which(df_raw$ranknow <= 100)
# names_top <- unique(df_raw$name[ind_top])
# n_names_top <- length(names_top)
# df_raw_top <- subset(df_raw, name %in% names_top)
# df_cp <- df_raw_top[, c("date", "name", "close")]
# df_vol <- df_raw_top[, c("date", "name", "volume")]
#--------------------------------
start_dates <- c()
for(i in 1:n_names)
{
  this_name <- namevec[i]
  df_this <- subset(df_raw, name == this_name)
  this_start_date <- df_this$date[1]
  start_dates[i] <- this_start_date
}
# class(df_raw$date[this_ind[1]])
# as.Date(start_dates[1:3])
# as.numeric(start_dates[1:3]) < 16000
# sum(start_dates <= 15900)
cut_off_date <- as.Date("2017-01-01")
ind_keep <- which(as.Date(start_dates) < cut_off_date)
keep_these <- namevec[ind_keep]
#n_names <- length(keep_these)
#as.Date(start_dates[ind_keep])
df_subset <- subset(df_raw, name %in% keep_these)
df_subset <- df_subset[which(df_subset$date > cut_off_date), c("date", "name", "close")]
namevec <- unique(df_subset$name)
n_names <- length(namevec)
#--------------------------------
df_ts <- df_subset %>% spread(name, close)
date_vec <- df_ts$date
df_ts$date <- NULL
xts_ts <- xts(df_ts, date_vec)
o <- apply(xts_ts, 2, function(x) length(which(is.na(x))))
table(o)
rmcols <- which(o > 2)
colnames(xts_ts)[rmcols]
xts_cp_mat <- xts_cp_mat[, -rmcols]
#xts_vol_mat <- xts_vol_mat[, -rmcols]
o <- apply(xts_ts, 1, function(x) length(which(is.na(x))))
table(o)
xts_ts <- na.spline(xts_ts)
#xts_vol_mat <- na.spline(xts_vol_mat)
o <- apply(xts_ts, 1, function(x) length(which(is.na(x))))
table(o)
#--------------------------------
n_ts <- ncol(xts_ts)
#--------------------------------
datevec <- index(xts_ts)
mat_diff <- diff(scale(xts_ts))
mat_diff <- mat_diff[-1, ]
date_vec <- date_vec[-1]
out_collModes <- collectiveModes(mat_diff, datevec, df_group = NULL,
                                 Contrib_as_ModeSq = F,
                                 AggregateContributions = F)
#--------------------------------


















#If you have to remove duplicates use this:
# df_list <- list()
# for(i in 1:n_names)
# {
#   this_name <- namevec[i]
#   df_this <- subset(df_subset, name == this_name)
#   df_this <- df_this[, c("date", "close")]
#   colnames(df_this)[2] <- this_name
#   ind_rm <- which(duplicated(df_this$date)) - 1
#   df_list[[i]] <- df_this[-ind_rm,]
# }
# df_cp <- join_all(df_list)
# df <- df %>% group_by(name) %>% mutate(which(duplicated(date)) - 1)
#------
df_wide <- df_subset %>% spread(name, close)
df_wide <- subset(df_wide, date > cut_off_date)
datevec <- df_wide$date
rownames(df_wide) <- as.character(datevec)
df_wide$date <- NULL
o <- c()
for(i in 1:n_names){x <- df_wide[,i]; o[i] <- length(which(is.na(x)))} 
max(o)
table(o)
ind_rm <- which(o > 63)
names_with_lots_NA <- colnames(df_wide)[ind_rm]
n_names <- n_names - length(ind_rm)
df_wide <- df_wide[, -ind_rm]
df_wide <- na.spline(df_wide)

df_cp <- as.data.frame(df_wide)

ts <- df_cp$Litecoin
period <- sapply(df_cp, findfrequency)
print(period)

ts <- as.ts(df_cp$Diamond)
trend <- ma(ts, order = 30, centre = T)
plot(ts)
lines(trend)
plot(trend)
detrend <- ts - trend
plot(detrend)

#---
per <- 89
df_cp_ma <- as.data.frame(sapply(df_cp, ma, per))
df_cp_dt <- df_cp - df_cp_ma
df_cp_dt <- df_cp_dt[-c(1:per), ]
df_cp_ma <- df_cp_ma[-c(1:per), ]
datevec2 <- datevec[-c(1:per)]
per <- 55
df_cp_dt_smooth <- as.data.frame(sapply(df_cp_dt, ma, per))
df_cp_dt_smooth <- df_cp_dt_smooth[-c(1:per), ]
datevec3 <- datevec2[-c(1:per)]
df_cp_ma <- df_cp_ma[-c(1:per), ]
#---
selectnames <- c("Bitcoin", "Litecoin", "WorldCoin", "Deutsche.eMark", "Ripple", "HoboNickels")
df_mnth <- data.frame(date = datevec3, df_cp_dt_smooth)
df_mnth$Year <- year(df_mnth$date)
df_mnth <- df_mnth %>% gather(name, cp, Anoncoin:Zetacoin)
df_mnth$Month <- month(df_mnth$date)
df_mnth$Day <- yday(df_mnth$date)
df_mnth <- subset(df_mnth, name %in% selectnames)
df_mnth <- df_mnth %>% group_by(Year, name) %>% mutate(scale(cp))

df_mnth$Year <- as.factor(df_mnth$Year)
gg <- ggplot(df_mnth, aes(x = Day, y = `scale(cp)`, group = Year, color = Year)) + geom_line()
gg <- gg + facet_wrap(~name, ncol = 2)
gg

df_sub <- subset(df_mnth, name == "WorldCoin")
gg <- ggplot(df_sub, aes(x = Day, y = `scale(cp)`, group = Year, color = Year)) + geom_line()
gg


df_sub <- df_cp_ma[, "WorldCoin"]
df_sub <- data.frame(date = datevec3, zcp = df_sub)
df_sub$Year <- year(df_sub$date)
df_sub$Month <- month(df_sub$date)
df_sub$Day <- yday(df_sub$date)
df_sub <- df_sub %>% group_by(Year) %>% mutate(scale(cp))

df_sub$Year <- as.factor(df_sub$Year)
gg <- ggplot(df_sub, aes(x = Day, y = zcp, group = Year, color = Year)) + geom_line()
gg




df_per <- df_mnth %>% spread(cp, name)
period <- df

df_per <- df_mnth %>% group_by(name, Year) %>% mutate(sapply, findfrequency)


# df_mnth$name_year <- paste(df_mnth$name, df_mnth$Year)
# df_mnth$date <- NULL
# df_mnth$cp <- NULL
# df_mnth$name <- NULL
# df_mnth$Year <- NULL
# colnames(df_mnth)[3] <- "zcp_pa"
# head(df_mnth)
# df_mnth <- as.data.frame(df_mnth)
# df_mnth <- df_mnth %>% spread(Month, zcp_pa)
# class(df_mnth$zcp_pa)
# df_mnth$Year



monthsvec <- months(as.Date(df_cp$date))
#unique(monthsvec)
zdf_cp <- scale(df_cp)
# ts <- as.ts(zdf_cp$Bitcoin)
# stl(ts, s.window = "periodic")

zdf_cp <- data.frame(Month = monthsvec, zdf_cp)
zdf_cp <- zdf_cp %>% group_by(Month) %>% summarize_all(mean)
class(zdf_cp)
zdf_cp <- as.data.frame(zdf_cp)
rownames(zdf_cp) <- zdf_cp$Month
zdf_cp$Month <- NULL
zdf_cp <- as.data.frame(t(zdf_cp))
zdf_cp$name <- rownames(zdf_cp)
#--
ind_num <- c(1:n_names)
res <- PCA(zdf_cp[, ind_num])
fviz_screeplot(res, ncp=5)
#fviz_pca_ind(res, repel = T)
fviz_pca_biplot(res)
#--
corMatMy <- cor(zdf_cp[, ind_num])
corrplot(corMatMy, order = "hclust")
#--
zdf_cp$name <- NULL
res <- PCA(t(zdf_cp))
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res)
corMatMy <- cor(t(zdf_cp))
corrplot(corMatMy, order = "hclust")


scale_this <- function(x) as.vector(scale(x))
df_cp2 <- data.frame(date = datevec, df_cp)
df_cp2$Year <- year(df_cp2$date)
df_cp2$Month <- month(df_cp2$date)
df_cp2 <- df_cp2 %>% gather(name, cp, Anoncoin:Zetacoin)
df_cp2 <- df_cp2 %>% group_by(Year, name) %>% mutate(scale(cp))



df_cp_trunc <- data.frame(date = datevec, df_cp)
df_cp_trunc <- subset(df_cp_trunc, date < as.Date("2017-06-15") & date > as.Date("2015-01-01"))
datevec_trunc <- df_cp_trunc$date
df_cp_trunc$date <- NULL
zdf_cp <- scale(df_cp_trunc)
df_plot <- data.frame(date = datevec_trunc, zdf_cp)


df_outPCA <- data.frame(res$ind$contrib)
get_clust_tendency(zdf_cp, ind_num, n = 10, gradient = list(low = "steelblue",  high = "white"))
# get_clust_tendency(zdf_cp, n = 50, gradient = list(low = "steelblue",  high = "white"))
# get_clust_tendency(df_cp, n = 50, gradient = list(low = "steelblue",  high = "white"))
#Optimum number of clusters
fviz_nbclust(zdf_cp, kmeans, method = "gap_stat")



pam.res <- pam(zdf_cp, 4)
fviz_cluster(pam.res)
df_pam <- data.frame(Clust = pam.res$clustering)
df_pam <- cbind(df, df_pam)
df_pam$pamClust <- as.factor(df_pam$pamClust)
out_rf <- randomForest(zdf_cp)
# str(out_rf)
k = 2
MDSplot(out_rf, k)
k = 3
MDSplot(out_rf, zdf_cp, k)
out_rf.mds <- stats::cmdscale(1 - out_rf$proximity, eig = TRUE, k)
plot(out_rf.mds$points)
str(out_rf.mds)
class(out_rf.mds$points)
rownames(out_rf.mds$points)
df_rf <- as.data.frame(out_rf.mds$points)
df_rf <- data.frame(df_rf, name = rownames(out_rf.mds$points))
class(df_rf$name)
df_rf$name <- as.character(df_rf$name)
head(df_rf)
ggplot(df_rf, aes(x = V1, y = V2)) + geom_point() + geom_text(aes(label = name))
df_rf_g0 <- df_rf[df_rf$V1 > 0, ]
df_rf_l0 <- df_rf[df_rf$V1 < 0, ]







mc <- Mclust(zdf_cp)
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



































































u <- df_raw$name
df_list <- list()
for(i in 1:n_names_top)
{
  print(i)
  this_name <- namevec[i]
  this_ind <- which(u == this_name)
  df <- df_cp[this_ind, c("date", "close")]
  #df$close <- df_cp$close
  colnames(df)[2] <- this_name
  df_list[[i]] <- df
  print(head(df))
  rm(df)
}
df_CP <- join_all(df_list)
rm(df_list); gc()


#u <- which(df_raw$name == "Bitcoin")
df <- df_raw[, c("date", "close", "name")]
#colnames(df)
# datevec <- df$date
# df$date <- NULL
# xts_cp <- xts(df, datevec)

o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
max(o)
rmcols <- which(o > 2300)
colnames(xts_cp)[rmcols]
xts_cp <- xts_cp[, -rmcols]
xts_vol <- xts_vol[, -rmcols]
o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
max(o)
rmrows <- c(1:800)
xts_cp <- xts_cp[-rmrows,]
# xts_cp <- na.spline(xts_cp)
# xts_vol <- na.spline(xts_vol)
o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
o
rmcols <- which(o > 0)
colnames(xts_cp)[rmcols]
xts_cp <- xts_cp[, -rmcols]
dim(xts_cp)
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





df$Month <- months(as.Date(df$date))
zdf <- scale(df[, "close"])
zdf <- data.frame(Month = df$Month, close = zdf)
zdf <- zdf %>% group_by(Month) %>% summarize_all(mean)
class(zdf_cp)
zdf_cp <- as.data.frame(zdf_cp)
rownames(zdf_cp) <- zdf_cp$Month
zdf_cp$Month <- NULL
zdf_cp <- as.data.frame(t(zdf_cp))

class(zdf)


#fromdate<-"2012-01-12"; todate <- "2017-10-08"
















#df_cp <- fortify(apply.weekly(xts_cp, mean))
#df_cp_w <- fortify(apply.monthly(df_cp, mean))
#df_cp_w <- fortify(apply.yearly(df_cp, mean))

df_cp <- fortify(xts_cp)
#df_cp$Month <- months(as.Date(df_cp$Index))
df_cp$Index <- NULL
#unique(df_cp$Month)
zdf_cp <- scale(df_cp)
zdf_cp <- data.frame(Month = df_cp$Month, zdf_cp)
zdf_cp <- zdf_cp %>% group_by(Month) %>% summarize_all(mean)
class(zdf_cp)
zdf_cp <- as.data.frame(zdf_cp)
rownames(zdf_cp) <- zdf_cp$Month
zdf_cp$Month <- NULL
zdf_cp <- as.data.frame(t(zdf_cp))
zdf_cp$Symb <- rownames(zdf_cp)
#--
ind_num <- c(1:12)
res <- PCA(zdf_cp[, ind_num])
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res)
#--
corMatMy <- cor(zdf_cp[, ind_num])
corrplot(corMatMy, order = "hclust")
#--
df_cp$Month <- NULL
res <- PCA(t(zdf_cp))
fviz_screeplot(res, ncp=5)
fviz_pca_biplot(res)


df_outPCA <- data.frame(res$ind$contrib)
get_clust_tendency(zdf_cp, n = 50, gradient = list(low = "steelblue",  high = "white"))
# get_clust_tendency(zdf_cp, n = 50, gradient = list(low = "steelblue",  high = "white"))
# get_clust_tendency(df_cp, n = 50, gradient = list(low = "steelblue",  high = "white"))
#Optimum number of clusters
fviz_nbclust(t(df_cp), kmeans, method = "gap_stat")


pam.res <- pam(t(zdf_cp), 8)
fviz_cluster(pam.res)
df_pam <- data.frame(Clust = pam.res$clustering)
df_pam <- cbind(df, df_pam)
df_pam$pamClust <- as.factor(df_pam$pamClust)
out_rf <- randomForest(t(df_cp))
# str(out_rf)
k = 2
MDSplot(out_rf, k)

mc <- Mclust(t(zdf_cp))
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
































# pkg_vec <- c("tidyr", "FactoMineR", "factoextra", "GGally", "cluster", "mclust", "randomForest", "corrplot", "quantmod", "lubridate")
# n <- length(pkg_vec)
# for(i in 1:n)
# {
#   install.packages(pkg_vec[i])
# }



















#========================================================
library(crypto)

df_crypto <- getCoins()

symbvec <- unique(df_crypto$symbol)
namevec <- unique(df_crypto$name)
n <- length(namevec)
df_cp <- df_crypto[, c("date", "name", "close")]
df_vol <- df_crypto[, c("date", "name", "volume")]

u <- df_crypto$name
df_list <- list()
for(i in 1:n)
{
  print(i)
  this_name <- namevec[i]
  df <- df_cp[u == this_name, c("date", "close")]
  df$close <- as.numeric(df$close)
  colnames(df)[2] <- this_name
  df_list[[i]] <- df
}
df_out <- join_all(df_list)
rm(df_list); gc()

#fromdate<-"2012-01-12"; todate <- "2017-10-08"
datevec <- df_out$date
df_out$date <- NULL
xts_cp <- xts(df_out, datevec)

o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
max(o)
rmcols <- which(o > 2300)
colnames(xts_cp)[rmcols]
xts_cp <- xts_cp[, -rmcols]
xts_vol <- xts_vol[, -rmcols]
o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
max(o)
rmrows <- c(1:800)
xts_cp <- xts_cp[-rmrows,]
# xts_cp <- na.spline(xts_cp)
# xts_vol <- na.spline(xts_vol)
o <- c()
for(i in 1:ncol(xts_cp)){x <- xts_cp[,i]; o[i] <- length(which(is.na(x)))} 
o
rmcols <- which(o > 0)
colnames(xts_cp)[rmcols]
xts_cp <- xts_cp[, -rmcols]
dim(xts_cp)
#---------------
#df_cp <- fortify(apply.weekly(xts_cp, mean))
#df_cp_w <- fortify(apply.monthly(df_cp, mean))
#df_cp_w <- fortify(apply.yearly(df_cp, mean))





#================================
#================================
#================================
#================================
#================================
#================================

currpairs <- c("GBP=X", "AUS=X", "JPY=X", "CHF=X")
getSymbols("EUR=X",src="yahoo",from="20010-01-01")
