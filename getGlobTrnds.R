getGlobTrnds<-function(fromdate, todate)
{  
  #Get Global Trends----------------------
  #Notes:
  #Have to load this package first: library(quantmod)
  #library(ggplot2)
  #library(tidyr)
  #To get symbols of all stocks on NYSE, NASDAQ, and AMEX use (eg. NYSE): stockSymbols("NYSE")$Symbol
  #!!!!!1)MAKE SURE TO REMOVE EXISTING "gtIDmat" BEFORE RUNNING.!!!!!
  #Eg.: fromdate<-"2012-01-12"; todate <- "2017-07-30"; getGlobTrnds(fromdate, todate)
  #2)If going to include ^VIX, must also include VXX in order to get volume proxy
  #3)fromdate can be string of the format "YYYY-MM-DD" or number of format YYYYMMDD
  #(If string, this program will convert the string to YYYYMMDD)  
  #4) TO get VVIX see getVVIX()
  #--Get ETFs representing:
  #1) Major Int Mkts (S&P 500, Nikkei, Dax, etc.)
  #2) main US mkt sectors
  #3) commodities
  #4) currencies
  #5) emerging mkts
  #--------------------------------------------  
  #setwd("C:/Users/Ben/Documents/R/myRfiles") -- my old computer
  #setwd("C:/Users/bensc/OneDrive/Documents")
  #setwd("D:/OneDrive - CGIAR/Documents")
  #--------------------------------------------  
  if(exists("gtIDmat")==TRUE){rm(gtIDmat,pos=".GlobalEnv")}
  print("Reading in golbal trend ID data.  Could take a moment...")
  #gtIDmat<-read.csv("GlobTrndsID.csv",as.is=TRUE)
  gtIDmat<-read.csv("GlobTrndsID.csv", stringsAsFactors = FALSE)
  gtIDmat<-gtIDmat[,2:ncol(gtIDmat)];gtIDmat<<-gtIDmat
  if(exists("GenSpecIDmat")==TRUE){print(exists("GenSpecIDmat"));print(dim(GenSpecIDmat));rm(GenSpecIDmat,pos=".GlobalEnv")}
  if(exists("cpgtetfmat")==TRUE){print(exists("cpgtetfmat"));print(dim(cpgtetfmat));rm(cpgtetfmat,pos=".GlobalEnv")}
  if(exists("lcpgtetfmat")==TRUE){print(exists("lcpgtetfmat"));print(dim(lcpgtetfmat));rm(lcpgtetfmat,pos=".GlobalEnv")}
  if(exists("volgtetfmat")==TRUE){rm(volgtetfmat,pos=".GlobalEnv")}
  if(exists("lvolgtetfmat")==TRUE){rm(lvolgtetfmat,pos=".GlobalEnv")}
  if(exists("datevec")==TRUE){rm(datevec,pos=".GlobalEnv")}
  #--------------------------------------------
  # if(class(fromdate)=="character")
  # {
  #   fromdate_str<-fromdate
  #   fromdate<-as.numeric(gsub('-', '', fromdate))
  # }
  #--------------------------------------------
  Nrowinstkmat <- nrow(gtIDmat)
  #--------------------------------------------  
  # cool progress bar to see the % of completion
  #pb <- txtProgressBar(min = 0, max = Nrowinstkmat, style=3)
  #--------------------------------------------
  cpgtetflist<-list();volgtetflist<-list()
  gtetfdone<-c();Gendone<-c();Subdone<-c();Specdone<-c()
  removed<-c();Nrem<-0
  t<-0;Nprocsd<-0;
  
  remove_vec <- c("UNG", "JJT", "LEDD", "FORX", "FUDD", "ICI", "NGE", "THD")
  
  for(i in 1:Nrowinstkmat)
  {
    Nprocsd<-Nprocsd+1; #Num processed (total, with or without download error)
    gtetf<-gtIDmat[i,1];Gen_type<-gtIDmat[i,2];Sub_type<-gtIDmat[i,3];Spec_track<-gtIDmat[i,4]
    
    if(gtetf %in% remove_vec){Nrem<-Nrem+1;removed[Nrem]<-gtetf;print(paste("Removed: ",gtetf)); next()}
    
    tryit <- try(getSymbols(gtetf, from = fromdate, to = todate, src = "yahoo", adjust = TRUE))
    if(inherits(tryit, "try-error"))# only true if an error occurs
    {Nrem<-Nrem+1;removed[Nrem]<-gtetf;print(paste("Removed: ",gtetf))
    next()
    }else{
      t<-t+1 #Num. retrieved (didn't encounter error)
      xts_ohlcv <- get(tryit)
      ucpvec <- xts_ohlcv[, 4]
      uvolvec <- xts_ohlcv[, 5]
      # ucpvec <- zoo(ucpvec)
      # uvolvec<-zoo(uvolvec)
      datevecstr<-as.Date(rownames(zoo(ucpvec)))
      print(head(datevecstr))
      
      if(t==1)
      {
        ucpmat <- ucpvec
        uvolmat <- uvolvec
      }
      else
      {
        #--
        ucpmat<-merge.xts(ucpmat,ucpvec)
        #--
        uvolmat<-merge.xts(uvolmat,uvolvec)
        #--
      }
      #----------------------
      #Plot cp and vol to visually check
      df_plot <- fortify(ucpvec)
      colnames(df_plot) <- c("Date", "CP Adj")
      #colnames(df_plot) <- c("Date", "OP", "HP", "LP", "CP", "Vol", "CP Adj")
      #df_plot <- gather(df_plot, Item, Value, OP:`CP Adj`)
      #df_plot <- subset(df_plot, Item %in% c("CP Adj", "Vol"))
      gg <- ggplot(df_plot, aes(x = Date, y = `CP Adj`)) + geom_line()
      gg <- gg + ggtitle(gtetf)
      print(gg)
      Sys.sleep(0.75)
      # par(mar=c(3,2.5,1,1.5));par(mfrow=c(2,1))
      # plot(ucpvec,type="l",main=paste("CP",gtetf))
      # plot(uvolvec,type="l",main=paste("Vol",gtetf))
      #----------------------
      #      cpgtetflist[[gtetf]]<-cpgtetf
      #      volgtetflist[[gtetf]]<-volgtetf
      gtetfdone[t]<-gtetf;Gendone[t]<-Gen_type;Subdone[t]<-Sub_type;Specdone[t]<-Spec_track
      #--PRINT
      print(paste("Global Trend: ",gtetf, " Gen. Type: ",Gen_type," Subtype: ",Sub_type," Specifically Tracks: ",Spec_track))
      print(paste("Num. processed so far: ",Nprocsd))
      print(paste("Num. passed so far: ",t))
      print(paste("Num. removed: ",Nrem))
      print(paste("Num. remaining: ",(Nrowinstkmat-Nprocsd)))
      print(paste("Percent complete: ",(round(100*Nprocsd/(Nrowinstkmat-Nrem), 1))))
    }
  }
  #cpgtetfmat<-do.call(cbind,cpgtetflist)
  #cpgtetfmat<-replacena(cpgtetfmat)
  cpgtetfmat<-ucpmat
  colnames(cpgtetfmat)<-gtetfdone
  #volgtetfmat<-do.call(cbind,volgtetflist)
  #volgtetfmat<-replacena(volgtetfmat)
  volgtetfmat<-uvolmat
  colnames(volgtetfmat)<-gtetfdone
  #gtetfdone<<-gtetfdone
  GenSpecIDmat<-cbind(gtetfdone,Gendone,Subdone,Specdone)
  colnames(GenSpecIDmat)<-c("Global Trend ETF","General Type","Sub-Type","Specific Track")
  #--------------------------------------
  if(is.element("^VIX",GenSpecIDmat[,1])==TRUE)
  {
    ind_VXX<-which(GenSpecIDmat[,1]=="VXX")
    ind_VIX<-which(GenSpecIDmat[,1]=="^VIX")
    volgtetfmat[,ind_VIX]<-volgtetfmat[,ind_VXX]  
  }
  #--------------------------------------
  #setwd("C:/Users/Ben/Documents/R/myRfiles")
  #setwd("C:/Users/bensc/OneDrive/Documents")
  write.zoo(as.xts(cpgtetfmat),"cpgtetf.csv",sep=",")
  write.zoo(as.xts(volgtetfmat),"volgtetf.csv",sep=",")
  #outmat<-list(cpgtetfmat,volgtetfmat,GenSpecIDmat)
  cpgtetfmat <<- cpgtetfmat; volgtetfmat <<- volgtetfmat; GenSpecIDmat <<- GenSpecIDmat
  datevec <<- time(cpgtetfmat)
  #  if(class(datevec)=="character")
  #  {
  #    datevec<-as.numeric(gsub('-', '', datevec))
  #    datevec<<-as.Date(as.character(datevec),format="%Y%m%d")
  #    fromdate_str<-fromdate
  #    fromdate<-as.numeric(gsub('-', '', fromdate))
  #  }
  #return(outmat)
  #setTxtProgressBar(pb, i)
}