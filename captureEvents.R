captureEvents<-function(indata,wind=13,graph1=1,graph2=1)
{
  datevec<-as.numeric(gsub('-', '', datevec))
  datevec<-as.Date(as.character(datevec),format="%Y%m%d")
  #-----------------------
  #   ema1<-EMA(indata,3)
  #   indata<-ema1
  out <- slopeOsc(indata,wind,graph1, getbinry = T)
  slope<-out[[1]]
  mu_slope<-out[[2]]
  #-----------------------
  #one way to do it...
  ind_ut<-which(slope>mu_slope)-round(wind/2)
  ind2_ut<-which(diff(ind_ut)>1)
  ind2_ut<-c(1,ind2_ut,length(ind_ut))
  #-----------------------
  #another way...slope of slope (2nd derivative of EMA)
  #   out2<-slopeOsc(slope,wind=3,graph=0);#out2<<-out2
  #   slope2<-out2[[1]]
  #   mu_slope2<-out2[[2]]
  #   ind_ut<-which(slope2>mu_slope2)
  #   ind2_ut<-which(diff(ind_ut)>1)
  #   ind2_ut<-c(1,ind2_ut,length(ind_ut))
  #-----------------------
  price_range<-max(indata,na.rm=TRUE)-min(indata,na.rm=TRUE)
  #print(price_range)
  n_ut<-0
  indxlist<-list()
  magvec<-c()
  magvecNormd <- c()
  magvecNormd2 <- c()
  #roughvec<-c()
  for(i in 2:length(ind2_ut))
  {
    #  print(i)
    ind<-c()
    ind<-ind_ut[(ind2_ut[i-1]+1):ind2_ut[i]]
    if(length(ind)==1){print("hello"); next}
    if(ind[1]<8)
    {
      ind<-c(1:(ind[length(ind)]+8))
    }
    else if((ind[length(ind)]+8)>length(indata))
    {
      ind<-c((ind[1]-8):length(indata))
    }
    else
    {    
      ind<-c((ind[1]-8):(ind[length(ind)]+8))
    }
    u_beg<-which(indata[ind]==min(indata[ind],na.rm=TRUE))
    u_end<-which(indata[ind]==max(indata[ind],na.rm=TRUE))
    ind_beg<-ind[u_beg]
    ind_end<-ind[u_end]
    #  print(paste(ind_beg,ind_end))
    if(ind_beg>ind_end){next}
    min_thisRange <- min(indata[ind],na.rm=TRUE)
    this_movemnt<-max(indata[ind],na.rm=TRUE)-min_thisRange
    #  print(this_movemnt)
    if(this_movemnt<.1*price_range){next}
    #--Passed
    n_ut<-n_ut+1
    indxlist[[n_ut]]<-c(ind_beg:ind_end)
    magvec[n_ut]<-this_movemnt
    frac_increase <- this_movemnt / min_thisRange
    magvecNormd[n_ut] <- frac_increase
    this_slope <- this_movemnt / length(ind)
    magvecNormd2[n_ut] <- frac_increase * this_slope
    #roughvec[n_ut]<-sd(indata[ind])/mean(indata[ind])
  }
  print(paste("Num. Up-trends: ",n_ut))
  #-----------------------
  if(graph2==1)
  {
    #     windows()
    par(mar=c(3,2.5,1,1.5));par(mfrow=c(2,1))
    #indata<-zoo(indata)
    indata<-xts(indata,datevec)
    indata<-zoo(indata)
    plot(indata,type="l",main="Data with trends")
    #lines(indata,col="coral")
    for(i in 1:length(indxlist))
    {
      ind<-c()
      ind<-indxlist[[i]]
      lines(indata[ind],col="green",lwd=2)
    }
    #--
    slope<-as.zoo(slope)
    plot(slope,type="l",main="EMA Slope Osc.")
    abline(h=mu_slope,col="green")
    #--
    #     slope2<-as.zoo(slope2)
    #     plot(slope2,type="l",main="2nd Order EMA Slope Osc.")
    #     abline(h=mu_slope2,col="green")
    
  }
  #return(list(indxlist,magvec,roughvec))
  return(list(indxlist,magvec, magvecNormd, magvecNormd2))
}