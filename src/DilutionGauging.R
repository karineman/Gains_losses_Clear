#OLD BASELINE START
baseline_start<-function(raw = rawSalt, choose.baseline.window=30, rise = riseTime){
  # determine index of inject time
  start.inj.idx<- which(raw$datetime==rise)
  # baseline calculation - baseline window define time peroid that you are pulling the median from
  # baseline window default is 30 dt (i.e.30 minutes if your dt is 60 seconds)
  baseline.start.idx<- ifelse(start.inj.idx-choose.baseline.window>0, start.inj.idx-choose.baseline.window,0)
  baseline <- median(raw$Cl[baseline.start.idx:(start.inj.idx-1)])
  output <- baseline
  return(output)
}

#OLD BASELINE END
baseline_end<-function(raw = rawSalt, choose.baseline.window=30, baseend = endTime){
  # determine index of end of BTC time
  end.inj.idx<-which(raw$datetime==baseend)
  # baseline calculation - baseline window define time peroid that you are pulling the median from
  # baseline window default is 30 dt (i.e.30 minutes if your dt is 60 seconds)
  # baseend has added 30 minutes for processing ****LOOK IN CONDBTC if baseend does not add minutes this won't work right***
  baseline.end.idx<- ifelse(end.inj.idx-choose.baseline.window>0, end.inj.idx-choose.baseline.window,0)
  baseline <- median(raw$Cl[baseline.end.idx:(end.inj.idx-1)])

  output <- baseline
  return(output)
}
# OLD BASELINE LINEAR
baseline_linear<-function(line = adjSalt, rtime = riseTime, etime = endTime, bstart = baselineStart, bend = baselineEnd){
  # define baseline start/end time and value
  basex<- c(rtime, etime)
  basey<- c(bstart, bend)
  #create linear model between points
  BaseFit <- lm(basey~basex)
  # predict new concentrations for each date time
  newdata <- data.frame(basex = line$datetime)
  baselineLine <- predict(BaseFit, newdata)
  # create dataframe with info
  baselineLinear <-data.frame(line$datetime, baselineLine)
  colnames(baselineLinear)<-c('datetime','BaseConc')

  output <- baselineLinear
  return(output)
}

#determine time breakthrough begins - does this by finding the start of the longest postive run in slopes after injection - not a full proof solution
breakthroughStart <- function(z=i){
  slope <- c(0,diff(adjSalt$Cl_rollavg[adjSalt$datetime>=inject_datetime & adjSalt$datetime<=peak.datetime]))
  if(length(slope) <2){
    error[z] <- 'injection/peak time error' 
    next
  }
  slope_dir = ifelse(slope>0,1,0)
  runs<- rle(slope_dir)
  
  #index of rle vector where positive runs occur and determine position of my run
  myruns <- which(runs$values == 1 & runs$length >1)
  myruns.max <- which.max(runs$lengths[myruns])
  #position of where runs end
  runs.lengths.cumsum <- cumsum(runs$lengths)
  ends <-runs.lengths.cumsum[myruns]
  #new index is rle position previous to myruns
  newindex <- ifelse(myruns>1, myruns-1, 0)
  #starting positions of my runs
  starts <- runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  
  #index of max run start
  start.max <- starts[myruns.max]
  peakStart <- adjSalt$datetime[start.max]
  return(peakStart)
}

#Integrate under the curve 
#data is the dataframe with parameter timestep. Para corrected parameter concentration
BTCint<-function(data = filterSalt, para = filterSalt$CorrCl, dt = dt){
  intCl<- mutate(data,IntC_s1 = para+lag(para,default=0))%>%
    mutate(IntC = cumsum(IntC_s1*.5*dt))%>%
    mutate(recovery = IntC/max(IntC))
  return(intCl)
}
