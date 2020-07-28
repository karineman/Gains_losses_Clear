# Functions accompaning 1_Data_Clean_GL: Analysis for Gains and Losses in Clear Creek

#---- Separate breakthrough curves from each round
# using time to peak plus 2 times time to peak
# currently does not separate the curves as needed
BTC_separation<- function (df=dis_meta_nest$data, df2=COND){
  data1<-inner_join(df2, df, by = c('date', 'site')) %>% 
  # convert to time since injection
  mutate(time_since_injection = as.numeric(difftime(datetime, inj_time, tz= 'CST6CDT',
                                         units = "min"))) %>% 
  filter(datetime>=rise_time & datetime<=end_time)
  peak_time<-data1$time_since_injection[which.max(data1$SC_uS)]
  twotime_peak<-peak_time+2*peak_time
  data<-inner_join(df2, df, by = c('date', 'site')) %>% 
  mutate(time_since_injection = difftime(datetime, inj_time, tz= 'CST6CDT',
                                           units = "min")) %>% 
  filter(datetime>=rise_time & time_since_injection<=twotime_peak)
  return(data)
}

#---- Separate breakthrough curves from each round
# using start and end time from visual inspection of Conductivity data
# currently separate the curves as needed
# base_time is time added to determine baseline concentrations
BTC_separation_set<- function (df=dis_meta_nest$data, df2=COND, base_time=10){
  data1<-inner_join(df2, df, by = c('date', 'site')) %>% 
    # convert to time since injection
    mutate(time_since_injection = as.numeric(difftime(datetime, inj_time, tz= 'CST6CDT',
                                           units = "min"))) %>% 
    mutate(begin=rise_time-base_time*60) %>% 
    mutate(end=end_time-base_time*60) %>% 
    filter(datetime>=begin & datetime<=end_time) %>% 
    mutate(time_to_peak=time_since_injection[which.max(SC_uS)])
  # twotime_peak<-peak_time+2*peak_time
  # data<-inner_join(df2, df, by = c('date', 'site')) %>% 
  #   mutate(time_since_injection = difftime(datetime, inj_time, tz= 'CST6CDT',
  #                                          units = "min")) %>% 
  #   filter(datetime>=rise_time & time_since_injection<=twotime_peak)
  return(data1)
}

#---- Determine baseline (same as used int BTC_Functions as part of the Conductivity_BTC_Clear Analysis)
# from time of injection through the baseline window (which is set by user)
# baseline window default is 120 dt (i.e.10 minutes if your dt is 5 seconds)
baseline <-
  function(df = BTC_GL_nest$data) {
    #BASELINE START
    # convert uS/cm to g/L NaCl
    df$NaCl = 0.51 * df$SC_uS / 1000
    # determine index of inject time
    start.inj.idx <- which(df$datetime == df$rise_time[1])
    # determine index of inject time plus baseline window
    baseline.start.idx <- which(df$datetime == df$begin[1])
    # determine median SC_uS during the window defined before the BTC
    startbase <-
      median(df$NaCl[baseline.start.idx:start.inj.idx], na.rm = TRUE)
    # apply start background correction subtract background conductivity value from data
    df$NaCl_BGS <- pmax(df$NaCl - startbase, 0)
    #BASELINE END
    #determine index
    end.inj.idx <- which(df$datetime == df$end_time[1])
    # determine index of inject time plus baseline window
    baseline.end.idx <- which(df$datetime == df$end[1])
    # determine median SC_uS during the window defined after the BTC
    endbase <-
      median(df$NaCl[baseline.end.idx:end.inj.idx], na.rm = TRUE)
    # apply end background correction subtract background conductivity value from data
    df$NaCl_BGE <- pmax(df$NaCl - endbase, 0)
    #BASELINELINEAR
    # define baseline start/end time and value
    basex <- c(df$rise_time[1], df$end_time[1])
    basey <- c(startbase, endbase)
    #create linear model between points
    BaseFit <- lm(basey ~ basex)
    # use model to predict linear baseline conductivity
    newdata <- data.frame(basex = df$datetime)
    baselineLine <- predict(BaseFit, newdata)
    # apply linear background correction subtract background conductivity value from data
    df$NaCl_BGL <- pmax(df$NaCl - baselineLine, 0)
    return(df)
  }

#---- Determine discharge (same as used int BTC_Functions as part of the Conductivity_BTC_Clear Analysis)
# integrate the background corrected NaCl breakthrough curve
# find Q with known mass NaCl assuming all the mass was recovered
# set dt to match the timestep of your data
discharge<-function(df=BTC_dis_nest$data){
  dt=as.numeric(df$datetime[2]-df$datetime[1])
  output <- df %>%
    #Integration
    mutate(IntC_s1_L = NaCl_BGL + lag(NaCl_BGL, default = 0)) %>%
    mutate(IntC_L = cumsum(IntC_s1_L * .5 * dt)) %>% 
    mutate(C_D = max(IntC_L)) %>% 
    # mutate(per_rec_L = IntC_L / max(IntC_L)) %>%
    mutate(Q_D = mass_NaCl[1]/C_D[1]) %>% 
    mutate(M_D= mass_NaCl)
    # halfRecovery.idx<- which.min(abs(output$per_rec_L-.50)) 
    # output2<-output %>% 
    # mutate(time_to_half = time_since_injection[halfRecovery.idx]*60) %>% 
    # mutate(Med_Vel = length/time_to_half)
  return(output)
}

#---- Determine Discharge in upstream reach station

upper_discharge<-function(df=Dis_sum_nest$data){
  output<-df %>% 
    mutate(M_U= lag(M_D, default = NA)) %>% 
    mutate(C_U= lag(C_D, default = NA)) %>% 
    mutate(Q_U= lag(Q_D, default = NA)) 
  return(output)
}

#---- Determine change in discharge over the reach
# 
delta_discharge<-function(df=delta_Q_nest$data){
  Q_1<- df$Q[1]
  Q_2<- df$Q[2]
  Q_3<- df$Q[3]
  Q_4<- df$Q[4]
output<- df %>% 
  #Change of Q from station 1 (in this case site 0)
  mutate(delta_Q_1= Q-Q_1) %>% 
  #Change of Q from station 2 (in this case site 150)
  mutate(delta_Q_2= ifelse(site>150, Q-Q_2, NA)) %>% 
  #Change of Q from station 3 (in this case site 350)
  mutate(delta_Q_3= ifelse(site>350, Q-Q_3, NA)) %>% 
  #Change of Q from station 4 (in this case site 600)
  mutate(delta_Q_4= ifelse(site>600, Q-Q_4, NA)) 
  return(output)
}


#---- Determine Mass Recovery (same as used in BTC_Functions as part of the Conductivity_BTC_Clear Analysis)
# integrate the background corrected NaCl breakthrough curve
# calculate mass recovery with Q calculated above for each station on specific dates
# set dt to match the timestep of your data

recovery<-function(df=BTC_recov_nest$data){
  dt=as.numeric(df$datetime[2]-df$datetime[1])
  output <- df %>%
    #Integration
    mutate(IntC_s1_L = NaCl_BGL + lag(NaCl_BGL, default = 0)) %>%
    mutate(IntC_L = cumsum(IntC_s1_L * .5 * dt)) %>% 
    mutate(C_UD = max(IntC_L)) %>% 
    mutate(Mrec =  Q_D*C_UD[1]) %>% 
    mutate(Mloss = Mrec-M_U) %>% 
    mutate(per_Mloss=Mloss/M_U) %>% 
    mutate(delta_Q = Q_D-Q_U) %>% 
    mutate(Qloss_min = Mloss/C_U) %>% 
    mutate(Qloss_max = Mloss/C_UD[1]) %>% 
    mutate(Qgain_min = delta_Q-Qloss_min) %>% 
    mutate(Qgain_max = delta_Q-Qloss_max)
  return(output)
}
