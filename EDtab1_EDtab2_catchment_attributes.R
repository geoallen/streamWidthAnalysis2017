# EDtable1_2_catchment_attributes.R
# George Allen, Sept 2015

# For each stream network, adds gauge and topographic
# information to table 2:

EDtab1_EDtab2_catchment_attributes <- function(inTabPaths, tabNames, csvOut, workingDir) {
  
  require(MASS)
  
  fieldDates = rbind(# physiographically diverse surveys:
                     c("2015-05-25", "2015-05-26", ''),
                     c("2015-05-29","2015-05-30",""),
                     c("2015-06-01","2015-06-02",""),
                     c("6/20/2015","",""),
                     c("28/06/15","29/06/15",""),
                     c("7/12/2015","7/13/2015","7/14/2015"),
                     c("09/27/2015","",""),
                     # repeat surveys: 
                     c("2015-10-27","",""),
                     c("2015-12-09","",""),
                     c("2016-02-02","",""),
                     c("2016-02-14","",""),
                     c("2016-03-04","",""),
                     c("2016-03-04","",""))
  
  # gauge upstream basin area:
  gauge_A_km2 = c(10.59, 27.2, 16.8, 104, 6350, 64, 1.3, rep(0.48, 6))
  
  # get paths of discharge records: 
  qTabPaths = paste0(sapply(strsplit(inTabPaths,'/field/'), '[[' ,1), 
                     '/discharge_records/daily')
  # type of gauge record:
  gaugeType = c('USGS', 'USGS', 'USGS', 'CPCRW', 'NZhourly', 'NZdaily', 'stony',
                rep('stony_sub', 6))
  # size of catchments (ha):
  basinA = c(1.8, 4.49, 3.62, 2.58, .06, 0.56, 1.28, rep(0.484, 6))
 
  # set up out table:
  tabHdr = c("wName", 
             "Date.Surveyed", 
             "gauge_A_km2",
             "QRecLength_yrs", 
             "Q_cms", 
             "runoff_mmPerDay",
             "Q_percentile", 
             "N_Width_Obs",	
             "ADN_length", 
             "drainageDensity", 
             "basinSA", 
             "lnMode", 
             "lnR2", 
             "kernMode", 
             "mean1OrdStrWidth_cm", 
             "med1OrdStrWidth_cm", 
             "mode1OrdStrWidth_cm",
             "ADN_relief")
  oTab = data.frame(array(NA, c(length(tabNames), length(tabHdr))))
  names(oTab) = tabHdr
  oTab$wName = tabNames
  
  oTab$gauge_A_km2 = gauge_A_km2
  
  # band width and type of kernel density estimation:
  bw = 10
  kernel = "gaussian"
  
  # histogram bin interval:
  int = 100
  
  # define variable for the mean width of all first order streams:
  fOw_all = vector()
  
  # for each catchment, generate table of statistics:
  for (i in 1:length(inTabPaths)){
    print(i)
    
    # add survey date to table:
    oTab$Date.Surveyed[i] = paste(range(fieldDates[i, fieldDates[i,] != ""]), collapse=" - ")
    
    ###############################
    ## Analyze discharge data:
    
    # option to skip the gauge data analysis calculation to save time:
    io = 1; if (io == 0){
      oTab$Q[i] = "skipped"
      oTab$QRecLength_yrs[i] =  "skipped"
    }else{
      
      # analyze USGS gauge records:
      if (gaugeType[i] == 'USGS'){ oTab = USGSgauge(qTabPaths, fieldDates, oTab, i)}  
      # analyze Caribou gauge records:
      if (gaugeType[i] == 'CPCRW'){ oTab = CPCRW(qTabPaths, fieldDates, oTab, i)}
      # analyze NZ hourly records from NZ:
      if (gaugeType[i] == 'NZhourly'){ oTab = NZhourly(qTabPaths, fieldDates, oTab, i)}
      # analyze NZ daily records:
      if (gaugeType[i] == 'NZdaily'){oTab = NZdaily(qTabPaths, fieldDates, oTab, i)}
      # analyze stony 5 min interval lower gauge records:
      if (gaugeType[i] == 'stony'){ oTab = stony(qTabPaths, fieldDates, oTab, i)}
      # analyze stony subcatchment 5 min interval upper gauge records:
      if (gaugeType[i] == 'stony_sub'){ oTab = stony_sub(qTabPaths, fieldDates, oTab, i)}
    }
    
    
    
    ###############################
    ## field data:
    tab = data.frame(read.csv(inTabPaths[i], skip=49))
    
    tab$percent_nonflow[is.na(tab$percent_nonflow)] = 0
    w_raw = tab$flowing_width*(1-tab$percent_nonflow/100)
    w_raw[grep('RF', tab$code)] = w_raw[grep('RF', tab$code)]*39.3701 # Rangefinder Convert
    notZero = w_raw!=0
    w = w_raw[notZero]
    w = w * 2.54 # inches to cm convert
    
    N = length(w)
    L = (N*5*1e-3)
    DD = L/basinA[i]
    SA = 100*(sum(w*1e-2)*5*1e-6)/basinA[i]
    
    # fit a log normal frequency distribution using MLE:
    dlnInt = 0.1
    dlnSeq = seq(0, max(w), dlnInt)
    
    fit = fitdistr(w, "log-normal")$estimate
    dln = dlnorm(dlnSeq, fit[1], fit[2])
    lnMode = dlnSeq[which.max(dln)]
    
    breaks = seq(0, floor(max(w)+int), int)
    h = hist(w, breaks=breaks, plot=F)
    
    Ei = dlnorm(h$mids, fit[1], fit[2])*sum(w)*int*max(dln)
    lnR2 = r2(h$counts, Ei)
    
    
    # bin width data:
    int = 10
    breaks = seq(0, floor(max(w)+int), int)
    h = hist(w, breaks=breaks, plot=F)
    
    # find mode of data with gaussian density kernel:
    # charactarize distribution with a density kernel:
    kern = density(w, bw=bw, kernel=kernel, na.rm=T)
    kernMode = kern$x[which.max(kern$y)]
    
    
    # find first order mean stream width:
    fOw_raw = w_raw[tab$stream_order == 1]
    notZero = (fOw_raw!=0 & !is.na(fOw_raw))
    fOw = fOw_raw[notZero]*2.54# inches to cm convert
    fOw_all = c(fOw_all, fOw)
    mFoW = mean(fOw, na.rm=T) 
    medFoW = median(fOw, na.rm=T)
    modeFoW = density(fOw, bw=bw, kernel=kernel, na.rm=T)
    kernModeFoW = modeFoW$x[which.max(modeFoW$y)]
    
    oTab$N_Width_Obs[i] = N
    oTab$ADN_length[i] = round(L, 3)
    oTab$drainageDensity[i] = round(DD, 3)
    oTab$basinSA[i] = round(SA, 3)
    oTab$lnMode[i] = round(lnMode, 3)
    oTab$lnR2[i] = round(lnR2, 3)
    oTab$kernMode[i] = round(kernMode, 3)
    oTab$mean1OrdStrWidth_cm[i] = round(mFoW, 3)
    oTab$med1OrdStrWidth_cm[i] = round(medFoW, 3)
    oTab$mode1OrdStrWidth_cm[i] = round(kernModeFoW, 3)
    oTab$ADN_relief[i] = round(diff(range(tab$elev_m, na.rm=T)))
  }
  
  print(paste0("mean first order width: ", round(mean(fOw_all),1),'±', round(sd(fOw_all),1)))
  print(paste0("median first order width: ", round(median(fOw_all),1),
               ' +', quantile(fOw_all)[[2]],
               ' -', quantile(fOw_all)[[4]]))
  
  modeFoW = density(fOw_all, bw=bw, kernel=kernel, na.rm=T)
  kernModeFoW = modeFoW$x[which.max(modeFoW$y)]
  print(paste0("mode first order width: ", round(kernModeFoW,1), " cm"))
  
  print(mean(oTab$ADN_relief))
  oTab = as.data.frame(t(oTab))
  colnames(oTab) = tabNames
  oTab = oTab[-1,]
  
  write.csv(oTab, csvOut)
  cmd = paste('open', csvOut)
  system(cmd)
  
  
  #####################################################
  # functions:
  
  # R2 statistical test:
  r2 <- function(Oi, Ei){ 1-(sum((Oi-Ei)^2))/(sum((Oi-mean(Oi))^2)) }
  
  # analyze gauge data: 
  USGSgauge <- function(qTabPaths, fieldDates, oTab, i){
    # analyze USGS gauge records:
    qTab = read.table(paste0(qTabPaths[i], '.txt'), sep="\t", fill=T, skip=14, header=T)[-1,]
    q = as.numeric(as.vector(qTab[, 4]))*0.0283 # cfs to cms convert
    
    fieldQ = q[match(fieldDates[i,], qTab$datetime)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
  }
  
  CPCRW <- function(qTabPaths, fieldDates, oTab, i){
    # analyze caribou-poker creek gauge data (data from Jay Jones):
    # converted original xls file to csv in excel and renamed as "daily.csv" 
    # caribou 15 min interval records:
    qTab = read.csv(paste0(qTabPaths[i], '.csv'), header=T)
    hourlyQ = qTab$AvgOfDischarge.L.s.*0.001 # L/s to cms convert
    
    # convert hourly Q to average daily Q:
    dates = as.vector(unique(qTab$Date.and.Time))
    q = dates
    for(j in 1:length(dates)){
      q[j] = mean(hourlyQ[dates[j] == qTab$Date.and.Time], na.rm=T)
    }
    q = as.numeric(q)
    fieldQ = q[match(fieldDates[i,], dates)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
  }
  
  
  NZhourly <- function(qTabPaths, fieldDates, oTab, i){
    # NZ hourly records:
    qTab = read.csv(paste0(qTabPaths[i], '.csv'), header=T)
    hourlyQ = qTab$flow_rate_m3.s
    split = strsplit(as.vector(qTab$date.time), ' ')
    allDates = rep(NA, nrow(qTab))
    for (j in 1:length(split)){
      allDates[j] = split[[j]][1]
    }
    # average daily Q:
    dates = unique(allDates)
    q = dates
    for(j in 1:length(dates)){
      q[j] = mean(hourlyQ[dates[j] == allDates], na.rm=T)
    }
    
    fieldQ = q[match(fieldDates[i,], dates)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
  }
  
  
  NZdaily <- function(qTabPaths, fieldDates, oTab, i){
    # analyze NZ daily records:
    qTab = read.csv(paste0(qTabPaths[i], '.csv'), header=T, skip=4)
    q = as.numeric(as.vector(qTab$Discharge..m3.s.))
    fieldQ = q[match(fieldDates[i,], qTab$Date)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
    
  }
  
  
  stony <- function(qTabPaths, fieldDates, oTab, i){
    # converted original xls file to csv in excel and renamed as "daily.csv" 
    qTab = read.csv(paste0(qTabPaths[i], '.csv'), header=T)
    hourlyQ = qTab$Q..L.s..conversion.using.rating.curve.
    
    split = strsplit(as.vector(qTab$Timestamp), ' ')
    allDates = rep(NA, nrow(qTab))
    for (j in 1:length(split)){
      allDates[j] = split[[j]][1]
    }
    
    # average daily Q:
    dates = unique(allDates)
    q = rep(NA, length(dates))
    for(j in 1:length(q)){
      q[j] = mean(hourlyQ[dates[j] == allDates], na.rm=T)
    }
    q = as.numeric(q)
    fieldQ = q[match(fieldDates[i,], dates)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)*1e-3 # L/s to cms
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
    
  }
  
  
  stony_sub <- function(qTabPaths, fieldDates, oTab, i){
    # converted original xls file to csv in excel and renamed as "daily.csv" 
    qTab = read.csv(paste0(qTabPaths[i], '.csv'), header=T)
    hourlyQ = qTab$Q
    
    split = strsplit(as.vector(qTab$date), ' ')
    allDates = rep(NA, nrow(qTab))
    for (j in 1:length(split)){
      allDates[j] = split[[j]][1]
    }
    
    # average daily Q:
    dates = unique(allDates)
    q = rep(NA, length(dates))
    for(j in 1:length(q)){
      q[j] = mean(hourlyQ[dates[j] == allDates], na.rm=T)
    }
    q = as.numeric(q)
    fieldQ = q[match(fieldDates[i,], dates)]
    fieldQrange = range(as.numeric(fieldQ), na.rm=T)*1e-3 # L/s to cms
    fieldQmean = mean(fieldQrange)
    
    cdf = ecdf(q)
    cdfRange = range(100*cdf(fieldQ), na.rm=T)
    cdfMean = mean(cdfRange)
    
    oTab$QRecLength_yrs[i] =  round(length(q)/365)
    oTab$Q_cms[i] = paste0(round(fieldQmean, 3), '±', round(fieldQmean-fieldQrange[1], 3))
    oTab$Q_percentile[i] = paste0(round(cdfMean, 3), '±', round(cdfMean-cdfRange[1], 3))
    oTab$runoff_mmPerDay[i] =  paste0(round(1e-3*60*60*25*fieldQmean/oTab$gauge_A_km2[i], 3), 
                                      '±', round(1e-3*60*60*25*(fieldQmean-fieldQrange[1])/oTab$gauge_A_km2[i], 3))# basin-averaged runoff
    
    return(oTab)
    
  }
  
  
  
}
