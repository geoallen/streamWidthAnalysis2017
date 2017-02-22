fig2_distributions <- function(inTabPaths, tabNames, pdfOut){
  
  ############################
  # fig2_distributions.R
  # George Allen, Sept 2016
  
  # description: plots figure 2 in allen et al., Nature
  
  # check for existence of all required files:
  if (F %in% file.exists(inTabPaths)){
    message("Field Data CSV files are missing")
  }
  
  # load requied library:
  require(MASS)
  
  ############################
  # user defined parameters:
  
  # colors pallets:
  lcols = c("#0000DB", "#009E88", "#C67A00", "#B03060", "#A0A000", 
            "#007DD1", "#C80000", rep("#000000", 6))
  fcols = c(rep("#C1C1C0", 7), "#737474", "#878787", "#9B9B9B", 
            "#B1B0B1", "#C4C4C3", "#DAD9D9")
  
  arrowLength = 0.01
  
  # histogram binning interval:
  int = 10
  # x resolution of line plots:
  dlnInt = 1
  # number of tick marks on the plot axis:
  nXtx = 10
  nYtx = 5
  
  # density kernel type and bandwidth:
  kernel="gaussian"
  bw = 10 # nrd0 # note: variance = (bw of gaussian kernel)^2
  
  ############################
  # create a table to store distribution fit parameters:
  names = c("catchment", "mu", "sigma", "maxW", "maxDensity")
  fitParams = as.data.frame(array(NA, c(length(tabNames), length(names))))
  names(fitParams) = names
  fitParams$catchment = tabNames
  
  ############################
  # open a PDF device to plot distributions:
  pdf(pdfOut, height=4.5, width=8)
  
  layoutMatrix = rbind(c(1,2,3,4,5), 
                       c(1,2,3,4,5), 
                       c(6,7,8,9,10), 
                       c(6,7,8,9,10), 
                       c(14,14,11,12,13), 
                       c(14,14,11,12,13), 
                       c(14,14,15,15,15))
  
  layout(layoutMatrix)
  op = par(oma=c(3,3,0,0), mar=c(1,1,3,1))
  
  ############################
  # Plot panels a-h:
  for (i in 1:length(inTabPaths)){
    
    print(i)
    
    # import and process data:
    w = tabReader(inTabPaths[i])
    
    # calculate mode width with a gaussian density kernel:
    kern = density(w, bw=bw, kernel=kernel)
    modeW = kern$x[which.max(kern$y)]
    
    # calculates the limits of plot and the fitted 
    # distribution curves:
    if (i < 8){
      limit = limitCalc(w, dlnInt)
    }else{
      limit = limitCalc(w, dlnInt, 350, 70)
    }
    
    # convert limit list items to objects: 
    for(j in 1:length(limit)){
      assign(names(limit)[j], limit[[j]])
    }
    
    ############################
    # Plot histograms, lognormal fits, and mode locations:
    
    h = hist(w, breaks, freq=T, plot=T, axes=F, main='', 
             xlab="W (cm)", ylab="N",
             ylim=c(0, maxY+maxY*0.1), xlim=c(0,maxX), 
             col=fcols[i], border=NA)
    
    tx = seq(0, floor(maxX/50)*50, length=nXtx)
    axis(1, tx, c(min(tx), rep(NA, nXtx-2), max(tx)), lwd=0.7)
    tx = seq(0, maxY, length=nYtx)
    axis(2, tx, c(min(tx), rep(NA, nYtx-2), max(tx)), lwd=0.7, las=1)
    
    lines(lineSeq, fln, col=lcols[i], lwd=1.4)
    
    # plot mode widths:
    arrows(modeW, max(fln)+max(fln)*arrowLength, modeW, max(fln), 
           length=.03, lwd=1, col=lcols[i])
    text(modeW+max(lineSeq)*.05, 
         max(fln, na.rm=T)+max(fln, na.rm=T)*arrowLength*2,
         paste0(round(modeW), " cm"), col=lcols[i], cex=0.96, adj=0)
  
    # add title:
    if (i <= 8){title(paste0(letters[i], '. ', tabNames[i]),
                     adj=0, line=0, cex=0.8, font=2, col.main=lcols[i])
    }
    
    if (i >= 8){
      mtext(tabNames[i], adj=0, font=3, cex=0.6, col.main=lcols[i])
    }
    
    # keep track of lognormal fit parameters:
    fitParams[i, c(2:5)] = c(lnfit[1], lnfit[2], max(w), max(dln))
  }
  
  ############################
  # plot panel i, composite density distribution functions:
  
  xlim = c(0, 250)
  ylim = c(0, max(fitParams$maxDensity))
  lineSeq = seq(0, ceiling(max(xlim)/dlnInt)*dlnInt, dlnInt)
  
  plot(xlim, ylim, type="n", axes=F, xlab="W (cm)", ylab="Density")
  tx = seq(0, xlim[2], length=nXtx)
  ylab = c(min(tx), rep(NA, nXtx-2), max(tx))
  axis(1, tx, ylab, lwd=0.7)
  
  tx = seq(0, ylim[2], length=nYtx)
  xlab = round(c(min(tx), rep(NA, nYtx-2), max(tx)), 2)
  axis(2, tx, xlab, lwd=0.7, las=1)
  
  ccols = c(lcols[1:7], fcols[8:13])
  for (i in length(inTabPaths):1){
    
    dln = dlnorm(lineSeq, fitParams$mu[i], fitParams$sigma[i])
    lines(lineSeq, dln, col=ccols[i], lwd=1.4)
    
  }
  
  # plot average distribution fit on top:
  dln = dlnorm(lineSeq, mean(fitParams$mu), mean(fitParams$sigma))
  lines(lineSeq, dln, col=1, lwd=2)
  
  text(mean(xlim), mean(ylim)+mean(ylim)/5, 
       paste("k =", round(mean(fitParams$mu),2)), pos=2)
  text(mean(xlim), mean(ylim), 
       paste("sigma =", round(mean(fitParams$sigma),3)), pos=2)
  
  title("i. All surveys", adj=0, line=1, cex=0.8, font=2)
  
  dev.off()
  cmd = paste('open', pdfOut); system(cmd)

  print(fitParams)
  #write.csv(fitParams, 'E:/misc/2015_09_01_Small_Stream_Width_Fieldwork/tables/lognormal_fitParameters.csv', row.names=F)
  
  ############################
  # functions:
  
  # reads in field data table file, extracts and processes data.
  # returns stream width measurements in cm:
  tabReader = function(filePath){
    tab = data.frame(read.csv(filePath, skip=49))
    tab$percent_nonflow[is.na(tab$percent_nonflow)] = 0
    w_raw = tab$flowing_width*(1-tab$percent_nonflow/100)
    w_raw[grep('RF', tab$code)] = w_raw[grep('RF', tab$code)]*39.3701 # Rangefinder Convert
    notZero = w_raw!=0 
    w = w_raw[notZero]
    w = w * 2.54 # inches to cm convert
    return(w)
  }
  
  # takes in width data, the X spacing interval of the lines
  # that will be plotted, and optionally the maximum X and Y 
  # values of plot. Returns the maximum X and Y values (if they
  # are not already definied), the X vector used to plot lines, 
  # and density and frequency vectors along this line:
  limitCalc = function(w, dlnInt, maxX=NA, maxY=NA){
    
    # find X limit:
    if(is.na(maxX)){
      maxX = ceiling(max(w)/dlnInt)*dlnInt
    }
    
    # calculate a sequence used to plot lines:
    lineSeq = seq(0, ceiling(maxX/dlnInt)*dlnInt, dlnInt)
    
    # caluclate liklihood at each point along the line:
    lnfit = fitdistr(w, "log-normal")$estimate
    dln = dlnorm(lineSeq, lnfit[1], lnfit[2])
    fln = dln*length(w)*int
    
    # bin widths:
    breaks = seq(0, floor(max(w)+int), int)
    h = hist(w, breaks, plot=F)
    
    # find y limit:
    if(is.na(maxY)){
      maxY = ceiling(max(c(h$counts, fln))/10)*10
    }
    
    return(list(maxX=maxX, maxY=maxY, lineSeq=lineSeq, 
                breaks=breaks, lnfit=lnfit, dln=dln, fln=fln))
    
  }
}