fig3_widthModel <- function(fieldtopoPaths, tabNames, csvOut, pdfOut){

  # widthModel.R
  # George Allen, June 27, 2017
  # georgehenryallen@gmail.com

  # This script generates artifical stream width data for each basin where DEM data are
  # available. It compares the distributins between these modeled stream widths
  # to the field observed stream widths.
  # Inputs:
  # 1. topoTable: table that includes field observations (stream widths)
  #    DEM-derived data (upstream drainage area, slope)
  # 2. Data table from Trampush et al., (2014); WRR; doi:10.1002/2014WR015597;
  #    Contains international bankful width &  drainage area measurements
  # Outputs:
  # 1. Fig. 3 showing the distributions of modeled and observated stream widths.


  ############################################
  # automatically install and load libraries:
  if (!"foreign" %in% rownames(installed.packages())) {install.packages("foreign")}; require(foreign)
  if (!"MASS" %in% rownames(installed.packages())) {install.packages("MASS")}; require(MASS)
  if (!"shapefiles" %in% rownames(installed.packages())) {install.packages("shapefiles")}; require(shapefiles)

  trampushCSVpath = here('auxiliary_data', 'Trampush_etal_2014_WRR_sup_Bankful_W_Q.csv')

  fieldtopoPaths = c(
    here('locationStreamSurveys','kings.csv'),
    here('locationStreamSurveys','sagehen.csv'),
    here('locationStreamSurveys','elder.csv'),
    here('locationStreamSurveys','caribou.csv'),
    here('locationStreamSurveys','blueduck.csv'),
    here('locationStreamSurveys','stony.csv')
  )


  ##############################################

  # set Manning's:
  MannN = 0.04 #s/m^(1/3)

  # Earth's gravitational acceleration:
  g = 9.81 #m/s^2

  # convert Manning's n to Parking's k:
  k = (8.1 * g^0.5 * MannN)^6

  # analyze Trampush etal (2014) bankful width drainage area relationship:
  trampush = trampush_WidthAreaRegression(trampushCSVpath, plot=F)

  # Surveyed catchment area (from ED Table 1):
  catchment_A_m2 = c(180, 449, 362, 258, 6, 56, 128, rep(48, 6))*1e4 #m2

  # find the median bankful width to depth ratio for drainage
  # areas overlapping with surveyed catchment areas:
  trampushW2H = trampush$W2Hratio[trampush$A_m2<max(catchment_A_m2)]

  # process stream gauge data to calculate basin averaged runoff:
  # gauge upstream basin areas. Data from ED Table #1:
  gauge_A_m2 = c(10.59, 27.2, 16.8, 104, 6350, 64, 1.3, rep(0.48, 6))*1e6 #m2

  #  Data from ED Table #1:
  Q_cms = c(0.303,	0.062, 0.057,	0.137,	418.627,	0.156,	0.019,	0.008,	0.014,	0.015,	0.009,	0.013,	0.013)
  #c(0.3029814, 0.04445615, 0.08919546, 0.137, 419, 0.155, 0.0075) #m3/s
  runoff_m = Q_cms/gauge_A_m2 # m/s; basin-averaged runoff

  # set up pdf:
  pdf(pdfOut, width=7, height=7)
  layoutMatrix = rbind(c(1,1,3,3),
                       c(1,1,3,3),
                       c(2,2,4,4),
                       c(2,2,4,4),
                       c(8,8,5,5),
                       c(8,8,5,5),
                       c(8,8,6,6),
                       c(8,8,6,6),
                       c(8,8,7,7),
                       c(8,8,7,7))

  layout(layoutMatrix)
  op = par(oma=c(3,3,0,0), mar=c(1,1,3,1))


  ############################################
  # plot the first panel, showing relationship between channel shape
  # and the channel shape parameter, r:
  chanShapePlot(w = seq(0, 1, length.out=50),
                w2dratio = 8,
                r = c(1, 2, 5, 10), # shape parameter
                waterLevel = 1.09,
                colfield = F)

  ############################################
  # plot model and observed width distributions:
  # for each basin, model widths and compare to observations:
  for (i in c(1,2,3,5,6)){

    print(paste0(i, '... ', tabNames[i]))
    # read in data table:
    tab = tableReader(fieldtopoPath = fieldtopoPaths[i])

    # generate a population of r from a given distribution:
    r = runif(nrow(tab), 1, 10)
    #r = rnorm(nrow(tab), 5, 2)

    # calculate median width-depth ratio from Trampush dataset:
    WDR = median(trampushW2H, na.rm=T)

    # stream width model (equation 3):
    wMod_raw = widthModel(
                      a = exp(trampush$model$coefficients[[1]]),
                      A = tab$acc_area_m2, # m2
                      b = trampush$model$coefficients[[2]],
                      Q = runoff_m[i]*tab$acc_area_m2, # m3/s,
                      g = 9.8, # m/s2
                      S = abs(tab$slope),
                      k = k, # m; roughness length scale from Parker (1991)
                      bfWidth2DepthRatio = WDR,
                      r = r
                    )

    # remove any infinite or zero width values. These are rare but
    # trip up fitting a lognormal distribution to the data:
    wMod_raw = wMod_raw * 100 # m to cm convert
    finiteWidths = wMod_raw > 0 & wMod_raw < 1e5 & is.finite(wMod_raw)
    modeledWidths = wMod_raw[finiteWidths]

    # add model data to table:
    tab$r_param = r
    zeroWidths = which(tab$w_model==0)
    tab$w_model = round(wMod_raw, 2)
    tab$w_model[zeroWidths] = 0
    csvOutPath = paste0(csvOut, '_', tabNames[i], '.csv')
    write.csv(tab, csvOutPath, row.names=T)

    # MLE fit lognormal distribution to modeled and observed width data:
    modelFit = distribAnalyzer(widths=modeledWidths, binInterval=10)
    observedFit = distribAnalyzer(widths=tab$flowing_width, binInterval=10)

    # Plot histograms and fits of modeled and observed width data:
    distribPlotter(modelFit, observedFit, binInterval,
                   plotXmax=300, tabNames[i], i)

    options(warn=-1)

    # Pearson's Chi-squared test for count data:
    chi = chisq.test(tab$flowing_width[finiteWidths], modeledWidths) # low p-value signifies a good fit
    print(paste(i, tabNames[i],
                "    X2:", round(chi$statistic[[1]], 4),
                "    p:", round(chi$p.value[[1]], 4)))

    # Kolmogorov-Smirnov test bewteen raw model and observed data:
    ks = ks.test(tab$flowing_width, modeledWidths)
    print(paste(i, tabNames[i],
                "    ks D:", round(ks$statistic[[1]], 4),
                "    ks p:", round(ks$p.value[[1]], 4)))
    options(warn=0)

    # plot a map of the modeled widths:
    # Figure 1 - stream width map generator:
    #source(paste0(wd, '/R/widthMap.R'))
    #pdfOut = 'E:/misc/2015_09_01_Small_Stream_Width_Fieldwork/misc/mikeLambModel/mapsOfModelWidths/modelMap.pdf'
    #widthMap(tabOutPath, "konza", pdfOut)

  }

  dev.off()
  # cmd = paste('open', pdfOut)
  # system(cmd)

}



##############################################
# functions:

chanShapePlot = function(w, w2dratio, r, waterLevel, colfield=F){

  # plot the width-height relationship with different shape parameters
  wbf = max(w)
  hbf = 2*wbf/w2dratio
  N = length(w)
  lev = round(N/waterLevel)

  # for color-field image:
  if (colfield ==T){
    x = seq(min(r), max(r), length.out=1e2)
    lx = x^2
    r = (lx-min(lx))/(max(lx)-min(lx))*(max(x)-min(x))+min(x)
  }
  n = length(r)
  cols = rainbow(n)

  # plot width-height relationship with different shape params:
  # first part:
  plot(c(-1*wbf, wbf), c(0, hbf), type='n', las=1,
       xlab = "Channel width", ylab = "Channel height",
       ann=F, axes=F,
       main="Modeled channel geometries")
  for (i in 1:n){
    h = hbf*(w/wbf)^r[i]
    lines(c(rev(-1*w), w), c(rev(h), h), col=cols[i])
  }
  if(n < 5){
    legVals = as.vector(quantile(r, seq(0,1,length.out=n)))
  }else{
    legVals = as.vector(quantile(r, seq(0,1,length.out=4)))
  }
  legend("top", paste("r =", floor(legVals)), col = rainbow(5),
         lty = c(1,1), merge = T, xjust=1, cex=0.8, box.col = NA)

  # second part:
  plot(c(-1*wbf, wbf), c(0, hbf), type='n', las=1,
       ann=F, axes=F)
  h = hbf*(w/wbf)^legVals[3]
  polygon(c(rev(-1*w[-c(lev:N)]), w[-c(lev:N)]),
          c(rev(h[-c(lev:N)]), h[-c(lev:N)]),
          col=rgb(0.85,0.95,1,1), border=1)
  lines(c(rev(-1*w), w), c(rev(h), h), col=1, lwd=1.4)
}


trampush_WidthAreaRegression <- function(trampushCSVpath, plot=FALSE){

  # fits a least squares regression on Trampush et al., (2014) data.
  csv = read.csv(trampushCSVpath, header=T)

  # exclude rows with no data:
  notNA = !is.na(csv$DA..km2.) & !is.na(csv$Wbf..m.)
  A = csv$DA..km2.[notNA]*1e6 # convert km2 to m2
  w = csv$Wbf..m.[notNA]

  # take natural log of width and area data and run regression:
  lA = log(A)
  lw = log(w)
  model = lm(lw~lA)

  # plot regression and print statistical summary:
  if(plot==T){
    plot(A, w, log='xy',
         xlab="Drainge Area (m2)",
         ylab="Bankfull Width (m)",
         pch=16, cex=.5)
    par(new=T)
    plot(lA, lw, type='n',
         xlab='', ylab='', axes=F,
         main="Trampush et al. (2014)")
    abline(model, col=2, lwd=1.8)
    print("Trampush area-width regression statistical summary:")
    print(summary(model))
  }

  list = list(
    model=model,
    W2Hratio=csv$Wbf..m./csv$Hbf..m.,
    A_m2=A
  )

  return(list)

}


tableReader <- function(fieldtopoPath){

  # read in and process input table.
  tab = read.csv(fieldtopoPath, header=T, skip = 49)
  tab$percent_nonflow[is.na(tab$percent_nonflow)] = 0
  w = tab$flowing_width*(1-tab$percent_nonflow/100)
  w[grep('RF', tab$code)] = w[grep('RF', tab$code)]*39.3701 # Rangefinder Convert
  tab$flowing_width = w * 2.54 # inches to cm convert
  notZero = w != 0
  tab = tab[notZero, ]

  return(tab)

}


widthModel <- function(a, A, b, Q, g, S, k, bfWidth2DepthRatio, r){

  # Equation 3 that models stream width:
  w = Q ^ (3/(5*r+3)) *
    (a*A^b) ^ ((r-1)/(r+3/5)) *
    (8.1 *
       (g * S)^(1/2) *
       k^(-1/6) *
       (1/bfWidth2DepthRatio)^(5/3) *
       (1-1/(r+1))) ^ (-3/(5*r+3))

  return(w)

}


binInterval = 10 # cm; histogram bin size

distribAnalyzer <- function(widths, binInterval){

  xMax = max(widths)+binInterval
  breaks = seq(0, xMax, binInterval)
  h = hist(widths, breaks, plot=F)
  fit = fitdistr(widths, "log-normal")$estimate # MLE lognormal distr
  lineSeq = seq(0, xMax, length=5e2) # sequence for line
  dln = dlnorm(lineSeq, fit[1], fit[2]) # densities
  fln = dln*length(widths)*binInterval # frequencies

  list = list(h=h, fln=fln, fit=fit, length=length(widths))

  return(list)

}





distribPlotter <- function(modelFit, observedFit, binInterval, plotXmax=NA, catchName, i){

  # calculate plot limits:
  if (is.na(plotXmax)){
    xlim = c(0, max(c(modelFit$h$breaks, observedFit$h$breaks)))
  }else{
    xlim = c(0, plotXmax)
  }
  ylim = c(0, max(c(modelFit$h$counts, modelFit$fln,
                    observedFit$h$counts, observedFit$fln)))

  # plot histograms:
  plot(observedFit$h,
       xlim=xlim, ylim=ylim,
       xlab="Stream width (cm)", ylab="N", las=1,
       col='gray', border=0, main='')
  title(paste0(letters[i], '. ', catchName), adj=0)
  par(new=T)
  plot(modelFit$h,
       xlim=xlim, ylim=ylim,
       axes=F, main="", xlab="", ylab="",
       border=2, col=2, density=20)

  # plot lognormal fit lines:
  lineSeq = seq(0, xlim[2], length=5e2)
  dln = dlnorm(lineSeq, observedFit$fit[1], observedFit$fit[2]) # densities
  fln = dln*observedFit$length*binInterval # frequencies
  lines(lineSeq, fln, lty=1, col=1, lwd=1.5)

  dln = dlnorm(lineSeq, modelFit$fit[1], modelFit$fit[2]) # densities
  fln = dln*modelFit$length*binInterval # frequencies
  lines(lineSeq, fln, col=2, lwd=1.8)

  if (i == 1){
    legend("topright", c("model", "observed"), col = c(2, 1),
           lty = c(1,1), merge = T, xjust=1, cex=0.8)
  }

}
