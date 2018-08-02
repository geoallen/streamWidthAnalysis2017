EDfig2_EDtab3_GOF <- function(inTabPaths, modTabDir, tabNames, pdfOut, csvOut) {

  options(warn=-1)

  # ED_fig2_tab3_GOF.R
  # George Allen, Feb 2016
  print(paste("Generating ED Fig. 2 - GOF of sample distributions"))

  #################################################################
  # required packages:
  require(MASS)


  #################################################################
  # custom graphing parameters:

  # histogram binning interval (cm):
  int = 10

  # x resolution of line plots:
  dlnInt = 0.1


  #################################################################
  # set up table to store all the fitted distribution and GOF statistics info:
  names = c("N", "N_BIN",
            "gshape", "rate",  "gX2", "gX2p", "gKS_D", "gKS_p",
            "meanln", "sdlog", "lnX2", "lnX2p", "lnKS_D", "lnKS_p",
            "wshape", "scale", "wX2", "wX2p", "wKS_D", "wKS_p",
            "pXmin", "pAlpha", "pX2", "pX2p", "pKS_D", "pKS_p",
            "modMean", "modSD", "modX2", "modX2p", "modKS_D", "modKS_p")
  statTab = data.frame(array(NA, c(13, length(names))))
  names(statTab) = names


  #################################################################
  # open PDF device to write out ED Fig. 2:

  pdf(pdfOut, height=9.8, width=7.3)
  layout(rbind(c(1,3,5), c(2,4,6),
               c(7,9,11), c(8,10,12),
               c(13,15,17), c(14,16,18),
               c(19,21,23), c(20,22,24),
               c(25,27,29), c(26,28,30)))


  op = par(oma = c(2,5,2,2), mar = c(2,2,2,2))

  #pdf(pdfOut, height=8, width=7.2)
  #op = par(oma = c(2,4,0,0), mar = c(2,1,0,0))


  #################################################################
  # for each stream survey:
  for (i in 1:length(inTabPaths)){

    print(paste(i, "running", tabNames[i], "..."))

    # import and process data:
    tab = data.frame(read.csv(inTabPaths[i], skip=49))
    tab$percent_nonflow[is.na(tab$percent_nonflow)] = 0
    w_raw = tab$flowing_width*(1-tab$percent_nonflow/100)
    w_raw[grep('RF', tab$code)] = w_raw[grep('RF', tab$code)]*39.3701 # Rangefinder Convert
    notZero = w_raw!=0 & !is.na(w_raw)
    w = w_raw[notZero]
    w = w * 2.54 # inches to cm convert

    # calculate mode first order stream width for Pareto fit:
    # find first order mean stream width:
    fOw_raw = w_raw[tab$stream_order == 1]
    notZero = (fOw_raw!=0 & !is.na(fOw_raw))
    fOw = fOw_raw[notZero]*2.54# inches to cm convert
    medFoW = median(fOw, na.rm=T)
    modeFoW = density(fOw, bw=10, kernel="gaussian", na.rm=T)
    modeFoW = modeFoW$x[which.max(modeFoW$y)]
    minW = modeFoW

    # set up x sequence to plot lines:
    lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)

    # get list of model table paths:
    modTabPaths = list.files(modTabDir, '.csv', full.names=T)

    #################################################################
    # Distribution fitting with maximum liklihood estimation:

    # gamma
    gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
    dg = dgamma(lineSeq, gfit[1], gfit[2])
    pg = pgamma(lineSeq, gfit[1], gfit[2])
    # logNormal
    lnfit = fitdistr(w, "log-normal")$estimate
    dln = dlnorm(lineSeq, lnfit[1], lnfit[2])
    pln = plnorm(lineSeq, lnfit[1], lnfit[2])
    # weibull
    wfit = suppressWarnings(fitdistr(w, "weibull")$estimate)
    dw = dweibull(lineSeq, wfit[1], wfit[2])
    pw = pweibull(lineSeq, wfit[1], wfit[2])
    # pareto
    parfit = pareto.mle(w[w>minW], minW)#modeW)
    dpar = dpareto(lineSeq, parfit[[1]], parfit[[2]])
    ppar = ppareto(lineSeq, parfit[[1]], parfit[[2]])

    # didn't try ftting beta, f, log-logistic, or chi squared distributions (but they prob don't fit)
    # poisson did not fit at all.

    #################################################################
    # Quantify goodness of fit:

    # Pearson's Chi-squared test for count data:
    breaks = seq(0, (max(w)+int), int)
    h = hist(w, breaks, plot=F)
    # gamma:
    f = suppressWarnings(dgamma(h$mids, gfit[1], gfit[2]))
    gChi = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=1e3)
    # lognormal:
    f = dlnorm(h$mids, lnfit[1], lnfit[2])
    lnChi = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=1e3) # low p-value signifies a good fit
    # weibull:
    f = suppressWarnings(dweibull(h$mids, wfit[1], wfit[2]))
    wChi = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=1e3)
    # pareto
    f = dpareto(h$mids[h$mids>minW], parfit[[1]], parfit[[2]])
    parChi = chisq.test(h$counts[h$mids>minW], p=f, rescale.p=T, simulate.p.value=T, B=1e3)

    # Two sided One sample KS GOF test:
    jw = jitter(w) # to remove ties
    # gamma
    gamks = ks.test(jw, "pgamma", gfit[1], gfit[2], alternative="two.sided")
    # logNormal
    lnks = ks.test(jw, "plnorm", lnfit[1], lnfit[2], alternative="two.sided")
    # weibull
    weibks = ks.test(jw, "pweibull", wfit[1], wfit[2], alternative="two.sided")
    # pareto
    parks = pareto.test(jw[jw>minW], minW, 2e2) #modeW


    #################################################################
    # Quantify differences between surveyed widths to modeled widths (Fig. 3):
    # if a file contained model data exists, read in:
    modTabPath = modTabPaths[grep(tabNames[i], modTabPaths, ignore.case=T)]

    if (length(modTabPath) > 0){
      modTab = read.csv(modTabPath, header=T)
      # remove any infinite or zero width values. These are rare but
      # trip up fitting a lognormal distribution to the data:
      finiteWidths = modTab$w_model > 0 & modTab$w_model < 1e5 & is.finite(modTab$w_model)
      # MLE fit lognormal distribution to modeled width data:
      modFit = distribAnalyzer(widths=modTab$w_model[finiteWidths], binInterval=10)$fit
      # Pearson's Chi-squared test for count data:
      modChi = chisq.test(modTab$flowing_width[finiteWidths], modTab$w_model[finiteWidths]) # low p-value signifies a good fit
      # Kolmogorov-Smirnov test bewteen raw model and observed data:
      modks = ks.test(modTab$flowing_width,  modTab$w_model)
    } else {
      modFit[[1]] = modFit[[2]] =
        modChi$statistic[[1]] = modChi$p.value[[1]] =
        modks$statistic[[1]] = modks$p.value[[1]][[1]] = NA
    }


    #################################################################
    # add distribution fit and GOF statastics to a table:
    statTab[i, ] = c(length(w), length(h$counts),
                     gfit[[1]], gfit[[2]], gChi$statistic[[1]], gChi$p.value[[1]], gamks$statistic[[1]], gamks$p.value[[1]],
                     lnfit[[1]], lnfit[[2]], lnChi$statistic[[1]], lnChi$p.value[[1]], lnks$statistic[[1]], lnks$p.value[[1]][[1]],
                     wfit[[1]],  wfit[[2]], wChi$statistic[[1]], wChi$p.value[[1]], weibks$statistic[[1]], weibks$p.value[[1]],
                     parfit[[1]], parfit[[2]], parChi$statistic[[1]], parChi$p.value[[1]], parks$D[[1]], parks$p[[1]],
                     modFit[[1]], modFit[[2]], modChi$statistic[[1]], modChi$p.value[[1]], modks$statistic[[1]], modks$p.value[[1]][[1]])

    #################################################################
    # Plot figure:

    if (i==8){par(mfg=c(7,1))}

    # PDF PLOTS (upper panel):

    # add histogram:
    maxX = max(lineSeq)
    maxY =  max(c(h$density, dln, dg))
    h = hist(w, breaks, freq=F, plot=T,
             ylim=c(0, maxY), xlab="",
             main='', xaxt='n',
             col="light gray", border=0, ylab="D(w)", las=1)

    # add title:
    if (i < 8){title(paste0(letters[i], '. ', tabNames[i]), adj=0, line=1, font=2)
    }else{
      title(tabNames[i], adj=0, line=1, font=2)
    }

    # add chi square GOF statistics:
#    text(maxX, maxY-2*maxY/6, paste0("X2=", round(gChi$statistic),
#                                      ",  p=", round(gChi$p.value, 3)), pos=2, cex=0.75, col=4)
#    text(maxX, maxY-3*maxY/6, paste0("X2=", round(lnChi$statistic),
#                                      ",  p=", round(lnChi$p.value, 3)), pos=2, cex=0.75, col=2)
#    text(maxX, maxY-4*maxY/6, paste0("X2=", round(wChi$statistic),
#                                      ",  p=", round(wChi$p.value, 3)), pos=2, cex=0.75, col="orange")
#    text(maxX, maxY-5*maxY/6, paste0("X2=", round(parChi$statistic),
#                                      ",  p=", round(parChi$p.value, 3)), pos=2, cex=0.75, col=rgb(0.4,0.4,0.4))

    # add fitted distributions:
    lines(lineSeq[dpar<maxY&dpar!=0], dpar[dpar<maxY&dpar!=0], col=rgb(.4,.4,.4))
    lines(lineSeq, dw, col="orange")
    lines(lineSeq, dln, col=2)
    lines(lineSeq, dg, col=4)

    # CDF PLOTS (lower panel):

    # generate emirical cdf:
    eCDF = ecdf(w)

    # generate fitted distribution CDF data:
    rgam = rgamma(1e4, gfit[1], gfit[2])
    rln = rlnorm(1e4, lnfit[1], lnfit[2])
    rweib = rweibull(1e4, wfit[1], wfit[2])
    rpar = rpareto(1e4, parfit[[1]], parfit[[2]])

    gCDF = ecdf(rgam)
    lnCDF = ecdf(rln)
    wCDF = ecdf(rweib)
    parCDF = ecdf(rpar)

    # add empirical CDF:
    plot(lineSeq, eCDF(lineSeq),
         xlim=c(0, max(w)+1), ylim=c(0, 1), axes=F, main='',
         xlab="w (cm)", ylab="P(w)", type='l', lwd=2.5)
    axis(1)
    axis(2, las=1)
    options(scipen = -5)

    # add KS GOF test statistics:
#    text(max(lineSeq), 0.5, paste0("D=", round(gamks$statistic,3),
#                                   ",  p=", round(gamks$p.value,3)), cex=0.75, pos=2, col=4)
#    text(max(lineSeq), 0.35, paste0("D=", round(lnks$statistic,3),
#                                   ",  p=", round(lnks$p.value,3)), cex=0.75, pos=2, col=2)
#    text(max(lineSeq), 0.2, paste0("D=", round(weibks$statistic,3),
#                                   ",  p=", round(weibks$p.value,3)), cex=0.75, pos=2, col="orange")
#    text(max(lineSeq), 0.05, paste0("D=", round(parks$D,3),
#                                   ",  p=", round(parks$p,3)), cex=0.75, pos=2, col=rgb(0.4,0.4,0.4))

    # add fitted distributions:
    options(scipen = 10)
    lines(lineSeq, parCDF(lineSeq), col=rgb(.4,.4,.4), lty=1, lwd=1)
    lines(lineSeq, wCDF(lineSeq), col="orange", lty=1, lwd=1)
    lines(lineSeq, lnCDF(lineSeq), col=2, lty=1, lwd=1)
    lines(lineSeq, gCDF(lineSeq), col=4, lty=1, lwd=1)

  }

  options(warn=0)

  # display PDF:
  dev.off()
  cmd = paste('open', pdfOut)
  system(cmd)

  # return stat Tab:
  oTab = as.data.frame(t(statTab))
  oTab = round(oTab, 3)
  colnames(oTab) = tabNames

  write.csv(oTab, csvOut)
  cmd = paste('open', csvOut)
  system(cmd)

}

#################################################################
# specialized functions for the pareto fit:
# code based on Elvis's reply on stackexhange forum:
# http://stats.stackexchange.com/questions/78168/how-to-know-if-my-data-fits-pareto-distribution

# distribution, cdf, quantile and random functions for Pareto distributions:
dpareto <- function(x, xm, a) ifelse(x > 0, a*xm**a/(x**(a+1)), 0)
ppareto <- function(q, xm, a) ifelse(q > 0, 1 - (xm/q)**a, 0 )
qpareto <- function(p, xm, a) ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/a))
rpareto <- function(n, xm, a) qpareto(runif(n), xm, a)

# fit pareto distribution to data with MLE:
pareto.mle <- function(x, xm=min(x)){
  a = length(x)/(sum(log(x))-length(x)*log(xm))
  return( list(xm = xm, a = a))
}

# compute the KS statistic, and uses parametric bootstrap to estimate the p-value.
pareto.test <- function(x, xm=min(x), B=2.5e2){
  a = pareto.mle(x, xm)

  # KS statistic
  D = ks.test(x, function(q) ppareto(q, a$xm, a$a))$statistic

  # estimating p value with parametric bootstrap
  n = length(x)
  emp.D = numeric(B)
  for(b in 1:B){
    xx = rpareto(n, a$xm, a$a);
    aa = pareto.mle(xx)
    emp.D[b] = ks.test(xx, function(q) ppareto(q, aa$xm, aa$a))$statistic
  }

  return(list(xm = a$xm, a = a$a, D = D, p = sum(emp.D > D)/B))
}

# used to analyze distribution of modeled widths:
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
