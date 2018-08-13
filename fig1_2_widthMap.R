fig1_2_widthMap <- function(inTabPaths, tabNames, pdfOut){

  # widthMap.R
  # by George Allen, Feb 2016

  # runs through each stream centerline, calculates the orthogonal
  # direction to the along stream direction at each vertex.
  # shapefiles need with UTM coordinates and also joined
  # field data for this to work.

  require(foreign)

  ############################
  # input parameters:
  n = 5 # N vertices overwhich to calculate direction (must be odd numbers >1)
  wt = c(5,5,3,1,1)/15 # weights for the weighted mean calculation
  wMult = 20 # multiplier for the displayed cross section multiplier
  setPlotLim = T
  rotation = c(120, -140, 90, -110, -160, -150, -90,
               rep(0, 6)) # angle to rotate figure
  panelScale = c(rep(1, 7), rep(0.6, 6)) # X size of each panel plot

  # stream shapefile DBF paths:
  # TODO: make a separate file with these gis data, and provide that as an arg to the function.
  inDbfPaths = sub('field/', 'field/shapefiles/', inTabPaths)
  inDbfPaths = sub('field_dat.csv', 'segments_pts.dbf', inDbfPaths)

  # watershed  outline shapefile DBF paths:
  # TODO: make a separate file with these gis data, and provide that as an arg to the function.
  inShedDbfPaths = sub('field/', 'field/shapefiles/', inTabPaths)
  inShedDbfPaths = sub('field_dat.csv', 'shedOutline.dbf', inShedDbfPaths)


  ############################
  pdf(pdfOut, width=4, height=11, fillOddEven=T)

  layoutTab = rbind(c(1,1,1,1,1,2,2,2,2),
                    c(1,1,1,1,1,2,2,2,2),
                    c(1,1,1,1,1,2,2,2,2),
                    c(1,1,1,1,1,2,2,2,2),
                    c(3,3,3,3,3,4,4,4,4),
                    c(3,3,3,3,3,4,4,4,4),
                    c(3,3,3,3,4,4,4,4,4),
                    c(5,5,6,6,6,7,7,7,7),
                    c(5,5,6,6,6,7,7,7,7),
                    c(5,5,6,6,6,7,7,7,7),
                    c(8,8,8,9,9,9,10,10,10),
                    c(8,8,8,9,9,9,10,10,10),
                    c(11,11,11,12,12,12,13,13,13),
                    c(11,11,11,12,12,12,13,13,13))

  layout(layoutTab)
  op = par(oma = c(0,0,2,0),
           mar = c(0,0,0,0))

  ############################

  # loop through each catchment and plot stream networks and widths:
  for (m in 1:length(inDbfPaths)){

    # import watershed outline dbf:
    shed = read.dbf(inShedDbfPaths[m])
    if ('dbf' %in% names(shed)){ shed = shed$dbf }

    # apply rotation transformation around center of watershed:
    rShed = rotator(shed$POINT_X, shed$POINT_Y, rotation[m])

    # determine plot bounds to scale widths:
    plot(rShed, type='n', asp=1, axes=F, ann=F)
    pltSz = par()$fin
    xR = range(rShed[,1])
    yR = range(rShed[,2])
    # find the limiting dimension:
    pScl = max(c(xR[2]-xR[1], yR[2]-yR[1])/pltSz)
    # multiply by a scaler to make widths a reasonable thickness:
    wScl = pScl*0.00032

    ############################
    # import and process stream data:
    tab = data.frame(read.dbf(inDbfPaths[m]))

    # standardize the header names:
    names(tab) = tolower(names(tab))
    names(tab) = sub('dbf.', '', names(tab))

    # calculate cross sectional direction at each vertex:
    x = tab$lon
    y = tab$lat
    s = tab$name
    l = nrow(tab)
    w =  tab$width_cm*0.01*wScl*wMult # convert to meters

    # smooth widths with a gaussian kernel:
    splF = splinefun(w)
    w = ksmooth(c(1:l), splF(1:l), "normal", bandwidth=5)$y

    # chop start and end of vectors calculate bearing between neighbors:
    p1x = x[-c((l-n+2):l)]
    p1y = y[-c((l-n+2):l)]
    p2x = x[-c(1:(n-1))]
    p2y = y[-c(1:(n-1))]

    # calculate centerline angle:
    a = atan2(p2y-p1y, p2x-p1x)*180/pi

    # make a original length
    a = c(rep(-999, floor(n/2)), a, rep(-999, floor(n/2)))

    #### handle start and end of segments (where angles get funky):

    # locate where new segments start:
    j = which(!duplicated(s))

    # insert NAs at start, end, and jump in vector:
    for (i in rev(1:length(j))){
      x = insertRow(x, NA, j[i])
      y = insertRow(y, NA, j[i])
      a = insertRow(a, NA, j[i])
    }

    x = c(x, NA)
    y = c(y, NA)
    a = c(a, NA)

    # get bounds of NA values:
    jNA = which(is.na(a))

    closeL = jNA[-1] - 1
    closeR = jNA[-length(jNA)] + 1
    farL = jNA[-1] - floor(n/2)
    farR = jNA[-length(jNA)] + floor(n/2)

    # use a linearly shrinking window to calculate bearing at ends of vectors:
    for (i in 1:(length(jNA)-1)){

      fL = farL[i]:closeL[i]
      rL = closeR[i]:farR[i]

      for (ii in 1:length(fL)){

        # calculate all points on left sides of jumps:
        L = c((fL[ii]-floor(n/2)), closeL[i])
        a[fL[ii]] = atan2((y[L[2]]-y[L[1]]), (x[L[2]]-x[L[1]]))*180/pi

        # handle all points on right sides of vectors:
        R = c(closeR[i], (rL[ii]+floor(n/2)))
        a[rL[ii]] = atan2((y[R[2]]-y[R[1]]), (x[R[2]]-x[R[1]]))*180/pi

      }
    }

    # remove NAs from vectors:
    x = na.omit(x)
    y = na.omit(y)
    a = na.omit(a)

    ############################
    # find XY of stream edges:
    q = 90-a
    q[q < 0] = q[q < 0] + 360

    o1x = x + cos(q*pi/180)*w
    o1y = y - sin(q*pi /180)*w
    o2x = x - cos(q*pi/180)*w
    o2y = y + sin(q*pi/180)*w

    # set XY coordinates with a zero width to NA:
    zW = tab$width_cm==0

    o1x[zW] = NA
    o1y[zW] = NA
    o2x[zW] = NA
    o2y[zW] = NA

    # add an NA between segments to plot centerline:
    for (i in rev(1:length(j))){
      x = insertRow(x, NA, j[i])
      y = insertRow(y, NA, j[i])
      o1x = insertRow(o1x, NA, j[i])
      o1y = insertRow(o1y, NA, j[i])
      o2x = insertRow(o2x, NA, j[i])
      o2y = insertRow(o2y, NA, j[i])
    }

    # organize XY for plotting polygons:
    ox = mapply(c, grouper(o1x), lapply(grouper(o2x), rev), NA)
    oy = mapply(c, grouper(o1y), lapply(grouper(o2y), rev), NA)

    ox = unlist(ox, use.names=F)
    oy = unlist(oy, use.names=F)


    # apply rotation transformation around center of watershed:
    cntr = colMeans(rShed, na.rm=T)
    rO = rotator(ox, oy, rotation[m], cntr)
    #clr0 = rotator(x, y, rotation[m], cntr)
    #seg1r0 = rotator(o1x, o1y, rotation[m], cntr)
    #seg2r0 = rotator(o2x, o2y, rotation[m], cntr)
    #############################
    #PLOT:

    # Watershed outline:
    #plot(rShed, type='n', asp=1, axes=F, ann=F)
    if (m > 7){
      polygon(rShed, density=0, border=T, lwd=0.5, col="red")
    }else{
      polygon(rShed, density=0, border=T, lwd=0.5)
    }

    # add subcatchment boundaries to panel g:
    if (m == 7){
      shed = read.dbf(inShedDbfPaths[m+1])
      if ('dbf' %in% names(shed)){ shed = shed$dbf }
      # apply rotation transformation around center of watershed:
      rShed = rotator(shed$POINT_X, shed$POINT_Y, rotation[m], cntr)
      polygon(rShed, density=0, border=T, lwd=0.5, col="red")
    }

    # add stream network:
    polygon(rO, col="blue", border=NA, fillOddEven=T)

    #lines(clr0, col=2, lwd=0.1)
    #segments(seg1r0[,1], seg1r0[,2],
    #         seg2r0[,1], seg2r0[,2], col=2, lwd=0.1)

    #lines(o1x, o1y, col='orange')
    #lines(o2x, o2y, col='yellow')

    # add scale bars (200 m):
    if (m == 5){
      sBar = rotator(c(cntr[1], cntr[1]), c(cntr[2]-100, cntr[2]+100), rotation[m])
    }else{
      sBar = rotator(c(xR[1], xR[1]), c(yR[1], yR[1]+200), rotation[m])
    }

    # add scale bar:
    if (m < 9){
      arrows(sBar[1], sBar[3], sBar[2], sBar[4], length=0.03)
    }

    # add width bar:
    if (m == 1){
      w500 = 1000*0.01*wScl*wMult
      wBarYbuf = yR[2] - 30

      wBarX = c(xR[1], xR[1]+500, xR[1]+500, xR[1])
      wBarY = c(wBarYbuf, wBarYbuf+w500/2, wBarYbuf-w500/2, wBarYbuf)
      polygon(wBarX, wBarY, col="blue", border=F)

      text(xR[1], wBarYbuf, "0", col="blue", pos=1, cex=0.8)
      text(xR[1]+500, wBarYbuf, ">1 m", col="blue", pos=1, cex=0.8)
    }
  }

  dev.off()
  # cmd = paste('open', pdfOut)
  # system(cmd)

}

############################
# functions:
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,length(existingDF)+1)] = existingDF[seq(r,length(existingDF))]
  existingDF[r] = newrow
  return(existingDF)
}

grouper <- function(x){
  idx = 1 + cumsum(is.na(x))
  nonNa = !is.na(x)
  split(x[nonNa], idx[nonNa])
}

rotator <- function(x, y, theta, rot_cntr=NULL){
  M = cbind(x, y)
  if (is.null(rot_cntr)){
    rot_cntr = colMeans(M, na.rm=T)
  }
  a = theta*(pi/180)
  rotM = matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol=2)
  return(t(rotM %*% (t(M) - rot_cntr) + rot_cntr))
}
