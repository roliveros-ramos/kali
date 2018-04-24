
.image.mapnl = function(lon, lat, z, center=0, hires=FALSE, add = FALSE, nlevel=1000, 
                        col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)),
                        land.col="darkolivegreen4", sea.col="aliceblue", boundaries.col = "black", 
                        grid.col="white", grid=FALSE, axes=TRUE, border=!axes, ...) {
  
  # lonData = .checkLongitude(lon)
  # if(!is.null(lonData$ind)) {
  #   z = z[lonData$ind, ]
  #   lon = lonData$lon
  # }
  
  image(x=lon, y=lat, z=z, col=col, axes=FALSE, add=add, xlab="", ylab="", ...)
  
  mapDetails(center=center, hires=hires,col=land.col, interior=FALSE, 
             axes=axes, border=border, boundaries.col=boundaries.col,
             grid=grid, grid.col=grid.col)  
  
  return(invisible())
}

.checkLongitude = function(lon, ...) {
  
  if(is.matrix(lon)) return(list(lon=lon))
  if(!any(lon>180)) return(list(lon=lon))
  lon[lon>180] = lon[lon>180] - 360
  xlon = sort(lon, index.return = TRUE)
  
  return(list(lon=xlon$x, ind=xlon$ix))
}

.plotSea = function(col="lightblue", x0=-300, y0=-150) {
  
  x1 = -x0
  y1 = -y0
  polygon(x=c(x0,x1,x1,x0), y=c(y0,y0,y1,y1), border=NA, col=col)
  
  return(invisible())
}

.getDomain = function(domain, axis) { 
  # create a data file
  
  if(!is.character(domain)) domain = deparse(substitute(domain))
  output = domains[[domain]][[axis]]
  
  return(output)
}

.rotate = function(x, revCol=TRUE) {
  
  # clockwise
  if(length(dim(x))==2) {
    x = t(x)
    if(revCol) x = x[, ncol(x):1]
  }
  if(length(dim(x))==3) {
    x = aperm(x, c(2,1,3))
    if(revCol) x = x[, ncol(x):1, ]
  }
  
  return(x)
}

.rotate2 = function(x, revCol=TRUE) {
  
  # anticlockwise
  if(length(dim(x))==2) {
    if(revCol) x = x[, ncol(x):1]
    x = t(x)
  }
  if(length(dim(x))==3) {
    if(revCol) x = x[, ncol(x):1, ]    
    x = aperm(x, c(2,1,3))
  }  
  
  return(x)
}

.coord2text = function(coord, type) {
  
  # write nicely coordinates
  degree = "\U00B0"
  hemi = if(coord==0) {
    rep(degree, 2)
  } else {
    if(coord>0) paste0(degree, c("N","E")) else paste0(degree, c("S","W"))
  }
  
  out = switch(type,
               lat = paste0(abs(coord), hemi[1]),
               lon = paste0(abs(coord), hemi[2]),
               as.character(coord)
  )
  
  return(out)
}


.axis.map = function(side, type, ...) {
  
  x = axTicks(side=side)
  axis(side=side, at=x, labels=coord2text(coord=x, type=type), ...)
  
  return(invisible(x))
}

# nice shape for frame plots
getmfrow = function(n) {  
  m1 = floor(sqrt(n))
  m2 = ceiling(n/m1)
  out = rev(sort(c(m1, m2)))
  return(out)
}

# rotate a image clock or anti-clockwise
rotate = function(x, direction="clockwise", revCol=TRUE) {
  
  out = switch(direction,
               clockwise     = .rotate(x, revCol=revCol),
               anticlockwise = .rotate2(x, revCol=revCol))
  
  return(out)
}

coord2text = function(coord, type) {
  
  if(!is.character(type)) type=deparse(substitute(type))
  out = sapply(coord, FUN=.coord2text, type=type)
  
  return(out)
}


