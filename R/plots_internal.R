
.image.mapnl = function(lon, lat, z, zlim, center=0, hires=FALSE, add = FALSE, nlevel=1000, 
                        col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)), land=TRUE,
                        land.col="darkolivegreen4", sea.col="aliceblue", boundaries.col = "black", 
                        grid.col="white", grid.lwd=0.5, grid=FALSE, axes=TRUE, border=!axes, labels=TRUE, ...) {

  pm = .findPrimeMeridian(as.numeric(lon))
  
  if(is.matrix(lon) & is.matrix(lat)) {
    poly.image(x=lon, y=lat, z=z, zlim=zlim, col=col, axes=FALSE, add=add, xlab="", ylab="", ...)
    map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE, 
                axes=axes, border=border, boundaries.col=boundaries.col, nx=nrow(z), ny=ncol(z),
                grid=FALSE, grid.col=grid.col, grid.lwd=grid.lwd, water=sea.col, land=land, labels=labels)  
  } else {
    image(x=lon, y=lat, z=z, zlim=zlim, col=col, axes=FALSE, add=add, xlab="", ylab="", ...)
    map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE, 
                axes=axes, border=border, boundaries.col=boundaries.col, nx=nrow(z), ny=ncol(z),
                grid=grid, grid.col=grid.col, grid.lwd=grid.lwd, water=sea.col, land=land, labels=labels)  
    
  }


  
  
  return(invisible())
}

.checkLongitude = function(lon, ...) {
  
  if(is.matrix(lon)) return(list(lon=lon))
  if(!any(lon>180)) return(list(lon=lon))
  lon[lon>180] = lon[lon>180] - 360
  xlon = sort(lon, index.return = TRUE)
  
  return(list(lon=xlon$x, ind=xlon$ix))
}

.plotSea = function(col="lightblue", x0=-400, y0=-150) {
  
  x1 = -x0
  y1 = -y0
  polygon(x=c(x0,x1,x1,x0), y=c(y0,y0,y1,y1), border=NA, col=col)
  
  return(invisible())
}

.getDomain = function(domain, axis) { 
  # create a data file
  if(is.null(domain)) return(NULL)
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


.axis.map = function(side, type, usr=NULL, n=5, ...) {
 
  if(is.na(side)) return(invisible())
  
  old_usr = par("usr")
  on.exit(par(usr=old_usr))
  if(is.null(usr)) usr = old_usr
  par(usr=usr)
  
  is.x <- side%%2 == 1
  
  x = pretty(usr[if(is.x) 1:2 else 3:4], n=n)
  xc = checkLongitude(x, "center")
  
  axis(side=side, at=x, labels=coord2text(coord=xc, type=type), ...)
  
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


findXlim = function(x) {
  xc = diff(range(checkLongitude(x, "center"), na.rm=TRUE))
  xl = diff(range(checkLongitude(x, "left"), na.rm=TRUE))
  pm = if(xc < xl) "center" else "left"
  out = range(pretty(checkLongitude(x, pm)), na.rm=TRUE)
  if(pm=="center") out = pmax(-180, pmin(out, 180))
  if(pm=="left") out = pmax(0, pmin(out, 360))
  attr(out, "pm") = pm
  return(out)
}

addPM = function(xlim) {
  pm = attr(xlim, "pm")
  if(!is.null(pm)) return(xlim)
  dd = diff(xlim)
  v0 = if(xlim[2]>=0) "P" else "N"
  v1 = if(xlim[1]>=0) "P" else "N"
  v2 = if(dd>=0) "C" else "D"
  if(v1=="P" & v2=="C") attr(xlim, "pm") = .findPrimeMeridian(xlim)
  if(v1=="P" & v2=="D") {
    if(xlim[1]>=180) {
      xlim = sort(checkLongitude(xlim, "center"))
      attr(xlim, "pm") = "center"
    } else {
      if(v0=="N") {
        xlim = checkLongitude(xlim, "left")
        attr(xlim, "pm") = "left"
      } else {
        xlim = c(-180, 180)
        attr(xlim, "pm") = "center"
      }
    }
  }
  if(v1=="N" & v2=="C") attr(xlim, "pm") = "center"
  if(v1=="N" & v2=="D") {
    xlim = c(-180, 180)
    attr(xlim, "pm") = "center"
  }
  return(xlim)
}
