# create color.themes environment
# addColorTheme() # names for the default color themes
# addDomain()
# setDefaultColorTheme()
# setDefaultDomain()

createGridAxes = function(lat,lon, dx, dy=dx) {
  # Create a rectangular grid given lat, lon and dxy.
  # No correction by Earth curvature
  if(dx <= 0 || dy <= 0) stop("dx and dy must be positive.")
  
  lats.rho = seq(from=min(lat),to=max(lat),by=dy)
  lons.rho = seq(from=min(lon),to=max(lon),by=dx)
  
  lat[which.min(lat)] = lat[which.min(lat)] - 0.5*dy
  lat[which.max(lat)] = lat[which.max(lat)] + 0.5*dy
  lon[which.min(lon)] = lon[which.min(lon)] - 0.5*dx
  lon[which.max(lon)] = lon[which.max(lon)] + 0.5*dx
  
  lats.psi = seq(from=min(lat),to=max(lat),by=dy)
  lons.psi = seq(from=min(lon),to=max(lon),by=dx)
  
  rho = list(lat=lats.rho, lon=lons.rho)
  psi = list(lat=lats.psi, lon=lons.psi)

  nlat = length(rho$lat)
  nlon = length(rho$lon)
  
  LAT = matrix(rho$lat, ncol=nlat, nrow=nlon, byrow=TRUE)
  LON = matrix(rho$lon, ncol=nlat, nrow=nlon)
  
  output = list(lon=lons.rho, lat=lats.rho, rho=rho, psi=psi, LON=LON, LAT=LAT)
  
  return(output)
}

createTimeAxis = function(start, end, frequency=12, center=FALSE, shift=TRUE) {
  times  = seq(from=start[1] + (start[2]-1)/frequency, 
               to= end[1] + end[2]/frequency, by=1/frequency)
  if(center) {
    out = list(bounds = times, center = times[-length(times)] + 0.5/frequency)
    out$year = floor(out$center)
    out$month = round(frequency*(out$center-out$year) + 0.5, 0)
    out$season = out$month
    seasons = matrix(if(shift) 1:12 else c(12, 1:11), ncol=3, byrow=TRUE)
    out$season[out$season %in% seasons[1,]] = "summer"
    out$season[out$season %in% seasons[2,]] = "fall"
    out$season[out$season %in% seasons[3,]] = "winter"
    out$season[out$season %in% seasons[4,]] = "spring"
    times = out
  } 
  return(times)
}

# createCoordinatesAxis = function(lat, lon, dx, dy) {
#   # create rho coordinate axis
#   lats = seq(from=min(lat),to=max(lat),len=dy+1)
#   lats = 0.5*(lats[-length(lats)] + lats[-1])
#   
#   lons = seq(from=min(lon),to=max(lon),len=dx+1)
#   lons = 0.5*(lons[-length(lons)] + lons[-1])
#   
#   lat_rho = matrix(lats, nrow=dx, ncol=dy, byrow=TRUE)
#   lon_rho = matrix(lons, nrow=dx, ncol=dy)
#   
#   dx = mean(diff(lons))
#   dy = mean(diff(lats))
#   
#   return(list(lat=lats, lon=lons, dx=dx, dy=dy, lat_rho=lat_rho, lon_rho=lon_rho))
# }


mapDetails = function(center, hires=FALSE, col="black", interior=FALSE, 
                      axes=TRUE, border=TRUE, boundaries.col="black",
                      grid=TRUE, grid.col="white", ...) {
  mapa = if (hires) {
    require(mapdata)
    "worldHires"
  }
  else {
    require(maps)
    "world"
  }
  
  if(isTRUE(grid)) grid(col=grid.col, lty=1)
  
  map2(mapa, center = center, fill = TRUE, col = col, add = TRUE, 
       interior=interior, border=boundaries.col, ...)
  if(axes) {
    map.axes2()
    mtext("LONGITUDE", 1, line = 1.8, cex = 0.9*par("cex"))
    mtext("LATITUDE", 2, line = 2.4, cex = 0.9*par("cex"))    
  } else {
    if(border) box()
  }
  return(invisible())
}


image.map2 = function(lon, lat, z, center=0, hires=FALSE, add = FALSE, nlevel=1000, 
                      col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)),
                      land.col="darkolivegreen4", sea.col="aliceblue", boundaries.col = "black", 
                      grid.col="white", grid=FALSE, axes=TRUE, border=!axes) {
  
  image(x=lon, y=lat, z=z, col=col, axes=FALSE, add=add)
  
  mapDetails(center=center, hires=hires,col=land.col, interior=FALSE, 
             axes=axes, border=border, boundaries.col=boundaries.col,
             grid=grid, grid.col=grid.col)  
  return(invisible())
}


image.map = function (lon, lat, z, center=0, hires=FALSE, add = FALSE, nlevel = 1000, horizontal = FALSE, 
                      legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
                                                                                   3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
                      bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = rev(rainbow(nlevel/10, 
                                                                                               start = 0/6, end = 4/6)), 
                      lab.breaks = NULL, axis.args = NULL, legend.args = NULL, axes=TRUE,
                      midpoint = FALSE, border = NA, lwd = 1, land.col="black",
                      sea.col="white", boundaries.col="grey", grid=FALSE, grid.col="white", ...) 
{
  old.par <- par(no.readonly = TRUE)
  info <- imageplot.info(x=lon, y=lat, z=z, ...)
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(x=lon, y=lat, z=z, add = add, col = col, axes=FALSE, 
            xlab="", ylab="", ...)
      mapDetails(center=center, hires=hires,col=land.col, interior=FALSE, 
                 axes=axes, border=border, boundaries.col=boundaries.col,
                 grid=grid, grid.col=grid.col)    
    }
    else {
      poly.image(x=lon, y=lat, z=z, add = add, col = col, midpoint = midpoint, 
                 border = border, lwd.poly = lwd, axes=FALSE, 
                 xlab="", ylab="",...)
      mapDetails(center=center, hires=hires,col=land.col, interior=FALSE, 
                 axes=axes, border=border, boundaries.col=boundaries.col,
                 grid=grid, grid.col=grid.col)  
    }
    big.par <- par(no.readonly = TRUE)
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  ix <- 1
  minz <- info$zlim[1]
  maxz <- info$zlim[2]
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  breaks <- list(...)$breaks
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  if (!horizontal) {
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    }
    else {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
  }
  else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    }
    else {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
  }
  do.call("axis", axis.args)
  box()
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), line = legend.mar - 2)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    par(mar = par("mar"))
    invisible()
  }
  return(invisible())
}

.plotSea = function(col="lightblue", x0=-300, y0=-150) {
  x1 = -x0
  y1 = -y0
  polygon(x=c(x0,x1,x1,x0), y=c(y0,y0,y1,y1), border=NA, col=col)
  return(invisible())
}

.getDomain = function(domain, axis) { # create a modifiable environment
  
  if(!is.character(domain)) domain = deparse(substitute(domain))
  domains = NULL
  domains$peru      = list(x=c(-90,-70), y=c(-20,0)) 
  domains$peps      = list(x=c(-100,-70), y=c(-40,10))
  domains$peruS     = list(x=c(-75,-70), y=c(-20,-15))
  domains$peruN     = list(x=c(-86,-78), y=c(-10, -2))
  domains$peruC     = list(x=c(-80,-74), y=c(-16,-10))
  domains$peruNC    = list(x=c(-87,-73), y=c(-16,-2))
  domains$peru2     = list(x=c(-88,-70), y=c(-20,-2))
  domains$peru3     = list(x=c(-93,-70), y=c(-20,6))
  domains$ESPacific = list(x=c(-100,-65), y=c(-50,10))
  
  output = domains[[domain]][[axis]]
  return(output)
  
}

plot.map = function(x, y=NULL, xlim=NULL, ylim=NULL, domain=NA, center=0, 
                    hires=FALSE, land.col="darkolivegreen4", sea.col="aliceblue", 
                    boundaries.col = "black", grid.col="white", grid=TRUE,
                    cex=0.5, pch=19, main=NULL, add=FALSE, axes=TRUE, 
                    border=!axes, asp=NA, axs=NULL, xaxs=axs, yaxs=axs, ...) {
  
  if(is.null(xaxs)) xaxs = if(is.null(xlim)&is.na(domain)) "r" else "i"
  if(is.null(yaxs)) yaxs = if(is.null(ylim)&is.na(domain)) "r" else "i"
  
  if(missing(x)) x = NA
  if(is.data.frame(x) & is.null(y)) {
    names(x) = tolower(names(x))
    y = x$lat
    x = x$lon
  }
  if(is.character(x)) {
    domain = x
    x      = NA
  }
  xy = xy.coords(x, y)
  xy$xlab = ""
  xy$ylab = ""
  
  if(is.null(xlim)) xlim = .getDomain(domain, "x")
  if(is.null(xlim)) {
    xlim = range(xy$x, na.rm=TRUE)
    xaxs = "r"
  }
  
  if(is.null(ylim)) ylim = .getDomain(domain, "y")
  if(is.null(ylim)) {
    ylim = range(xy$y, na.rm=TRUE)
    yaxs = "r"
  }
  
  if(!add) {
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs, asp=asp)
    .plotSea(col=sea.col)
    mapDetails(center=center, hires=hires,col=land.col, interior=FALSE, 
               axes=axes, border=border, boundaries.col=boundaries.col,
               grid=grid, grid.col=grid.col)    
    title(main=main)
  }
  points(xy, cex=cex, pch=pch, ...)
  return(invisible())
}

plot.domain = function(domain=NA, xlim=NULL, ylim=NULL, col="red", lwd=1, ...) {
  
   
  if(is.null(xlim)) xlim = .getDomain(domain, "x")
  if(is.null(xlim)) stop("latitudinal domain limits has not been specified.")
  
  if(is.null(ylim)) ylim = .getDomain(domain, "y")
  if(is.null(ylim)) stop("longitudinal domain limits has not been specified.")

  xpol = c(xlim, rev(xlim))
  ypol = c(ylim[1],ylim[1],ylim[2],ylim[2])
  polygon(x=xpol, y=ypol, border=col, lwd=lwd, ...)
  
  return(invisible())
}

map.axes2 = function(cex.axis=0.75, line=-0.4) {
  .axis.map(1, "lon", las=1, cex.axis=cex.axis, line=line, tick=FALSE)
  .axis.map(2, "lat", las=1, cex.axis=cex.axis, line=line, tick=FALSE)
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  axis(3, labels=FALSE)
  axis(4, labels=FALSE)
  box()
  return(invisible(NULL))
}

.sub = function(map) {
  map = gsub(";","",gsub(";;",",",map))
  map = strsplit(map, ",")
  map = sapply(map, as.numeric)
  return(map)
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

rotate = function(x, direction="clockwise", revCol=TRUE) {
  out = switch(direction,
               clockwise     = .rotate(x, revCol=revCol),
               anticlockwise = .rotate2(x, revCol=revCol))
  return(out)
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

coord2text = function(coord, type) {
  if(!is.character(type)) type=deparse(substitute(type))
  out = sapply(coord, FUN=.coord2text, type=type)
  return(out)
}

.axis.map = function(side, type, ...) {
  x = axTicks(side=side)
  axis(side=side, at=x, labels=coord2text(coord=x, type=type), ...)
  return(invisible(x))
}

rename = function(data, old, new) {
  names(data)[names(data) == old] = new
  return(data)
}

table2grid = function(data, var, lat, lon, dx=dy, dy=dx, FUN=sum, ...) {
  
  FUN = match.fun(FUN)
  
  if(length(lat)< 2 | length(lon) < 2) stop("At lest a pair of lat and lon values must be provided.")
  
  latL = length(lat) 
  lonL = length(lon)
  
  if(latL == 2) {  
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
    latCut = coords$psi$lat
    cnames = coords$rho$lat
  } else {
    latCut = sort(lat)
    cnames = latCut[-1] - diff(latCut)
  }

  if(lonL == 2) {  
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
    lonCut = coords$psi$lon
    rnames = coords$rho$lon
  } else {
    lonCut = sort(lon)
    rnames = lonCut[-1] - diff(lonCut)
  }
  
  latAsFactor = cut(data[,"lat"], latCut, labels=FALSE)
  lonAsFactor = cut(data[,"lon"], lonCut, labels=FALSE)
  
  map = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor),
               FUN=FUN, ...)
  
  rows = seq_along(cnames)
  cols = seq_along(rnames)
  
  complete.rows = as.numeric(rownames(map))
  complete.cols = as.numeric(colnames(map))
  
  missing.rows = rows[!(rows %in% complete.rows)]
  missing.cols = cols[!(cols %in% complete.cols)]
  
  map = cbind(map, matrix(ncol=length(missing.cols), nrow=nrow(map)))
  map = rbind(map, matrix(ncol=ncol(map), nrow=length(missing.rows)))
  
  irows = sort(c(complete.rows, missing.rows), index.return=TRUE)$ix
  icols = sort(c(complete.cols, missing.cols), index.return=TRUE)$ix
  
  map = map[irows,icols]
  map = map[nrow(map):1, ] # correct lat direction
  map = .rotate(map) # write map for ncdf file, origin in down left corner
  
  rownames(map) = rnames
  colnames(map) = cnames
  
  return(map)
  
}


countMap = function(data, var, lat, lon, dx, dy=dx, ...) {
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  
  latAsFactor = cut(data[,"lat"], coords$psi$lat, labels=FALSE)
  lonAsFactor = cut(data[,"lon"], coords$psi$lon, labels=FALSE)
  
  map = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor),
               FUN=length)
  
  rows = seq_along(coords$rho$lat)
  cols = seq_along(coords$rho$lon)
  
  complete.rows = as.numeric(rownames(map))
  complete.cols = as.numeric(colnames(map))
  
  missing.rows = rows[!(rows %in% complete.rows)]
  missing.cols = cols[!(cols %in% complete.cols)]
  
  map = cbind(map, matrix(ncol=length(missing.cols), nrow=nrow(map)))
  map = rbind(map, matrix(ncol=ncol(map), nrow=length(missing.rows)))
  
  irows = sort(c(complete.rows, missing.rows), index.return=TRUE)$ix
  icols = sort(c(complete.cols, missing.cols), index.return=TRUE)$ix
  
  map = map[irows,icols]
  map = map[nrow(map):1, ] # correct lat direction
  map = .rotate(map) # write map for ncdf file, origin in down left corner
  
  rownames(map) = coords$rho$lon
  colnames(map) = coords$rho$lat
  
  image.map(coords$lon, coords$lat, map, ...)
  return(map)
  
}



table2array = function(data, var, lat, lon, start, end, 
                       dx, dy=dx, frequency=12, FUN=sum, toPA=FALSE) {
  
  FUN = match.fun(FUN)
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  latAsFactor  = cut(data[,"lat"], coords$psi$lat, labels=FALSE)
  lonAsFactor  = cut(data[,"lon"], coords$psi$lon, labels=FALSE)
  timeAsFactor = cut(data[,"time"], times$bounds, labels=FALSE)
  
  map0 = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor, timeAsFactor),
                FUN=FUN, na.rm=TRUE)
  
  rows    = seq_along(coords$rho$lat)
  cols    = seq_along(coords$rho$lon)
  slices  = seq_along(times$center)
  
  complete.rows   = as.numeric(rownames(map0))
  complete.cols   = as.numeric(colnames(map0))
  complete.slices = as.numeric(dimnames(map0)[[3]])
  
  missing.rows   = rows[!(rows %in% complete.rows)]
  missing.cols   = cols[!(cols %in% complete.cols)]
  missing.slices = slices[!(slices %in% complete.slices)]
  
  map = array(dim=c(length(rows), length(cols), length(slices)))
  map[seq_along(complete.rows), seq_along(complete.cols), 
      seq_along(complete.slices)] = map0 
  
  irows = sort(c(complete.rows, missing.rows), index.return=TRUE)$ix
  icols = sort(c(complete.cols, missing.cols), index.return=TRUE)$ix
  itime = sort(c(complete.slices, missing.slices), index.return=TRUE)$ix
  
  map = map[irows,icols, itime]
  
  map = map[nrow(map):1, , ] # correct lat direction
  map = .rotate(map) # write map for ncdf file, origin in down left corner
  
  dimnames(map) = list(coords$rho$lon, coords$rho$lat, round(times$center,2))
  
  if(toPA) map=.toPA(map)
  
  return(map)
  
}

rotateMap =  function(x) {
  # Check this function!!
  out = x$lon>0
  out[is.na(out)] = FALSE
  x$lon[out] = x$lon[out] - 360
  return(x)
}

.toPA = function(x) {
  return(!!x + 0)
}

toPA= function(x, thr=0, prob=FALSE) {
  out = x
  out[x>thr] = 1 
  out[x<=thr] = 0
  if(prob) out = out*x
  return(out)
}

toPAm= function(x, thr=c(0,0.5)) {
  thr = sort(thr)
  out = x
  out[x>thr[1]]  = 0.5 
  out[x>thr[2]]  = 1 
  out[x<=thr[1]] = 0
  return(out)
}

map2data = function(z, x, y) {
  # Transform a matrix-like object z to a data.frame
  # object
  # z is a matrix of dimensions length(x) times
  # length(y) 
  if(nrow(z) !=length(x) || ncol(z)!=length(y)) {
    stop("dimension of matrix z must be length(x) times
         # length(y)")
  }
  X = matrix(x, ncol=length(y), nrow=length(x))
  Y = matrix(y, ncol=length(y), nrow=length(x), byrow=T)
  
  X = .array2vector(X)
  Y = .array2vector(Y)
  Z = .array2vector(Z)
  
  out = data.frame(x=X, y=Y, z=z)
  return(out)
}


.array2vector = function(z) {
  dim(z) = c(prod(dim(z)),1)
  return(z)
}

map2 =  function(database,center,...){
  # From stackoverflow
  Obj <- map(database,...,plot=FALSE)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=TRUE)
  polygons <- apply(id, 1, function(i) {coord[i[1]:i[2],]} )
  
  # split up polygons that differ too much
  polygons <- lapply(polygons, function(x) {
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    return(x)
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

extractValidData = function(files, control, var="traj", 
                            output=NULL) {
  # files is a set of ncdf files paths
  # extract data from variables with the same dimension of control, and
  # matrix with the same lat and lon dimensions. 
  cat("Reading control variable...\n")
  control = open.ncdf(control)
  control.var = get.var.ncdf(control, var)
  control.dim = dim(control.var)
  ix = which(!is.na(control.var))
  close.ncdf(control)
  gc(verbose=FALSE)
  base = NULL
  for(file in files) {
    cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      isValidArray  = identical(data1$var[[var]]$varsize, control.dim)
      isValidMatrix = identical(data1$var[[var]]$varsize, control.dim[1:2])
      if( !isValidArray & !isValidMatrix ) {
        cat("Skipping variable", var,"\n")  
        next
      } else {
        
        cat("Adding variable", var,"\n")
        base.names = colnames(base)
        if(isValidMatrix) {        
          base = cbind(base, 
                       array(get.var.ncdf(data1, var), dim=control.dim)[ix])
        }
        if(isValidArray) {
          base = cbind(base, get.var.ncdf(data1, var)[ix])
        }
        colnames(base) = c(base.names, var)
        
      }
      
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  cat("Writing final data base.\n")
  base = data.frame(base)
  if(!is.null(output)) write.csv(base, output)
  return(base)
}


ncdf2data = function(files, slices, control, var, 
                     control.dim=NULL) {
  # files is a set of ncdf files paths
  # extract data from variables with the same dimension of control, and
  # matrix with the same lat and lon dimensions. 
  cat("Reading control variable...\n")
  if(is.null(control.dim)) {
    control = open.ncdf(control)
    control.dim = control$var[[var]]$varsize
    cat("Extracting data with dimensions", control.dim, "\n")
    close.ncdf(control)
    gc(verbose=FALSE)    
  }
  base = NULL
  for(file in files) {
    cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      isValidArray  = identical(data1$var[[var]]$varsize, control.dim)
      isValidMatrix = identical(data1$var[[var]]$varsize, control.dim[1:2])
      if( !isValidArray & !isValidMatrix ) {
        cat("Skipping variable", var,"\n")  
        next
      } else {
        cat("Adding variable", var,"\n")
        base.names = colnames(base)
        if(isValidMatrix) {
          newvar = array(get.var.ncdf(data1, var), dim=control.dim)[,,slices]
        }
        if(isValidArray) {
          newvar = get.var.ncdf(data1, var)[,,slices]
        }
        dim(newvar) = NULL
        base = cbind(base, newvar)
        colnames(base) = c(base.names, var)        
      }
      
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  cat("Writing final data base.\n")
  base = data.frame(base)
  return(base)
}


extractData = function(files, data, lat, lon, start, end, 
                       dx, dy=dx, frequency=12, 
                       dim.names = c("lat", "lon", "time")) {
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  latAsFactor  = cut(data[, dim.names[1]], coords$psi$lat, labels=FALSE)
  lonAsFactor  = cut(data[, dim.names[2]], coords$psi$lon, labels=FALSE)
  timeAsFactor = cut(data[, dim.names[3]], times$bounds, labels=FALSE)
  
  ix  = cbind(lonAsFactor, latAsFactor, timeAsFactor)
  ix2 = cbind(lonAsFactor, latAsFactor)
  
  control.dim = c(length(coords$rho$lon),
                  length(coords$rho$lat),
                  length(times$center))                  
  
  for(file in files) {
    cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      if(!identical(data1$var[[var]]$varsize, control.dim)) {
        if(identical(data1$var[[var]]$varsize, control.dim[1:2])) {
          cat("Adding variable", var,"\n")
          data[, var] = get.var.ncdf(data1, var)[ix2]  
        } else {
          cat("Skipping variable", var,"\n")  
          next  
        }
      } else {
      cat("Adding variable", var,"\n")
      data[, var] = get.var.ncdf(data1, var)[ix]
      }
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  cat("Writing final data base.\n")
  return(data)
}


fill.map = function(x, mask=1, FUN="max") {
  # fill a map with the max (min, mean, median) value of the map
  FUN = match.fun(FUN)
  maximo = FUN(x, na.rm=TRUE)
  x[is.na(x)] = maximo
  x = x*mask
  return(x)
}

convexHull.map = function(x, mask, excludeZeros=TRUE) {
  
  if(all(is.na(x))) {
    return(x)
  } else {
    ind = if(excludeZeros) (!!x & !is.na(x)) else !is.na(x) 
    xc = col(x)[ind]
    xr = row(x)[ind]
    coord = cbind(xr, xc)
    xvar = x[coord]
    
    x.int = interp(x=xr, y=xc, z=xvar, 
                   xo=seq(len=nrow(x)), yo=seq(len=ncol(x)),
                   duplicate="strip", linear=TRUE,
                   extrap=FALSE)$z * mask
    x.int[mask==0] = NA
    return(x.int)    
    
  }
}

convexHull = function(x, mask) {
  dim.old = dim(x)
  out = apply(x, 3, convexHull.map, mask=mask)  
  dim(out) = dim.old
  return(out)
}

regridMap = function(old, new, normalize=TRUE, ...) {
  # generalize to any kind of input and output data
  
  new$mask[new$mask==0] = NA # correct mask
  
  ind = !apply(old, 1, function(x) any(is.na(x)))
  old = old[ind, ] # remove NAs
  
  newp   = interpp(old$lon, old$lat, old$z, 
                   xo=as.numeric(new$lon), 
                   yo=as.numeric(new$lat),
                   ...)$z
  
  newmap = matrix(newp, 
                  ncol=ncol(new$mask), 
                  nrow=nrow(new$mask))
  newmap[is.na(newmap)] = 0
  newmap = newmap*new$mask
  if(normalize) newmap = newmap/sum(newmap, na.rm=TRUE)
  
  return(newmap)
}

table2ncdf = function(file, species, 
                      lat, lon, dx, dy,
                      start, end, frequency=12,
                      toPA=FALSE, longNames=NULL,
                      FUN=mean, colClasses=NA) {
  
  
  cat("Reading data from file", paste0(file,"...", "\n"))
      
  precision = if(toPA) "integer" else "float"
  units     = if(toPA) "1/0" else "abundance"
  unitsName = if(toPA) "Presence/Absence" else "Abundance"

  output = sub(pattern=".csv",".nc", file)
  if(is.null(longNames)) longNames = paste(species, unitsName)
      
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  dim1 = dim.def.ncdf("lon",  "degrees", coords$rho$lon)
  dim2 = dim.def.ncdf("lat", "degrees", coords$rho$lat)
  dim3 = dim.def.ncdf("time", "months", times$center)
      
  biology = read.csv(file=file, colClasses=colClasses)
  
  checkNames = all(c("year", "month") %in% names(biology))
  if(!checkNames) stop(paste("Input file has no", 
                             sQuote(year), "or", sQuote(month),
                             "fields."))
  
  biology$time = biology$year + biology$month/12 - 1/24
  biology$traj = 1
  
  species = c(species, "traj")
  longNames = c(longNames, "Survey trajectory")
  
  cat("Creating maps from data...\n")
          
  species.data = list()
  ncdf.vars = list()
      
  for(var in seq_along(species)) {
    cat("Species:", paste0(species[var],"... "))
    species.data[[var]] = table2array(biology, 
                                      var=species[var], 
                                      lat=lat, lon=lon, 
                                      start=start, end=end, 
                                      dx=dx, FUN=FUN, 
                                      toPA=toPA)
    
    ncdf.vars[[var]]    = var.def.ncdf(name=species[var], units=units, 
                                       dim=list(dim1,dim2,dim3), missval=-9999, 
                                       longname=longNames[var], 
                                       prec=precision)
    
    cat("DONE.\n")
  }
  
  cat("Creating ncdf file...\n")
      
  nc.new = create.ncdf(output, ncdf.vars)

  for(var in seq_along(species)) {
    # put values to vars
    cat("\t Writing", paste0(species[var],"...\n"))
    put.var.ncdf(nc.new, ncdf.vars[[var]], species.data[[var]])
  }

  cat("Saving data.\n")
      
  close.ncdf(nc.new)
      
  cat("File", output, "has been written.","\n")     
  
  return(invisible(output))
  
}


trainingDatabase = function(files, control, var="traj", 
                                output=NULL) {
  extractValidData(files=files, control=control, var=var, 
                              output=output)  
}

  
writePredictionData = function(files, control, var="traj", path,
                      start, end, frequency=12) {

  times  = createTimeAxis(start=start, end=end, frequency=frequency, 
                          center=TRUE)$center
  
  control = open.ncdf(control)
  control.dim = control$var[[var]]$varsize
  cat("Extracting data with dimensions", control.dim, "\n")
  close.ncdf(control)
  gc(verbose=FALSE)
  
  n = ceiling(log(length(times)))

  for(t in seq(along=times)) {
    
    cat("----------------------\n")
    cat("----- time step =", t, "\n")
    cat("----------------------\n")
    base = ncdf2data(files = files, slice=t,
                     control = control, 
                     var = var,
                     control.dim=control.dim)
    
  file = file.path(path, sprintf(paste0("base_env_%0",n,"d.csv"),t))
  write.csv(base, file)
    
  }
  
  return(invisible(control.dim[1:2]))
  
}

setDataBase = function(base, factors) {
  base = setFactors(base, factors)
  base = base[complete.cases(base),]  
  return(base)
}

setFactors = function(base, factors) {
  for(f in factors) {
    base[, f] = as.factor(base[, f])
  }
  return(base)
}


makePredictions = function(input, output, run, factors) {
  
  root = getwd()
  
  pred.files = dir(path=input)
  
  for(t in seq_along(pred.files)) {
    cat("Running map ", t, ": ", as.character(Sys.time()),"\n", sep="")
    data.file = file.path(input, pred.files[t])
    newdata = read.csv(data.file, row.names=1)
    newdata = setFactors(newdata, factors)
    complete = complete.cases(newdata)
    new.data = newdata[complete,]
    
    setwd(run)
    Projection(Proj=new.data, Proj.name=t, 
               BinRoc=TRUE, BinKappa=TRUE, BinTSS=TRUE)
    
    dir.complete = file.path(getwd(),paste0("proj.", t))
    file.complete = paste("Proj_", t, "_complete", sep="")
    save(complete, file=file.path(dir.complete, file.complete))
    
    setwd(root)
    
    file.copy(from=file.path(dir.complete, dir(path=dir.complete)), 
              to=file.path(output,dir(path=dir.complete)), overwrite=TRUE)
    
    unlink(dir.complete, recursive=TRUE, force=TRUE)
    
  }
  
}


fitModels = function(run, model.session,
                     GLM = FALSE, TypeGLM = "simple", Test = "AIC", GBM = FALSE, 
                     No.trees = 5000, GAM = FALSE, Spline = 3, CTA = FALSE, CV.tree = 50, 
                     ANN = FALSE, CV.ann = 5, SRE = FALSE, quant = 0.025, FDA = FALSE, 
                     MARS = FALSE, RF = FALSE, NbRunEval = 1, DataSplit = 100, 
                     NbRepPA = 0, strategy = "sre", coor = NULL, distance = 0, 
                     nb.absences = NULL, Yweights = NULL, VarImport = 0, Roc = FALSE, 
                     Optimized.Threshold.Roc = FALSE, Kappa = FALSE, TSS = FALSE, 
                     KeepPredIndependent = FALSE) {
  
  cat("Running BIOMOD.\n")
  cat("Starting at", as.character((Sys.time())), "\n")
  
  root = getwd()
  on.exit(setwd(root))
  setwd(run)
  
  Models(GLM = GLM, TypeGLM = TypeGLM, Test = Test, GBM = GBM, 
         No.trees = No.trees, GAM = GAM, Spline = Spline, CTA = CTA, CV.tree = CV.tree, 
         ANN = ANN, CV.ann = CV.ann, SRE = SRE, quant = quant, FDA = FDA, 
         MARS = MARS, RF = RF, NbRunEval = NbRunEval, DataSplit = DataSplit, 
         NbRepPA = NbRepPA, strategy = strategy, coor = coor, distance = distance, 
         nb.absences = nb.absences, Yweights = Yweights, VarImport = VarImport, 
         Roc = Roc, Optimized.Threshold.Roc = Optimized.Threshold.Roc, 
         Kappa = Kappa, TSS = TSS, KeepPredIndependent = KeepPredIndependent)
  
  cat("Ended at", as.character((Sys.time())), "\n")
  save.image(model.session)
  
  setwd(root)
  
  return(invisible())
}


getPredictionMap = function(output, time, species) {
  root = getwd()
  setwd(file.path(root,output))
  i = time
  proy = paste0("Proj_", i, "_",species)
  comp = paste0("Proj_", i, "_complete")
  load(proy)
  load(comp)
  proy = get(proy)
  pred = array(dim=c(pred.dim, 9))
  pred[complete] = as.vector(proy[,,1,1])/1000 
  return(pred)
}


.getDistance = function(xy, ref, ref2abs) {
  
  .getDist = function(xy, ref, ref2abs) {
    # calculate the distance from xy to ref. 
    # Distance is always positive.
    # ncol(ref) = 4, length(xy)=4, length(ref2abs) = nrow(ref)
    grados=acos(ref[,4]*xy[4] + ref[,3]*xy[3]*cos(ref[,2]-xy[2]))
    ind  = which.min(grados)
    dist = c(grados[ind], ref2abs[ind]) 
    return(dist)
  }
  
  out = t(apply(xy, 1, .getDist, ref=ref, ref2abs=ref2abs))
  
  return(out)
}

.getDistance.old = function(data, ref, ref2abs=NULL, lon="lon", lat="lat") {
  
  x    = as.matrix(data[,c(lat, lon)])
  cx   = as.matrix(ref[,c(lat, lon)])
  
  grados2km  = 1.852*60*180/pi
  grados2rad = pi/180
  
  if(!is.null(ref2abs)) ref2abs = ref2abs/grados2km
  
  ref = cx*grados2rad
  ref = cbind(ref, cos(ref[,1]), sin(ref[,1])) # c(lat, lon, cos(lat), sin(lat))
  xy  = x*grados2rad  
  xy  = cbind(xy, cos(xy[,1]), sin(xy[,1]))
  
  output = array(dim=c(nrow(data), 1 + !is.null(ref2abs)))
  ind = complete.cases(xy)
  output[ind,] = .getDistance(xy=xy[ind, ], ref=ref, ref2abs=ref2abs)*grados2km
  return(output)
}


getDistance = function (data, ref, ref2abs = NULL, lon = "lon", lat = "lat") 
{
  data     = data[,c(lat,lon)]
  dat      = data[!duplicated(data),]
  data$ind = 1:nrow(data)
  
  x = as.matrix(dat[, c(lat, lon)])
  cx = as.matrix(ref[, c(lat, lon)])
  grados2km = 1.852 * 60 * 180/pi
  grados2rad = pi/180
  if (!is.null(ref2abs)) ref2abs = ref2abs/grados2km
  ref = cx * grados2rad
  ref = cbind(ref, cos(ref[, 1]), sin(ref[, 1]))
  xy = x * grados2rad
  xy = cbind(xy, cos(xy[, 1]), sin(xy[, 1]))
  output = array(dim = c(nrow(dat), 1 + (!is.null(ref2abs))))
  ind = which(complete.cases(xy))
  output[ind, ] = .getDistance(xy = xy[ind, ], ref = ref, ref2abs = ref2abs) * 
    grados2km
  xnames = if(!is.null(ref2abs)) c("dist", "ref2abs") else "dist"
  colnames(output) = xnames
  output = data.frame(dat, output)
  output = merge(data, output, all=TRUE, sort=FALSE)
  ind    = sort(output$ind, index=TRUE)$ix
  output = as.matrix(output[ind, xnames])
  return(output)
}

getSignedDistance = function(data, ref, abs, lon="lon", lat="lat") {

  grados2km  = 1.852*60*180/pi
  
  DateStamp("Calculating distance from data to abs.\t")
  data2abs = getDistance(data=data, ref=abs, lon=lon, lat=lat)
  DateStamp("Calculating distance from ref to abs.\t")
  ref2abs  = getDistance(data=ref, ref=abs, lon=lon, lat=lat)
  DateStamp("Calculating distance from data to ref.\t")
  data2ref = getDistance(data=data, ref=ref, ref2abs=ref2abs, lon=lon, lat=lat)
  
  sign = -2*(data2abs < data2ref[,2]) + 1 # inside is negative
  
  dist = sign*data2ref[,1]
  
  output = list(dist=dist, abs=data2abs)
  return(output)
}

read.Herve = function(file, skip=9, units=TRUE, k=1000, output=".") {
  
  .parseLine = function(line) {
    x = unlist(strsplit(line," "))
    x = x[x!=""]
    return(x)
  }
  
  .parseHeader = function(line) {
    cat("\t Parsing header...\n")
    var = paste(strsplit(file, "\\.")[[1]][c(5:4)],collapse=".")
    var = paste0(var,".")
    header = .parseLine(dat[skip+1])
    header = sub("\\.", "",tolower(header))
    header = sub("longitude", "lon", header)
    header = sub("latitude", "lat", header)
    header = sub("val", "aver", header)
    header = sub("stddev", "sd", header)
    ind = header=="aver"
    rmax = sum(ind)
    header[ind] = paste0(var, "avg_r", seq_len(rmax)-1)
    ind = header %in% c("min", "max", "sd", "mode")
    header[ind] = paste0(var, header[ind], "_r", rmax-1)
    return(header)
  }
  
  .readHerve = function(line, ncols) {
    x = .parseLine(line)
    if(length(x)!=ncols) {
      out = c(x[2:6], rep(NA, ncols-7))
    } else {
      x[x=="********"] = NA
      x[is.nan(x)] = NA
      out = x[-c(1,ncols)]
    }
    out = as.numeric(out)
    names(out) = NULL
    return(out)
  }
  
  cat("\t Reading file", file, "\n")
  dat = readLines(file)
  header = .parseHeader(dat[skip+1])
  ncols = length(header)
  cat("\t\t",header,"\n")
  if(units) unit = .parseLine(dat[skip+1+units])
  dat = dat[-seq_len(skip+1+units)]
  
  cat("\t Parsing data...\n")
  x4 = array(dim=c(length(dat), ncols-2))
  pp = as.integer(seq(from=0, to=nrow(x4), len=k+1))
  
  for(j in seq_len(k)) {
    ind = seq.int(pp[j]+1, pp[j+1])
    x4a = lapply(dat[ind], .readHerve, ncols=ncols)
    x4[ind,] = t(matrix(unlist(x4a), nrow=ncols-2))
  }
  
  cat("\t Creating data frame...\n")
  x4 = as.data.frame(x4)
  names(x4) = header[-c(1,ncols)]
  cat("\t Done.\n")
  
  file.out = paste(strsplit(file, "\\.")[[1]][c(1,4:6,8,2)],collapse=".")
  file.out = file.path(output, file.out)
  cat("\t Saving", sQuote(file.out), "...\n")
  write.csv(x4, file.out, row.names=FALSE)
  
  return(invisible(file.out))
}

DateStamp = function(...) cat(..., "\t\t | ", date(), "\n\n")


getIsoline = function(lon, lat, z, level, cutoff=0) {
  .getLine = function(x) rbind(cbind(x$x, x$y), NA)
  out = grDevices::contourLines(x=lon, y=lat, z=z, nlevels=1, levels=level)
  ll = unlist(lapply(out, function(x) length(x$x)))
  ns = length(ll)
  ll2 = sort(ll)
  cat("Creating ", ns, " isolines of lengths ", paste(ll2[-ns], collapse=", "),
      " and ", ll2[ns], ".\n", sep="")
  cat("Removing ", sum(ll<=cutoff)," isolines of length lower than cutoff=", cutoff, ".\n", sep="")
  out = out[ll>cutoff]
  out = do.call(rbind, lapply(out, FUN=.getLine))
  out = out[-nrow(out),]
  out = as.data.frame(out)
  names(out) = c("lon", "lat")
  return(out)
}

AUC = function(data, coordNames = c("lat", "lon"), obs="observed",
               models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  output = NULL
  for(i in seq_len(sum(models))) {
    output = rbind(output, auc(DATA,st.dev=st.dev, which.model=i,na.rm=na.rm))
  }
  rownames(output) = names(data)[models] 
  return(output)
}


kappa = function(data, coordNames = c("lat", "lon"), obs="observed",
                 models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  
  thr = as.numeric(optimal.thresholds(DATA=DATA, opt.methods="MaxKappa")[-1])
  output = NULL
  for(i in seq_len(sum(models))) {
    output = rbind(output, Kappa(cmx(DATA,threshold=thr[i], which.model=i,
                                     na.rm=na.rm), st.dev=st.dev))
  }
  rownames(output) = names(data)[models]
  return(output)
}


TSS = function(data, coordNames = c("lat", "lon"), obs="observed",
               models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  
  thr = as.numeric(optimal.thresholds(DATA=DATA, opt.methods="MaxSens+Spec")[-1])
  output = NULL
  for(i in seq_len(sum(models))) {
    CMX = cmx(DATA,threshold=thr[i], which.model=i, na.rm=na.rm)
    output = rbind(output, cbind(specificity(CMX, st.dev=st.dev), 
                                 sensitivity(CMX, st.dev=st.dev)))
  }
  if(st.dev) {
    tss = cbind(output[,1]+output[,3]-1, sqrt(output[,2]^2+output[,4]^2))
    colnames(tss) = c("TSS", "TSS.sd")
  } else {
    tss = cbind(output[,1]+output[,2]-1)
    colnames(tss) = c("TSS")
  }
  output = cbind(tss, output)
  rownames(output) = names(data)[models]
  return(output)
}

PredictivePerformance = function(data, coordNames = c("lat", "lon"), obs="observed",
                                 models=NULL, st.dev=TRUE, na.rm=TRUE) {
  out1 = AUC(data=data, coordNames = coordNames, obs=obs, 
             models=models, st.dev=TRUE, na.rm=na.rm)
  out2 = kappa(data=data, coordNames = coordNames, obs=obs, 
               models=models, st.dev=TRUE, na.rm=na.rm)
  out3 = TSS(data=data, coordNames = coordNames, obs=obs, 
             models=models, st.dev=TRUE, na.rm=na.rm)
  
  output1 = cbind(AUC=out1[,1], Kappa=out2[,1], out3[,c(1,3,5)])
  output2 = cbind(AUC.sd=out1[,2], Kappa.sd=out2[,2], out3[,c(2,4,6)])
  
  output = if(st.dev) {
    list(statistic=output1, sd=output2)
  } else {
    output1
  }
  return(output)
  
}

plotROC = function(data, coordNames = c("lat", "lon"), obs="observed",
                   models=NULL, opt.thresholds=TRUE, opt.methods=c(4,10),
                   req.sens=0.95, ...) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  auc.roc.plot(DATA=DATA, opt.thresholds=opt.thresholds, 
               opt.methods=opt.methods, req.sens=req.sens, ...)
  return(invisible())
}

plotThreshold = function(data, coordNames = c("lat", "lon"), obs="observed",
                         models=NULL, opt.thresholds=TRUE, opt.methods=c(4,10),
                         req.sens=0.95, which.model=1, ...) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  error.threshold.plot(DATA=DATA, which.model=which.model, opt.thresholds=opt.thresholds, 
                       opt.methods=opt.methods, req.sens=req.sens, ...)
  return(invisible())
}

calculateThresholds = function(data, coordNames = c("lat", "lon"), 
                               obs="observed", 
                               models=NULL, opt.methods=2:12,
                               req.sens=0.95, req.spec=0.5,
                               FPC=1, FNC=10, ...) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  out = optimal.thresholds(DATA=DATA, opt.methods=opt.methods, 
                           req.sens=req.sens, req.spec=req.spec,
                           FPC=FPC, FNC=FNC, ...)
  output = as.matrix(out[,-1])
  rownames(output) = out[,1]
  return(output)
}

plotMAPt = function(lon, lat, map, thr, ...) {
  image.map(lon=lon, lat=lat, 0 + map>thr, ...)
}

cleanBIOMOD = function() {
  if(exists("Biomod.material", env=globalenv())) rm(list="Biomod.material", envir=globalenv())
  if(exists("DataBIOMOD", env=globalenv())) rm(list="DataBIOMOD", envir=globalenv())
  if(exists("DataEvalBIOMOD", env=globalenv())) rm(list="DataEvalBIOMOD", envir=globalenv())
  return(invisible())
}

is.inside = function(x, range) {
  ind = (x >= range[1]) & (x <= range[2])
  return(ind)
}

geq = function(x, thr) {
  ind = which(x>=thr & !is.na(x))
  out = x[ind]
  names(out) = ind
  return(out)  
}

leq = function(x, thr) {
  ind = which(x<=thr & !is.na(x))
  out = x[ind]
  names(out) = ind
  return(out)  
}

.insider = function(data, control, var, alpha, lowerOnly) {
  # TODO: handle factors
  delta = alpha/2
  alpha = if(!lowerOnly) c(delta, 1 - delta) else c(1-alpha, 1)
  alpha = quantile(control[,var], prob=alpha, na.rm=TRUE)
  ind = (data[,var] > alpha[1]) & (data[,var] < alpha[2])
  return(ind)
}

clearAbsences = function(data, species, vars, alpha=0.05, lowerOnly=FALSE) {
  dataP    = data[data[,species]==1,]
  dataN    = data[data[,species]==0,]
  dataN_PA = cleanZeros(data=dataN, control=dataP, vars=vars, 
                        alpha=alpha, lowerOnly=lowerOnly)
  data_PA  = rbind(dataP, dataN_PA)
  return(data_PA)
}

cleanZeros = function(data, control, vars, alpha=0.05, lowerOnly=FALSE) {
  
  if(length(alpha)==1) alpha = rep(alpha, len=length(vars))
  if(length(alpha)!=length(vars)) stop("You must provide as many confidence thresholds as control vars.")
  if(length(lowerOnly)==1) lowerOnly = rep(lowerOnly, len=length(vars))
  if(length(lowerOnly)!=length(vars)) stop("You must provide as many confidence thresholds type as control vars.")
  ind = rep(TRUE, len=nrow(data))
  for(i in seq_along(vars)) {
    ind = ind & .insider(data, control, vars[i], alpha[i], lowerOnly[i])
  }
  ind = !ind & !is.na(ind)
  output = data[ind, ]
  return(output)
}

.getmfrow = function(n) {
  m1 = floor(sqrt(n))
  m2 = ceiling(n/m1)
  out = rev(sort(c(m1, m2)))
  return(out)
}

getmfrow = function(n) .getmfrow(n=n)

table.summary = function(..., sp) {
  names  = as.character(substitute(list(...))[-1L])
  object = list(...)
  output = NULL
  tnames = NULL
  for(i in seq_len(length(object))) {
    tnames = c(tnames, levels(as.factor(object[[i]][,sp])))
  }
  tnames = unique(tnames)
  output = array(0, dim=c(length(object), length(tnames)))
  colnames(output) = tnames
  rownames(output) = names
  for(i in seq_len(length(object))) {
    tnames = levels(as.factor(object[[i]][,sp]))
    output[i,tnames] = table(object[[i]][,sp])
  }
  total = apply(output, 1, sum, na.rm=TRUE)
  percent = 100*round(output/total,3)
  output = cbind(output, total, percent)
  return(output)
}

balancePA = function(data, PA, sp, rpa=0.1, seed=771104) {
  if(rpa<0 | rpa>1) error("rpa must be in [0,1].")
  ind.pa = which(PA[,sp]==0 & !is.na(PA[,sp]))
  ind.a  = which(data[,sp]==0 & !is.na(data[,sp]))
  ind.p  = which(data[,sp]==1 & !is.na(data[,sp]))
  
  np  = length(ind.p)
  nn  = length(ind.a)
  npa = max(np - nn, 0)
  ppa = npa/np
  
  if(ppa >= rpa) {
    set.seed(seed)
    ind.pa = sample(ind.pa, npa)
  } else {
    npa = floor(rpa*np)
    set.seed(seed)
    ind.pa = sample(ind.pa, npa)
    ind.a  = sample(ind.a, np-npa)
  }
  cat("PA ratio = ", npa/np, " (critic = ", rpa,")\n", sep="")
  print( c(presence=np, absence=np-npa, PA=npa) )
  output = rbind(data[ind.p,], 
                 data[ind.a,],
                 PA[ind.pa,])
  return(output)
}

splitDataSet = function(data, var, factor=0.15, seed=771104) {
  if(class(data)=="niche.models") {
    output = list(train=data$train, val=data$val)
  }
  if(is.data.frame(data)) {
    data = data[!is.na(data[,var]),]
    dataN = data[data[,var]==0,]
    dataP = data[data[,var]==1,]
    nn = nrow(dataN)
    np = nrow(dataP)
    set.seed(seed)
    indn = sample(nn, factor*nn)
    indp = sample(np, factor*np)
    training = rbind(dataN[-indn,], dataP[-indp,])
    validation = rbind(dataN[indn,], dataP[indp,])
    output = list(train=training, val=validation)
  }
  class(output) = c("niche.models", class(output))
  return(output)    
}

.gam.fmla = function(y, var, spline=NULL, loess=NULL, factors=NULL) {
  
  if(is.null(loess)) {
    loess = rep(FALSE, length(var))
    names(loess) = var
  }
  
#   ss = c(sst=2, sss=2, chl=2)
#   lo = c(sst=FALSE, sss=FALSE, chl=TRUE)
#   
#   ft = rep("s", length(var))
#   ft[lo] = "lo"
#   lt = rep(",%02d", length(var))
#   lt[lo] = ""
#   x = paste0(ft, "(",var, lt, ")")
#   x = paste(x, collapse=" + ")
  
  if(!is.null(spline)) {
    fmla = as.list(spline[var])
    fmla = append(fmla, paste0("s(", var, ",%d)", collapse=" + "), after=0)
    fmla = do.call(sprintf, fmla)
  } else {
    fmla = paste0(var, collapse=" + ")
  }
  
  if(!is.null(factors)) {
    
   fmla = paste0(c(fmla, factors), collapse=" + ") 
   
  }
  
  fmla = paste0(y, " ~ ",fmla)
  return(fmla)
}


gam.fmla = function(y, var, spline=NULL, factors=NULL) {
  fmla = .gam.fmla(y=y, var=var, spline=spline, factors=factors)
  fmla = as.formula(fmla, env=.GlobalEnv)
  return(fmla)
}

.fmla2txt = function(fmla) {
  out = paste(as.character(fmla)[c(2,1,3)], collapse=" ")
  return(out)
}

.getObjectName = function(x) {
  out = deparse(substitute(x))
  out = deparse(eval(substitute(quote(x))))
  return(out)
}

fakeGAM = function(object) {
  
  obj = deparse(eval(substitute(quote(object), env=as.environment(-1))))
  DateStamp("Fitting models for", sQuote(obj), "dataset.")
  return(obj)
}

fakeGAM2 = function(object) {
  
  obj = deparse(substitute(object, env=as.environment(-1)))
  DateStamp("Fitting models for", sQuote(obj), "dataset.")
  return(obj)
}

fitGAMs = function(object, formulas, FUN=identity, 
                   name=deparse(substitute), link="logit") {
  
  FUN = match.fun(FUN)
  object$transform = FUN
  
  gc(verbose=FALSE)
  
#   obj1 = deparse(eval(substitute(quote(object))))
#   obj2 = deparse(substitute(object, env=as.environment(-1)))
#   
  DateStamp("Fitting models for", sQuote(name), "dataset.")
#   DateStamp("Fitting models for", sQuote(obj2), "dataset.")
#   fakeGAM(object=object)
#   fakeGAM2(object=object)

  train = FUN(object$train)
  val   = FUN(object$val)
  
  if(!is.list(object$predicted)) object$predicted = NULL
  
  models = names(formulas)
  
  aic = object$fit[,"AIC"]
  bic = object$fit[,"BIC"]
  var = as.character(formulas[[1]])[2]
  
  object$predicted$lon = object$train[,"lon"]
  object$predicted$lat = object$train[,"lat"]
  object$predicted$observed = as.numeric(as.character(object$train[,var]))
  object$predicted = as.data.frame(object$predicted)
  
  object$validation$lon = object$val[,"lon"]
  object$validation$lat = object$val[,"lat"]
  object$validation$observed = as.numeric(as.character(object$val[,var]))
  object$validation = as.data.frame(object$validation)
  
  for(i in seq_along(formulas)) {
    model.name = models[i]
    model.formula = formulas[[i]]
    DateStamp("Training model ", model.name, ":\n", .fmla2txt(formulas[[i]]), sep="")
    object$formulas[[model.name]] = model.formula 
    #model.vars = .getModelVars2(model.formula, train)
    # TO_DO: filter complete cases
    model = gam(model.formula, family = binomial(link=link), data = train)
    gc(verbose=FALSE)
    model$anova = anova(model)
    object$models[[model.name]] = model
    # object$preplot[[models[i]]] = with(object$train, preplot(model))
    
    aic[model.name] = AIC(model)
    bic[model.name] = BIC(model)
    object$predicted[, model.name] = model$fitted
    object$validation[, model.name] = predict(model, newdata=val, type="response")
  }
  
  DateStamp("Computing Predictive Performance...")
  
  object$fit = cbind(AIC=aic, BIC=bic)
  
  object$performance$training   = PredictivePerformance(object$predicted, st.dev=FALSE)
  object$performance$validation = PredictivePerformance(object$validation, st.dev=FALSE)
  
  object$threshold$training  = calculateThresholds(object$predicted)
  object$threshold$validation = calculateThresholds(object$validation)
  
  DateStamp("DONE.")
  
  class(object$formulas) = c("niche.models.formulas", class(object$formulas))
  return(object)
}


getPredictions = function(data) {
  if(is.list(data$predicted)) return(do.call(cbind, data$predicted)) else return(invisible(NULL))
}

print.niche.models.formulas = function(x, ...) {
  out = matrix(sapply(x, .fmla2txt), ncol=1)
  colnames(out) = ""
  rownames(out) = names(x)
  print(out)
  return(out)
}

print.niche.models = function(x, ...) {
  model.group = deparse(substitute(x))
  best = .getBestModel(x)
  bfml = .fmla2txt(x$formulas[[best]])
  cat(model.group, "\n\tBest Model:", best, ":\n", bfml)
  return(invisible())
}


setPredictionFiles = function(prefix, start, end, aux=NULL, dir=".") {
  times = createTimeAxis(start=start, end=end, frequency=12, center=TRUE)
  files = paste0("Y",times$year, "M", sprintf("%02d", times$month))
  files = paste(prefix, files, sep="-")
  files = paste0(files, ".csv")
  files = file.path(dir, files)
  output = list(files=files, aux=aux, time=times)
  files.check = file.exists(files)
  aux.check = if(!is.null(aux)) file.exists(aux) else TRUE
  if(!all(files.check)) warning("Some prediction files don't exist")
  if(!all(aux.check)) warning("Auxiliar prediction files don't exist")
  return(output)
}

.getModelVars2 = function(formula, data) {
  ale  = as.character(formula)[3]
  vars = names(data)
  vars  = vars[unlist(lapply(vars, grepl, x=ale))]
  return(vars)  
}

.getModelVars = function(model) {
  ale  = as.character(model$formula)[3]
  vars = names(model$data)
  vars  = vars[unlist(lapply(vars, grepl, x=ale))]
  return(vars)  
}

predict.niche.models = function(object, model=NULL, pred.info, mask=NULL, radius=1,
                                precision=3, interp=!is.null(mask), FUN=identity,...) {

  
  if(!is.null(object$transform)) FUN0 = match.fun(object$transform)
  FUN  = match.fun(FUN)
  
  factor = 10^max(trunc(precision), 3)
  if( !is.null(model) && !is.na(model) ) {
    model.name = model
    model = getModel(object, model)
  } else {
    model.name = .getBestModel(object)
    model = getBestModel(object)
  }
  
  
  if(!is.null(pred.info$aux)) aux = read.csv(pred.info$aux)
  
  for(i in seq_along(pred.info$files)) {
    
    DateStamp("Prediction on ", pred.info$files[i], "...", sep="")
    
    newdata = read.csv(pred.info$files[i])
    
    if(i==1) {
      
      nlon = length(unique(newdata$lon))
      nlat = length(unique(newdata$lat))
      n = nlon*nlat
      
      if(nrow(newdata)!=n) stop("Predictions cannot be plotted in a grid") 
      pred = array(dim=c(nlon, nlat, length(pred.info$files)))
      
      lon = sort(unique(newdata$lon))
      lat = sort(unique(newdata$lat))
      
      coords = list(lon=lon, lat=lat)
    }
    
    if(!is.null(pred.info$aux)) {
      
#       newdata$lat = round(newdata$lat, 2)
#       newdata$lon = round(newdata$lon, 2)
#       aux$lat = round(aux$lat, 2)
#       aux$lon = round(aux$lon, 2)
      
      newdata = merge(newdata, aux, by=c("lon", "lat"), all=TRUE, sort=FALSE)
      
    }
    
    if(!is.null(object$transform)) newdata = FUN0(newdata)
    newdata = FUN(newdata)
    
    model.vars = .getModelVars(model) 
    ind = complete.cases(newdata[, model.vars])
    newdata = newdata[ind, ]
    
    if(nrow(newdata)==0) next
    
    temp = predict(model, type="response", newdata=newdata)
    temp = as.integer(round(factor*temp,3))
    pred[,,i][ind] = temp 
    
  }
  
  if(!identical(dim(pred)[1:2], dim(mask)) & interp) {
    warning("'map' and 'mask' dimensions don't agree. Cannot perform interpolation")
    interp = FALSE
    mask   = NULL
  }
  
  if(interp) {
    DateStamp("Performing bilinear spatial interpolation over mask...")
    pred = fillMap(pred, mask=mask, radius=radius)
    DateStamp("Performing spline temporal interpolation...")
    pred = interpolateMap(pred, anomalies=FALSE)
  }
  
  thr = object$threshold$validation[,model.name]
  fit = object$fit[model.name,]
  per = rbind(training=object$performance$training[model.name,],
              validation=object$performance$validation[model.name,])
  fml = object$formulas[model.name]
  
  info = list(coords=coords, time=pred.info$time, factor=factor,
              threshold=thr, fit=fit, performance=per, formula=fml,
              model=model.name, transform=FUN, mask=mask)
  
  DateStamp("Computing climatologies...")
  pred.mean   = apply(pred, 1:2, median, na.rm=TRUE)
  pred.sd     = apply(pred, 1:2, sd, na.rm=TRUE)
  pred.clim   = climatology(pred, info$time$month)
  pred.season = climatology(pred, info$time$season)
  
  output = list(prediction=pred, info=info, mean=pred.mean, sd=pred.sd,
                climatology=pred.clim, season=pred.season)
  
  DateStamp("DONE.")
  
  class(output) = "prediction.niche.models"
  return(output)
}

window.prediction.niche.models = 
  function(x, start=NULL, end=NULL, frequency=NULL, 
           deltat=NULL, extend=FALSE, ...) {
  
    if(is.null(start) & is.null(end)) return(x)
    
    year  = x$info$time$year
    month = x$info$time$month
    
    if(is.null(start)) start = c(year[1], month[1])
    if(is.null(end))   end   = c(tail(year, 1), tail(month, 1))
    
    oldtime = x$info$time$center
    newtime = createTimeAxis(start=start, end=end, frequency=12, center=TRUE)  
    date.min = min(newtime$center)
    date.max = max(newtime$center)
    
    time.ind = oldtime >= date.min & oldtime <= date.max 
    
    x$prediction = x$prediction[, , time.ind]
    x$info$time = newtime
    
    DateStamp("Computing climatologies for new time window...")
    x$mean        = apply(x$prediction, 1:2, median, na.rm=TRUE)
    x$sd          = apply(x$prediction, 1:2, sd, na.rm=TRUE)
    x$climatology = climatology(x$prediction, x$info$time$month)
    x$season      = climatology(x$prediction, x$info$time$season)
    

    
    return(x)
    
}



getModel = function(object, model.name) {
  
  if(!(model.name %in% names(object$models))) stop("Model ", sQuote(model.name), " not found")  
  model = object$models[[model.name]]
  dataset = deparse(substitute(object))
  model$call$data    = eval(parse(text=paste0("get(\"", dataset,"\")$train")))
  model$call$formula = as.formula(.fmla2txt(model$formula))
  
  return(model)
}

# how to reference to data inside the object

getFormula = function(object, model.name) {
  
  if(!(model.name %in% names(object$models))) stop("Model ", sQuote(model.name), " not found")  
  fmla = object$formulas[[model.name]]
  
  return(fmla)
}

.getBestModel =function(object, criteria="AUC") {
  n = which.max(object$performance$validation[,criteria])
  model = rownames(object$performance$validation)[n]
  return(model)
}

getBestModel = function(object, criteria="AUC") {
  model = getModel(object=object, model=.getBestModel(object=object, 
                                                      criteria=criteria))
  return(model)
}

getThreshold = function(x, criteria="MinROCdist") {
  thr = x$info$threshold[criteria]
  return(thr)
}

getMap = function(object, ...) {
  UseMethod("getMap")
}

getMap.prediction.niche.models = function(object, date, toPA=FALSE, prob=FALSE, 
                                          criteria="MinROCdist", ...) {
  thr = getThreshold(x=object, criteria=criteria)
  dates = data.frame(year=object$info$time$year, month=object$info$time$month)
  if(!is.null(date)) {
    if(date[2]>12 | date[2] < 1 | date[2]%%1!=0) {
      stop("month has to be an integer between 1 and 12")
    }
    slice = which(dates$year==date[1] & dates$month==date[2])
    if(length(slice)==0) stop("date is out of range")
    z = object$prediction[,,slice]/object$info$factor
    
    if(toPA) z = toPA(z, thr, prob=prob)
    output = list(lon=object$info$coords$lon, lat=object$info$coords$lat, z=z)
    return(output)
  }
}
  
plot.prediction.niche.models = function(x, y=NULL, date=NULL, 
                                        slice=NULL, type=NULL, mfclim=c(3,4),
                                        toPA=FALSE, prob=FALSE, criteria="MinROCdist", ...) {
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  thr = getThreshold(x=x, criteria=criteria)
  dates = data.frame(year=x$info$time$year, month=x$info$time$month)
  if(!is.null(date)) {
    if(date[2]>12 | date[2] < 1 | date[2]%%1!=0) {
      stop("month has to be an integer between 1 and 12")
    }
    slice = which(dates$year==date[1] & dates$month==date[2])
    if(length(slice)==0) stop("date is out of range")
    time.lab = paste(month.abb[date[2]],date[1], sep=" - ")
  } else {
    if(is.null(slice)) {
      time.lab = paste0("average ", min(dates$year), "-", max(dates$year))
    } else {
      time.lab = paste(month.abb[dates$month[slice]],
                       dates$year[slice], sep=" - ")
    }
  }
  if(is.null(date) & is.null(slice)) {
    z = x$mean/x$info$factor
    if(toPA) z = toPA(z, thr, prob=prob)
  } else {
    z = x$prediction[,,slice]/x$info$factor
    if(toPA) z = toPA(z, thr, prob=prob)
  }
  time.lab = toupper(time.lab)
  if(is.null(type)) {
    try(image.map(x$info$coords$lon, x$info$coords$lat, z, zlim=c(0,1), ...),
        silent=TRUE)
    points(-180,-180)
    mtext(time.lab, 3, line=1, adj=1, cex=0.75)    
  } else {
    if(type=="climatology") {
      par(mfrow=mfclim, oma=c(1,1,1,1), mar=c(4,4,1,1))
      z = x$climatology/x$info$factor
      if(toPA) z = toPA(z, thr, prob=prob)
      for(i in 1:12) {
        try(image.map(x$info$coords$lon, x$info$coords$lat, z[,,i], zlim=c(0,1),
                      ...),
            silent=TRUE)
        points(-180,-180)
        mtext(month.name[i], 3, line=1, adj=1, cex=0.75)            
      }
    }
    if(type=="seasonal") {
      par(mfrow=c(2,2), oma=c(1,1,1,1), mar=c(4,4,1,1))
      z = x$season/x$info$factor
      if(toPA) z = toPA(z, thr, prob=prob)
      slab = toupper(dimnames(x$season)[[3]])
      for(i in 1:4) {
        try(image.map(x$info$coords$lon, x$info$coords$lat, z[,,i], zlim=c(0,1), 
                      ...),
            silent=TRUE)
        points(-180,-180)
        mtext(slab[i], 3, line=1, adj=1, cex=0.75)            
      }
    }
    if(!(type %in% c("seasonal", "climatology"))) {
      stop("Invalid 'type' argument, it can be 'seasonal' or 'climatology'.")
    }
  }
  return(invisible())
}

saveAnimation = function(object, ...) {
  UseMethod("saveAnimation")
}

saveAnimation.prediction.niche.models = 
  function(object, file, dir=getwd(), interval=0.5, 
           ...) {
    
    n = dim(object$prediction)[3]
    DateStamp("\nCreating animation (",n," time steps).", sep="")
    
    try(saveGIF(
{
  for(i in seq_len(n)) {
    cat(i,"")
    plot(object, slice=i, ...)
  }
},
movie.name="temp.gif",
img.name="niche_model", clean=TRUE, verbose=FALSE, 
interval=interval, loop=1, check=TRUE, autobrowse=FALSE),
        silent=TRUE)
    tmp = file.path(ani.options("outdir"), "temp.gif")
    x = file.copy(from=tmp, file.path(dir, file), overwrite=TRUE)
    DateStamp("DONE.")
    return(invisible(x))
  }

saveAnimation.default = 
  function(object, file, dir=getwd(), interval=0.5, 
           ...) {
    
    n = dim(object)[3]
    DateStamp("\nCreating animation (",n," time steps).", sep="")
    
    try(saveGIF(
{
  for(i in seq_len(n)) {
    cat(i,"")
    image.plot(object[,,i], ...)
  }
},
movie.name="temp.gif",
img.name="slice", clean=TRUE, verbose=FALSE, 
interval=interval, loop=1, check=TRUE, autobrowse=FALSE),
        silent=TRUE)
    tmp = file.path(ani.options("outdir"), "temp.gif")
    x = file.copy(from=tmp, file.path(dir, file), overwrite=TRUE)
    DateStamp("DONE.")
    return(invisible(x))
  }


climatology = function(object, index, FUN="median") {
  FUN = match.fun(FUN)
  dims = dim(object)
  ind  = unique(index)
  uind = length(ind)
  if(is.numeric(ind)) ind = sort(ind)
  out = array(0L, dim=c(dims[1:2], uind))
  for(i in seq_len(uind)) {
    out[,,i] = as.integer(apply(object[,, which(index == ind[i])], 1:2, FUN, na.rm=TRUE))
  }
  dimnames(out)[[3]] = ind
  return(out)
}


anova.niche.models = function(object, model=NULL, criteria="AUC", ...) {
  if(is.null(model)) model=.getBestModel(object=object, criteria=criteria)
  print(object$models[[model]]$anova)
  return(invisible(NULL))
}

printMaps = function(object, dir, prefix, width=800, height=800, ...) {
  png(filename=file.path(dir, paste0(prefix,"-avg.png")),
      width=width, height=height)
  plot(object, ...)
  dev.off()
  png(filename=file.path(dir, paste0(prefix,"-season.png")),
      width=width, height=height)
  plot(object, type="seasonal", ...)
  dev.off()
  png(filename=file.path(dir, paste0(prefix,"-clim.png")), 
      width=2*width, height=1.5*height)
  plot(object, type="climatology", ...)
  dev.off()
  return(invisible())
}

plot.niche.models = function(x, model=NULL, 
                             mar=c(5,5,1,1), oma=c(0.5,1,4,1.5), ...) {
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  model.name = if(is.null(model)) .getBestModel(x) else model 
  dataset = deparse(substitute(x))
  if(model.name!="all") {
    
    if(!(model.name %in% names(x$models))) stop("Model ", sQuote(model.name), " not found")
    .plot.niche.models(x=x, model.name=model.name, dataset=dataset, mar=mar, oma=oma, ...)
    
  } else {
    
    models = names(x$models)
    
    for(model.name in models) {
    .plot.niche.models(x=x, model.name=model.name, dataset=dataset, mar=mar, oma=oma, ...)  
    }
    
  }
  
  return(invisible())
  
}

.plot.niche.models = function(x, model.name, dataset, 
                              mar=c(5,5,1,1), oma=c(0.5,1,4,1.5), ...) {
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  
  model = x$models[[model.name]]
  fml = .fmla2txt(model$formula)
  model$call$data    = eval(parse(text=paste0("get(\"", dataset,"\")$train")))
  model$call$formula = as.formula(fml)
  
  nvars = length(labels(model$terms))
  var = as.character(model$formula)[2]
  par(mfrow=.getmfrow(nvars), mar=mar, oma=oma)
  
  plot(model, se=TRUE, ...)
  mtext(paste0(var, " (", dataset, ")"), 3, line=2, outer=TRUE)
  mtext(fml, 3, line=1, outer=TRUE)
  
}

getCoordinates = function(object, index) {
  if(is.logical(index)) index=which(index)
  if(!is.numeric(index)) stop("index must be logical or numerical")
  cols = col(object)[index]
  rows = row(object)[index]
  out = list(x=rows, y=cols)
  return(out)
}


map2coord = function(object) {
  out = data.frame(x = as.numeric(row(object)),
                   y = as.numeric(col(object)),
                   z = as.numeric(object)
  )
  out = out[complete.cases(out),]
}

.getIndex = function(x, max, radius) {
  x = trunc(x)
  ind = seq(from=x-radius, to=x+radius, by=1)
  ind = ind[ind>0 & ind<=max]
  return(ind)
}



fillSquare = function(x, map, radius=1) {
  x = x[!complete.cases(x),]
  .fillSquare = function(x, map, radius) {
    ix = .getIndex(x[1], max=nrow(map), radius=radius)
    iy = .getIndex(x[2], max=ncol(map), radius=radius)
    z  = median(map[ix, iy], na.rm=TRUE)
    return(z)
  }
  out = apply(x, 1, .fillSquare, map=map, radius=radius)
  return(out)
}

fillMap = function(object, mask, radius=1, ...) {
  UseMethod("fillMap")
}


fillMap.matrix = function(object, mask, radius=1) {
  
  if(!all(is.na(object))) {
    if(!identical(dim(object), dim(mask))) stop("'object' and 'mask' dimension must agree.")
    coords = map2coord(object)
    ind    = is.na(object)&!is.na(mask)
    miss   = getCoordinates(object, ind)
    fill   = interpp(coords$x, coords$y, coords$z, miss$x, miss$y)
    map2   = object
    map2[ind] = fill$z
    fill   = as.data.frame(fill)
    fill$z[is.na(fill$z)] = fillSquare(x=fill, map=map2, radius=radius)
    fill$z[is.na(fill$z)] = 0 # just in case!
    map2[ind] = fill$z
    map2[is.na(mask)] = NA
  } else {
    map2 = object
  }
  return(map2)
}

fillMap.array = function(object, mask, radius=1) {
  dmap = dim(object)
  out = apply(object, 3, fillMap, mask=mask, radius=radius)
  dim(out) = dmap
  return(out)
}


fillMap.prediction.niche.models = function(object, mask, radius=1) {
  if(!identical(dim(object$prediction)[1:2], dim(mask))) 
    stop("'object' and 'mask' dimension must agree.")
  
  DateStamp("Performing spatial interpolation over mask...")
  object$prediction = fillMap(object$prediction, mask=mask, radius=radius)
  object$info$mask = mask
  
  DateStamp("Performing temporal interpolation...")
  object$prediction = interpolateMap(object$prediction, anomalies=FALSE)
  
  DateStamp("Recalculating climatologies...")
  object$mean         = apply(object$prediction, 1:2, median, na.rm=TRUE)
  object$sd           = apply(object$prediction, 1:2, sd, na.rm=TRUE)
  object$climatology  = climatology(object$prediction, object$info$time$month)
  object$season       = climatology(object$prediction, object$info$time$season)
  
  DateStamp("DONE.")
  return(object)
  
}

.interpolateMap = function(x, anomalies=FALSE, ...) {
  
  if(all(is.na(x)) | all(!is.na(x))) return(x)
  
  out = x
  t = seq_along(x)
  
  if(anomalies) {
    xclim = tapply(x, INDEX=rep(1:12, len=length(x)), FUN=mean, na.rm=TRUE)
    xclim = rep(xclim, len=length(out))
    out = x - xclim
  }
  
  xs = splinefun(x=t, y=out)
  out[is.na(out)] = pmin(pmax(xs(t[is.na(out)]),min(out, na.rm=TRUE)),
                         max(out, na.rm=TRUE))
  
  if(anomalies) out = out + xclim
  
  return(out)
}

interpolateMap = function(object, anomalies=FALSE, ...) {
  pred = apply(object, 1:2, .interpolateMap, anomalies=anomalies, ...)
  pred = aperm(pred, c(2,3,1))
  return(pred)
}


addLHT = function(object, sp, lifespan=1, 
                  ages=seq_len(ceiling(lifespan))-1, ...) {
  pars = list(sp=sp, lifespan=lifespan, ages=ages, ...)
  object$info$LHT = pars
  return(invisible(object))
}


normalize = function(x) {
  
  if(length(dim(x))> 3) stop("'x' cannot be an array of dimension greater than 3.")
  
  if(length(dim(x))< 3) out = x/sum(x,na.rm=TRUE) 
  if(length(dim(x))==3) out = apply(x, 3, normalize)
  dim(out) = dim(x)
  return(out)
  
}



.inArea = function(x, lat, lon) {
  cnull = c(-Inf, Inf)
  if(is.null(lat)) lat = cnull
  if(is.null(lon)) lon = cnull
  ilat = x$lat >= lat[1] & x$lat <= lat[2]
  ilon = x$lon >= lon[1] & x$lon <= lon[2]
  ind = ilat & ilon
  return(ind)
}

.ponderate = function(w, x) {
  ind = complete.cases(cbind(as.numeric(w), as.numeric(x)))
  w = w[ind]
  x = x[ind]
  out = weighted.mean(x=x, w=w, na.rm=TRUE)
  return(out)
}

ponderate = function(x, aux, vars=names(aux), lat=NULL, lon=NULL) {
  
  DateStamp("Starting calculations...")
  vv = vars %in% names(aux)
  if(sum(vv)==0) stop("No valid variables 'vars' in auxiliar file.")
  vars = vars[vv]
  
  dw = dim(x)
  if(length(dw)<2) stop("'w' must be an array or matrix")
  ncol = if(length(dw)==2) 1 else dw[3]
  w = matrix(normalize(x), ncol=ncol)
  
  if(!is.null(lat) | !is.null(lon)) {
    if(!is.null(lat)) lat = sort(lat) 
    if(!is.null(lon)) lon = sort(lon) 
    ind = .inArea(aux, lat=lat, lon=lon)
    aux = aux[ind,]
    w   = w[ind,]    
  }
  
  output = NULL
  for(var in vars) {
    DateStamp("Computing statistics for ", sQuote(var), "...", sep="")
    output$center[[var]] = apply(w*aux[,var], 2, .na.sum)
    output$density[[var]] = getDensity(x=aux[,var], w=w)
    
  }
  
  output$center = as.data.frame(output$center)
  DateStamp("DONE.")
  return(output)
}

.na.sum = function(x) if(!all(is.na(x))) sum(x, na.rm=TRUE) else NA

.getDensity = function(w, x, n=256, adjust=1.5) {
  if(all(is.na(w)) || all(is.na(x))) return(rep(NA, length=n))
  ran = range(pretty(x), na.rm=TRUE)
  ind = complete.cases(data.frame(w,x))
  x = x[ind]
  w = normalize(w[ind])
  out = density(x=x, weights=w, n=n, adjust=adjust,
                from=ran[1], to=ran[2])
  out = y=out$y
  return(out)
}

getDensity = function(x, w, n=256, adjust=1.5) {
  out = apply(w, 2, .getDensity, x=x, n=n, adjust=adjust)
  ale = range(pretty(x))
  ale = seq(from=ale[1], to=ale[2],len=n)
  out = list(xaxis=ale, densities=out)
  return(out)
}


getRleIndex = function(x) {
  x = rle(x)
  nslice = length(x$length)
  nrow   = cumsum(x$lengths)
  nrow   = cbind(begin=c(1, nrow[-nslice] + 1), end=nrow)
  return(nrow)
}

getValuesinRange = function(index, x, FUN, ...) {
  FUN = match.fun(FUN)
  .FUN = function(index, x, ...) FUN(x[index[1]:index[2]], ...)
  output = apply(index, 1, .FUN, x=x, ...)
  return(output)
}

getYearMin = function(object, index) {
  ind = getRleIndex(object$info$time[[index]])
  out = getValuesinRange(index=ind, x=object$info$time$year, FUN=min)
  out = out - min(out, na.rm=TRUE)
  return(out)
}

getYearMax = function(object, index) {
  ind = getRleIndex(object$info$time[[index]])
  out = getValuesinRange(index=ind, x=object$info$time$year, FUN=max)
  out = out - min(out, na.rm=TRUE) + 1
  return(out)
}


getSesonalMap = function(object) {
  nrow = getRleIndex(object$info$time$season)
  map  = object$prediction
  nslice = nrow(nrow)
  output = array(dim=c(dim(map)[1:2], nslice))
  
  for(i in seq_len(nslice)) {
    ali = nrow[i,1]:nrow[i,2]
    ale = map[,, ali]
    output[,,i] = apply(ale, 1:2, mean, na.rm=TRUE)*object$info$mask
  }
  return(output)
}

year2month = function(data, year=seq_along(data), month) {
  
  
  var.name = rev(strsplit(deparse(substitute(data)), "\\$")[[1]])[1]
  years = rep(year, each=12)
  months = rep(1:12, len=length(years))
  new.data = rep(NA, len=length(years))
  new.data[months==month] = data
  output = data.frame(year=years, month=months, new.data)
  names(output)[3] = var.name
  return(output)
  
}

.removeSector3 = function(x, coords, lon, lat) {
  
  xlat = is.inside(coords$lat, lat)
  xlon = is.inside(coords$lon, lon)
  xpred = x
  
  if(!is.null(lat) & !is.null(lon)) {
    xpred[xlon, xlat, ] = 0
  } else {
    if(is.null(lon)&!is.null(lat)) {
      xpred[,xlat,] = 0
    }
    if(!is.null(lon)&is.null(lat)) {
      xpred[xlon,,] = 0      
    }
  }
  xpred[is.na(x)] = NA
  return(xpred)
}

.removeSector2 = function(x, coords, lon, lat) {
  
  xlat = is.inside(coords$lat, lat)
  xlon = is.inside(coords$lon, lon)
  xpred = x
  
  if(!is.null(lat) & !is.null(lon)) {
    xpred[xlon, xlat] = 0
  } else {
    if(is.null(lon)&!is.null(lat)) {
      xpred[,xlat] = 0
    }
    if(!is.null(lon)&is.null(lat)) {
      xpred[xlon,] = 0      
    }
  }
  xpred[is.na(x)] = NA
  return(xpred)
}

removeSector = function(x, coords, lon, lat) {
  ndim = length(dim(x))
  if(ndim==2) out = .removeSector2(x, coords, lon, lat)
  if(ndim==3) out = .removeSector3(x, coords, lon, lat)
  return(out)
}

