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


mapDetails = function(center, hires=FALSE, col="black", interior=FALSE, ...) {
  mapa = if (hires) {
    require(mapdata)
    "worldHires"
  }
  else {
    require(maps)
    "world"
  }
  map2(mapa, center = center, fill = TRUE, col = col, add = TRUE, 
       interior=interior, ...)
  map.axes2()
  mtext("LONGITUDE", 1, line = 1.8, cex = 0.9)
  mtext("LATITUDE", 2, line = 2.4, cex = 0.9)
  return(invisible())
}

image.map = function (lon, lat, z, center=0, hires=FALSE, add = FALSE, nlevel = 1000, horizontal = FALSE, 
                      legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
                                                                                   3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
                      bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = rev(rainbow(nlevel/10, 
                                                                                               start = 0/6, end = 4/6)), 
                      lab.breaks = NULL, axis.args = NULL, legend.args = NULL, 
                      midpoint = FALSE, border = NA, lwd = 1, land.color="black", ...) 
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
      mapDetails(center=center, hires=hires, col=land.color)
    }
    else {
      poly.image(x=lon, y=lat, z=z, add = add, col = col, midpoint = midpoint, 
                 border = border, lwd.poly = lwd, axes=FALSE, 
                 xlab="", ylab="",...)
      mapDetails(center=center, hires=hires, col=land.color)
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

.plotSea = function(col="lightblue", x0=-200, y0=-100) {
  x1 = -x0
  y1 = -y0
  polygon(x=c(x0,x1,x1,x0), y=c(y0,y0,y1,y1), border=NA, col=col)
  return(invisible())
}

plot.map = function(x, y=NULL, xlim=NULL, ylim=NULL, center=0, hires=FALSE, 
                    land.col="black", sea.col="lightblue", main=NULL, ...) {
  xy = xy.coords(x, y)
  xy$xlab = ""
  xy$ylab = ""
  if(is.null(xlim)) xlim = range(xy$x, na.rm=TRUE)
  if(is.null(ylim)) ylim = range(xy$y, na.rm=TRUE)
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  .plotSea(col=sea.col)
  points(xy, ...)
  mapDetails(center=center, hires=hires,col=land.col, interior=FALSE)
  title(main=main)
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

.rotate = function(x) {
  # clockwise
  if(length(dim(x))==2) {
    x = t(x)
    x = x[, ncol(x):1]
  }
  if(length(dim(x))==3) {
    x = aperm(x, c(2,1,3))
    x = x[, ncol(x):1, ]
  }  
  return(x)
}

.rotate2 = function(x) {
  # anticlockwise
  if(length(dim(x))==2) {
    x = x[, ncol(x):1]
    x = t(x)
  }
  if(length(dim(x))==3) {
    x = x[, ncol(x):1, ]    
    x = aperm(x, c(2,1,3))
  }  
  return(x)
}

rotate = function(x, direction="clockwise") {
  out = switch(direction,
               clockwise     = .rotate(x),
               anticlockwise = .rotate2(x))
  return(out)
}

.coord2text = function(coord, type) {
  # write nicely coordinates
  hemi = if(coord==0) {
    rep("º", 2)
  } else {
    if(coord>0) c("ºN","ºE") else c("ºS","ºW")
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


table2grid = function(data, var, lat, lon, dx, dy=dx, FUN=sum) {
  
  FUN = match.fun(FUN)
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  
  latAsFactor = cut(data[,"lat"], coords$psi$lat, labels=FALSE)
  lonAsFactor = cut(data[,"lon"], coords$psi$lon, labels=FALSE)
  
  map = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor),
               FUN=FUN, na.rm=TRUE)
  
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
  
  return(map)
  
}

createTimeAxis = function(start, end, frequency=12, center=FALSE) {
  times  = seq(from=start[1] + (start[2]-1)/frequency, 
               to= end[1] + end[2]/frequency, by=1/frequency)
  if(center) {
    out = list(bounds = times, center = times[-length(times)] + 0.5/frequency)
    out$year = floor(out$center)
    out$month = round(frequency*(out$center-out$year) + 0.5, 0)
    times = out
  } 
  return(times)
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

toPA= function(x, thr=0) {
  x[x>thr] = 1 
  x[x<=thr] = 0
  return(x)
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
    if(all(!is.na(xy))) {
      grados=acos(ref[,4]*xy[4] + ref[,3]*xy[3]*cos(ref[,2]-xy[2]))
      ind  = which.min(grados)
      dist = grados[ind] # busco la distancia minima en Km en este intervalo
      if(!is.null(ref2abs)) {
        dist = c(dist, ref2abs[ind])
      }
    } else {
      dist = NA
    }
    return(dist)
  }
  
  out = t(apply(xy, 1, .getDist, ref=ref, ref2abs=ref2abs))
  
  return(out)
}

getDistance = function(data, ref, ref2abs=NULL, lon="lon", lat="lat") {
  
  x    = as.matrix(data[,c(lat, lon)])
  cx   = as.matrix(ref[,c(lat, lon)])
  
  grados2km  = 1.852*60*180/pi
  grados2rad = pi/180
  
  if(!is.null(ref2abs)) ref2abs = ref2abs/grados2km
  
  ref = cx*grados2rad
  ref = cbind(ref, cos(ref[,1]), sin(ref[,1])) # c(lat, lon, cos(lat), sin(lat))
  xy  = x*grados2rad  
  xy  = cbind(xy, cos(xy[,1]), sin(xy[,1]))
  output = .getDistance(xy=xy, ref=ref, ref2abs=ref2abs)*grados2km
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

DateStamp = function(...) cat(..., date(), "\n\n")


getIsoline = function(lon, lat, z, level, cutoff=0) {
  .getLine = function(x) rbind(cbind(x$x, x$y), NA)
  out = grDevices::contourLines(x=lon, y=lat, z=z, nlevels=1, levels=level)
  ll = unlist(lapply(out, function(x) length(x$x)))
  ns = length(ll)
  cat("Creating ", ns, " isolines of lengths ", paste(ll[-ns], collapse=", "),
      " and ", ll[ns], ".\n", sep="")
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
