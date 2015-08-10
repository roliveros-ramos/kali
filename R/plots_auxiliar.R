# create grid axes (lat, lon)
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

# create time axis
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

# rename map.axes an overwrite map's map.axes?
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


# this function change log coordinates from -180,+180 to 0-360
# change name
rotateMap =  function(x) {
  out = x$lon>0
  out[is.na(out)] = FALSE
  x$lon[out] = x$lon[out] - 360
  return(x)
}


# change the center of a map database
map2 =  function(database, center, ...){
  # change the center of the map (e.g. to show complete Pacific Ocean)
  # From stackoverflow (check reference and how to cite)
  Obj = map(database, ..., plot=FALSE)
  coord = cbind(Obj[[1]],Obj[[2]])
  # split up the coordinates
  id = rle(!is.na(coord[,1]))
  id = matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=TRUE)
  polygons = apply(id, 1, function(i) {coord[i[1]:i[2],]} )
  # split up polygons that differ too much
  polygons = lapply(polygons, function(x) {
    x[,1] = x[,1] + center
    x[,1] = ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id = x[,1] < 0
      x = rbind(x[id,],c(NA,NA),x[!id,])
    }
    return(x)
  })
  
  # reconstruct the object
  polygons = do.call(rbind,polygons)
  Obj[[1]] = polygons[,1]
  Obj[[2]] = polygons[,2]
  
  map(Obj,...)
}
 
# create method fill??
fill.map = function(x, mask=1, FUN="max") {
  
  # fill a map with the max (min, mean, median) value of the map
  
  FUN = match.fun(FUN)
  maximo = FUN(x, na.rm=TRUE)
  x[is.na(x)] = maximo
  x = x*mask
  
  return(x)
}

# create method convexHull, this will be convexHull.array
convexHull = function(x, mask) {
  
  dim.old = dim(x)
  out = apply(x, 3, convexHull.map, mask=mask)  
  dim(out) = dim.old
  
  return(out)
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

# change resolution of a map
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

