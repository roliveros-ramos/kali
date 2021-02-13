.longitude2Center = function(x, ...) {
  if (!any(x > 180, na.rm=TRUE)) 
    return(x)
  x[which(x > 180)] = x[which(x > 180)] - 360
  return(x)
}

.longitude2Left = function(x, ...) {
  if (!any(x < 0, na.rm=TRUE)) 
    return(x)
  x[which(x < 0)] = x[which(x < 0)] + 360
  return(x)
}

checkLongitude = function(x, primeMeridian="center", ...) {
  if(all(is.na(x))) return(x)
  out = switch(primeMeridian, 
               center = .longitude2Center(x, ...),
               left = .longitude2Left(x, ...))
  return(out)
}

.findPrimeMeridian = function(x) {
  if(all(is.na(x))) return("center")
  if(any(x<0)) return("center")
  if(any(x>180)) return("left")
  return("center")
}


validPoints = function(lon, lat, land=FALSE) {
  
  coords = cbind(lon=lon, lat=lat)
  xind = which(complete.cases(coords))
  coords = coords[xind, ]
  
  land = map(database="worldHires", fill=TRUE)
  land = map2SpatialPolygons(land, 
                             IDs=sapply(strsplit(land$names, ":"), FUN="[", i=1), 
                             proj4string=CRS("+proj=longlat"))
  sets = SpatialPoints(cbind(lon=coords[, "lon"], lat=coords[, "lat"]), 
                       proj4string=CRS("+proj=longlat"))
  ind = is.na(over(sets, land))
  
  if(isTRUE(land)) ind = !ind
  ind = xind[which(ind)]
  
  return(ind)
}
