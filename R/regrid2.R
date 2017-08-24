

regrid2 = function(object, old, new, mask, linear=TRUE, extrap=FALSE, ...) {
  UseMethod("regrid2")
}

regrid2.matrix = function(object, old, new, mask=NULL, linear, extrap, ...) {
  
  if(is.null(mask)&!is.null(new$mask)) mask=new$mask
  stopifnot(exists("lat", where=old), exists("lon", where=old))
  stopifnot(length(old$lat)==ncol(object), length(old$lon)==nrow(object))
  
  stopifnot(exists("lat", where=new), exists("lon", where=new)) 
  stopifnot(is.numeric(new$lat), !is.matrix(new$lat), 
            is.numeric(new$lon), !is.matrix(new$lon))
  
  old$x = rep(old$lon, ncol(object))
  old$y = rep(old$lat, each=nrow(object))
  old$z = as.numeric(object)

  old = old[c("x","y","z")]
  
  old = as.data.frame(old)
  old = old[complete.cases(old), ]
  
  new$x = new$lon
  new$y = new$lat
  
  newp = interp(x=old$x, y=old$y, z=old$z, xo=new$x, yo=new$y, 
                linear=linear, extrap=extrap, ...)$z
  newmap = if(!is.null(mask)) newp*mask else newp
  
  return(newmap)
}


regrid2.array = function(object, old, new, mask=NULL, linear, extrap, ...) {
  # new grid
  if(exists("LAT", where=new) & exists("LON", where=new)) {
    stopifnot(is.matrix(new$LAT), is.matrix(new$LON), dim(new$LAT)==dim(new$LON))
    nLAT = new$LAT
    nLON = new$LON
    nlat = as.numeric(nLAT)
    nlon = as.numeric(nLON)
  } else {
    if(exists("lat", where=new) & exists("lon", where=new)) {
      stopifnot(is.numeric(lat), !is.matrix(lat), is.numeric(lon), !is.matrix(lon))
      nLAT = matrix(new$lat, ncol=length(new$lat), nrow=length(new$lon), byrow=TRUE)
      nLON = matrix(new$lon, ncol=length(new$lat), nrow=length(new$lon))
      nlat = as.numeric(nLAT)
      nlon = as.numeric(nLON)
    } else {
      stop("'new' must contain latitude and longitude information.")
    }
  }
  
  ndim = seq_along(dim(object))[-c(1,2)]
  newmap = apply(object, ndim, regrid2, old=old, new=new, mask=mask, 
                 linear=linear, extrap=extrap, ...)
  dim(newmap) = c(dim(nLAT), dim(object)[-c(1,2)])
  return(newmap)
}
