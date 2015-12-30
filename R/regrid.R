

regrid = function(object, old, new, mask, ...) {
  UseMethod("regrid")
}

regrid.matrix = function(object, old, new, mask=NULL, ...) {
  # check mask
  if(is.null(mask)&!is.null(new$mask)) mask=new$mask
  # old grid
  if(exists("LAT", where=old) & exists("LON", where=old)) {
    stopifnot(dim(old$LAT)==dim(old$LON), dim(old$LAT)==dim(object))
    ilat = as.numeric(old$LAT)
    ilon = as.numeric(old$LON)
  } else {
    if(exists("lat", where=old) & exists("lat", where=old)) {
      stopifnot(length(old$lat)==ncol(object), length(old$lon)==nrow(object))
      iLAT = matrix(old$lat, ncol=ncol(object), nrow=nrow(object), byrow=TRUE)
      iLON = matrix(old$lon, ncol=ncol(object), nrow=nrow(object))
      ilat = as.numeric(iLAT)
      ilon = as.numeric(iLON)
    } else stop("'old' must contain latitude and longitude information.")
  }  
    # new grid
    if(exists("LAT", where=new) & exists("LON", where=new)) {
      stopifnot(is.matrix(new$LAT), is.matrix(new$LON), dim(new$LAT)==dim(new$LON))
      nLAT = new$LAT
      nLON = new$LON
      nlat = as.numeric(nLAT)
      nlon = as.numeric(nLON)
    } else {
      if(exists("lat", where=new) & exists("lat", where=new)) {
        stopifnot(is.numeric(new$lat), !is.matrix(new$lat), 
                  is.numeric(new$lon), !is.matrix(new$lon))
        nLAT = matrix(new$lat, ncol=length(new$lat), nrow=length(new$lon), byrow=TRUE)
        nLON = matrix(new$lon, ncol=length(new$lat), nrow=length(new$lon))
        nlat = as.numeric(nLAT)
        nlon = as.numeric(nLON)
      } else {
        stop("'new' must contain latitude and longitude information.")
      }
    }
      if(!is.null(mask)) {
        stopifnot(dim(mask)==dim(nLAT))
      } else mask=1
      
      iold = data.frame(lat=ilat, lon=ilon, z=as.numeric(object))
      iold = iold[complete.cases(iold), ]
      
      newp   = interpp(iold$lon, iold$lat, iold$z, 
                       xo=nlon, yo=nlat, ...)$z
      
      newmap = matrix(newp, ncol=ncol(nLAT), nrow=nrow(nLAT))*mask
      
      return(newmap)
}


regrid.array = function(object, old, new, mask=NULL, ...) {
  # new grid
  if(exists("LAT", where=new) & exists("LON", where=new)) {
    stopifnot(is.matrix(new$LAT), is.matrix(new$LON), dim(new$LAT)==dim(new$LON))
    nLAT = new$LAT
    nLON = new$LON
    nlat = as.numeric(nLAT)
    nlon = as.numeric(nLON)
  } else {
    if(exists("lat", where=new) & exists("lat", where=new)) {
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
  newmap = apply(object, ndim, regrid, old=old, new=new, mask=mask, ...)
  dim(newmap) = c(dim(nLAT), dim(object)[-c(1,2)])
  return(newmap)
}


    