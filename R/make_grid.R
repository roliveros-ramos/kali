
#' Make a rectangular grid
#'
#' @param lon Longitude of the left and right corners of the grid. 
#' @param lat Latitude of the upper and lower corners of the grid. 
#' @param dx Resolution in the longitude axis (in degrees).
#' @param dy Resolution in the latitude axis (in degrees).
#' @param n Number of points within a square to test ocean/land mask. See details.
#' @param hires Use high resolution coast line? Default to FALSE, TRUE will take considerly longer time to compute. 
#'
#' @return A list, including the grid and land/ocean mask.
#' @export
#'
#' @examples 
#' grid = make_grid(lon=c(2, 9), lat=c(40,45), dx=1/12, dy=1/12, n=1)
#' plot(grid)
make_grid = function(lon=lon, lat=lat, dx=dx, dy=dy, n=1, thr=0.8, hires=FALSE, mask=TRUE) {
  
  lon = range(lon, na.rm=TRUE)
  lat = range(lat, na.rm=TRUE)
  
  grid = .create_grid(lon=lon, lat=lat, dx=dx, dy=dy)
  if(isTRUE(mask)) {
    mm = .create_mask(lon=lon, lat=lat, dx=dx, dy=dy, n=n, thr=thr, hires=hires)
    grid$mask  = mm$mask
    grid$prob  = mm$ocean
    grid$n     = mm$n
    grid$hires = hires
  } else {
    grid$mask  = NULL
    grid$prob  = NULL
    grid$n     = 0
    grid$hires = FALSE
  }
  
  class(grid) = c("grid", class(grid))
  return(grid)
  
}

#' Update grid
#'
#' @param grid An object of class 'grid' 
#' @param thr A value between 0 and 1, used as threshold to define an ocean cell.
#'
#' @return A grid object.
#' @export
#'
update_grid = function(grid, thr) {
  if(!inherits(grid, "grid")) stop("This is not a valid grid object.")
  grid$mask = 0 + (grid$prob >= thr)
  return(grid)
}
# Methods -----------------------------------------------------------------

#' @exportS3Method plot grid
plot.grid = function(x, land.col="darkolivegreen4", sea.col="aliceblue", prob=FALSE,
                     boundaries.col="black", grid=TRUE, grid.lwd=1, grid.col="lightgray", lwd=2, ...) {
  
  if(is.null(x$mask)) {
    z = x$LAT 
    col = sea.col
    } else {
    if(isTRUE(prob)) {
      z = x$prob 
      col = colorful::divergencePalette(n=x$n, zlim=c(0,1), col=c(land.col, sea.col), center = 0.5, p=0.4)
    } else {
      z = x$mask
      col = c(land.col, sea.col)
      } 
  }  
    
    image.map(x$psi$lon, x$psi$lat, z, land=FALSE, legend=prob, zlim=c(0,1),
              col=col, grid=grid, grid.lwd=grid.lwd, grid.col=grid.col)
    map_details(fill=FALSE, col=boundaries.col, lwd=lwd, hires=x$hires, ...)
    return(invisible(NULL))
  
}

# Auxiliar functions ------------------------------------------------------

#' Are the points in the ocean (or land)?
#'
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes. 
#' @param hires Boolean. Use high resolution mask to compute ocean/land?
#'
#' @return A logical vector. For is_ocean, TRUE when the point is on water, FALSE for land or NA. Similiar for is_land.
#' @export
#'
is_ocean = function(lon, lat, hires=FALSE) {
  
  coords = cbind(lon=lon, lat=lat)
  xind = which(complete.cases(coords))
  coords = coords[xind, ]
  
  out = rep(FALSE, length=nrow(coords))
  
  layer = if(isTRUE(hires)) land_hr else land_lr

  crs = suppressWarnings(CRS(proj4string(layer)))
  sets = SpatialPoints(cbind(lon=coords[, "lon"], lat=coords[, "lat"]), 
                       proj4string=crs)
  ind = xind[which(is.na(over(sets, layer)))]
  out[ind] = TRUE
  return(out)
}

#' @rdname is_ocean
#' @export
is_land = function(lon, lat, hires=FALSE) return(!is_ocean(lon, lat, hires))
  
# Internal functions ------------------------------------------------------



.create_grid = function(lon, lat, dx, dy, center=FALSE) {
  
  # Create a rectangular grid given lat, lon and dxy.
  # No correction by Earth curvature
  if(dx <= 0 || dy <= 0) stop("dx and dy must be positive.")
  if(diff(lon) <= 0) stop("Longitude values must not be decreasing.")
  if(diff(lat) <= 0) stop("Latitude values must not be decreasing.")
   
  nx = floor(diff(lon)/dx) 
  ny = floor(diff(lat)/dy)
  
  dlon = dx*nx - diff(lon)
  dlat = dy*ny - diff(lat)
  
  if(abs(dlon) > 1e-6*dx) {
    nx  = nx + 1
    dlon = dx*nx - diff(lon)
    lon = lon + 0.5*dlon*c(-1, +1)
  }  
  if(abs(dlat) > 1e-6*dy) {
    ny  = ny + 1
    dlat = dy*ny - diff(lat)
    lat = lat + 0.5*dlat*c(-1, +1)
  }  
  
  if(isTRUE(center)) {
    lat[which.min(lat)] = lat[which.min(lat)] + 0.5*dy
    lat[which.max(lat)] = lat[which.max(lat)] - 0.5*dy
    lon[which.min(lon)] = lon[which.min(lon)] + 0.5*dx
    lon[which.max(lon)] = lon[which.max(lon)] - 0.5*dx
    nx = nx + 1
    ny = ny + 1
  }
  
  lats.psi = seq(from=min(lat),to=max(lat), length=ny + 1)
  lons.psi = seq(from=min(lon),to=max(lon), length=nx + 1)
  
  lats.rho = seq(from=min(lat) + 0.5*dy, to=max(lat) - 0.5*dy, length=ny)
  lons.rho = seq(from=min(lon) + 0.5*dx, to=max(lon) - 0.5*dx, length=nx)
  
  rho = list(lat=lats.rho, lon=lons.rho)
  psi = list(lat=lats.psi, lon=lons.psi)
  
  nlat = length(rho$lat)
  nlon = length(rho$lon)
  
  LAT = matrix(rho$lat, ncol=nlat, nrow=nlon, byrow=TRUE)
  LON = matrix(rho$lon, ncol=nlat, nrow=nlon)
  
  area = (111*dy)*(111*cos(LAT*pi/180)*dx)
  
  output = list(lon=lons.rho, lat=lats.rho, rho=rho, psi=psi, 
                LON=LON, LAT=LAT, area=area, 
                df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))
  
  return(output)
  
}

.grid_offset = function(n=0, dx, dy) {
  if(n==0) return(NULL)
  n = 2*n + 1
  x = seq(from=0, to=dx, length=n) - 0.5*dx
  y = seq(from=0, to=dy, length=n) - 0.5*dy
  return(expand.grid(y=y, x=x)[, 2:1])
}

.create_mask = function(lon, lat, dx, dy, n=2, thr=0.8, hires=FALSE) {
  
  grid = .create_grid(lon=lon, lat=lat, dx=dx, dy=dy)
  off  = .grid_offset(n=n, dx=dx, dy=dy)
  
  if(is.null(off)) {
    out = grid$df
  } else {
    out = list()
    for(i in seq_len(nrow(off))) {
      tmp = grid$df
      tmp$lon = tmp$lon + off$x[i] 
      tmp$lat = tmp$lat + off$y[i] 
      out[[i]] = tmp 
    }
    out = do.call(rbind, out)
  }
  
  ind = is_ocean(lon=out$lon, lat=out$lat, hires=hires)
  
  gg = array(0 + ind, dim=c(dim(grid$area), nrow(off)))
  gg = apply(gg, 1:2, mean)
  
  mask = 0 + (gg >= thr)
  
  output = list(mask=mask, ocean=gg, n=nrow(off)+1)
  
  return(output)
} 

