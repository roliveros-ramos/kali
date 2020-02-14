
#' Assign areas to geographical coordinates.
#'
#' @param data A data frame including spatial coordinates
#' @param what the database with the geographic areas. 
#' Currently, only "longhurst" is accepted.
#' @param longitude The name of the variable in \code{data} with longitud
#' information, by default "lat"
#' @param latitude The name of the variable in \code{data} with latitude
#' information, by default "lat"
#' @return
#' @export
#'
#' @examples
assign_areas = function(data, what, longitude="lon", latitude="lat") {
  
  layer = switch(what,
                 longhurst = kali::longhurst,
                 stop(sprintf(what, "database not available.")))
  
  key   = switch(what,
                 longhurst = "ProvCode",
                 stop(sprintf(what, "database not available.")))
  
  crs = proj4string(layer)
  pp = SpatialPoints(data[, c(longitude, latitude)], proj4string = CRS(crs))
  ind = over(pp, layer)[, key]
  return(ind)
  
}

