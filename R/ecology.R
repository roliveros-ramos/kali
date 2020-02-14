
assign_areas = function(data, what, latitude="lat", longitude="lon") {
  
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

