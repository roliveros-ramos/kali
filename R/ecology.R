
#' Assign region to geographical coordinates.
#'
#' @param data A data frame including spatial coordinates
#' @param what the database with the geographic areas. 
#' Currently, only "longhurst", "LME", "FAO" and "GFMC" are accepted.
#' @param longitude The name of the variable in \code{data} with longitud
#' information, by default "lat"
#' @param latitude The name of the variable in \code{data} with latitude
#' information, by default "lat"
#' @return A vector with the corresponding information (what) for the coordinates in data.
#' @export
#'
assign_region = function(data, what, longitude="lon", latitude="lat") {
  
  layer = switch(what,
                 longhurst = kali::longhurst,
                 LME       = kali::LME,
                 FAO       = kali::FAO,
                 GFMC      = kali::gsa,
                 stop(sprintf(what, "database not available.")))
  
  key   = switch(what,
                 longhurst = "ProvCode",
                 LME       = "LME_NUMBER",
                 GFMC      = "F_GSA_LIB",
                 FAO       = "zone",
                 stop(sprintf(what, "database not available.")))
  
  crs = proj4string(layer)
  pp = SpatialPoints(data[, c(longitude, latitude)], proj4string = CRS(crs))
  ind = over(pp, layer)[, key]
  return(ind)
  
}

#' Add known areas to a map
#'
#' @param what the database with the geographic areas. 
#' Currently, only "longhurst", "LME", "FAO" and "GFMC" are accepted.
#' @param subset The polygon within the database especified in 'what'.
#' @param col The color of the polygons.
#' @param ... Additional arguments to the plot function (currently not used)
#' @param list Boolean, list the available areas in the dataset.
#' @param plot 
#'
#' @return Invisibly, the polygon to plot.
#' @export
add_region = function(what, region=NULL, subset=NULL, col="red", list=FALSE, plot=TRUE, ...) {
  
  layer = switch(what,
                 longhurst = longhurst,
                 LME       = LME,
                 FAO       = FAO,
                 GFMC      = gsa,
                 stop(sprintf(what, "database not available.")))
  
  if(isTRUE(list)) return(layer@data)
  
  if(is.null(region)) {
    plot(x=layer, col=col, border=col, add=TRUE)
    return(invisible(layer))
  }
  
  ind = which(layer@data[, "region"] %in% region)
  plot(layer[ind, ], col=col, border=col, add=TRUE)
  
  return(invisible(layer[ind, ]))
  
}
