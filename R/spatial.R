.longitude2Center = function(x, ...) {
  if (!any(x > 180)) 
    return(x)
  x[x > 180] = x[x > 180] - 360
  return(x)
}

.longitude2Left = function(x, ...) {
  if (!any(x < 0)) 
    return(x)
  x[x < 0] = x[x < 0] + 360
  return(x)
}

checkLongitude = function(x, primeMeridian="center", ...) {
  out = switch(primeMeridian, 
               center = .longitude2Center(x, ...),
               left = .longitude2Left(x, ...))
  return(out)
}

.findPrimeMeridian = function(x) {
  if(any(x<0)) return("center")
  if(any(x>180)) return("left")
  return("center")
}

