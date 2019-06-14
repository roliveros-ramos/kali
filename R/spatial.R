.longitude2Center = function(x, ...) {
  if (!any(x > 180)) 
    return(x)
  x[which(x > 180)] = x[which(x > 180)] - 360
  return(x)
}

.longitude2Left = function(x, ...) {
  if (!any(x < 0)) 
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

