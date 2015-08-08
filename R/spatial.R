checkLongitude = function(x, ...) {
  if (!any(x > 180)) 
    return(x)
  x[x > 180] = x[x > 180] - 360
  return(x)
}