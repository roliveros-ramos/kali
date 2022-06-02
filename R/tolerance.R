


.getTolerance = function(x, what, tolerance=NULL) {
  x = subset(x, select=what, drop=TRUE)
  x = na.omit(x)
  q = quantile(x, prob=c(0.1, 0.25, 0.75, 0.9))
  iqd = q[3] - q[2]
  if(!is.null(tolerance)) {
    mina = min(tolerance)
    maxa = max(tolerance)
  } else {
    mina = min(q[2]-1.5*iqd, min(x))
    maxa = max(q[3]+1.5*iqd, max(x))
  }
  out = c(mina, q[1], mean(x), median(x), q[4], maxa)
  names(out) = c("minimum", "prefered_min","prefered_mean", "prefered_med", "prefered_max", "maximum")
  return(out)
}

#' Estimate Tolerance to an environmental factor
#'
#' @param x 
#' @param what 
#' @param split 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
getTolerance = function(x, what, split=NULL, ...) {
  
  tolerance = .getTolerance(x=x, what=what)
  
  if(is.null(split)) return(tolerance)
  
  f = subset(x, select=split, drop=TRUE)
  
  dat = split(x, f = f)
  
  out = lapply(dat, FUN=.getTolerance, what=what, tolerance=tolerance)
  
  out = c(out, '_all_levels'=list(tolerance))
  return(out)
  
}
