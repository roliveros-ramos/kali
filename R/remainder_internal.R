

# methods for 'rename' ----------------------------------------------------

rename.default = function(x, old, new, ...) {
  xname = deparse(substitute(x))
  if(is.null(names(x))) {
    stop(sprintf("Object %s has no names, nothing to do.", xname))
    return(x)
  }
  if(!(old %in% names(x))) {
    stop(sprintf("Object %s has no names, nothing to do.", xname))
  }
  names(x)[which(names(x) == old)] = new
  return(x)
}


.array2vector = function(z) {
  
  dim(z) = c(prod(dim(z)),1)
  
  return(z)
}


