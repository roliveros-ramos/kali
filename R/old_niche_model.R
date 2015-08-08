
# remove ;;, to be deprecated
.sub = function(map) {
  
  map = gsub(";", "", gsub(";;", ",", map))
  map = strsplit(map, ",")
  map = sapply(map, as.numeric)
  
  return(map)
}

