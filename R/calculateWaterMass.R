calculateWaterMass <- function(data, method = "oliveros", sst = "sst", sss = "sss", lon = "lon", 
                               lat = "lat", month = "month", depth = "depth", dc = "dc", 
                               asFactors = TRUE){
  
  output <- switch(tolower(method),
                   oliveros = .watermass_oliveros(data = data, sst = sst, sss = sss, lon = lon, 
                                                  lat = lat, month = month, dc = dc, 
                                                  asFactors = asFactors),
                   swartzman = .watermass_swartzman(data = data, sst = sst, sss = sss, lat = lat, 
                                                    month = month, dc = dc, asFactors = asFactors),
                   zuta = .watermass_zuta(data = data, sst = sst, sss = sss, lat = lat, 
                                          depth = depth, asFactors = asFactors))
  
  return(output)
}