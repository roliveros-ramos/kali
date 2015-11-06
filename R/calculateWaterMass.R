calculateWaterMass <- function(data, method = "oliveros", sst = "sst", sss = "sss", lon = "lon", 
                               lat = "lat", month = "month", dc = "dc", asFactors = TRUE){
  
  source("watermass_methods.R")
  
  output <- switch(tolower(method),
                   oliveros = .watermass_oliveros(data = data, sst = sst, sss = sss, lon = lon, lat = lat,
                                                 month = month, dc = dc, asFactors = asFactors),
                   swartzman = .watermass_swartzman(data = data, sst = sst, sss = sss, lat = lat, 
                                                   month = month, dc = dc, asFactors = asFactors))
  
  return(output)
}