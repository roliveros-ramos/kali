# # Calculate Watermass category
# # ASA  (1): Aguas Superficiales Antárticas
# # ATS  (2): Aguas Tropicales Superficiales
# # AMSA (3): Aguas de Mezcla Subtropical-Antártica
# # AES  (4): Aguas Ecuatoriales Superficiales
# # ACF  (5): Aguas Costeras Frías
# # ASS  (6): Aguas Subtropicales Superficiales
# 
# # INPUTS
# # Lon and Lat in Decimal degrees
# # Distance to coast in Km
# 
# calculateWaterMass = function(data, labels = c("ASA", "ATS", "AMSA", "AES", "ACF", "ASS"),
#                               asFactors=TRUE,
#                               sst="sst", sss="sss", dc="dc", lat="lat", lon="lon", month="month") {
#   
#   data = data[, c(sst, sss, dc, lat, lon, month)]
#   output = rep(NA, nrow(data))
#   noNA = complete.cases(data)
#   if(sum(noNA)==0) return(output)
#   
#   data = data[noNA, ]
#   
#   # Define breaks values
#   dcBreak = 350*1.852 # 350 nautical miles to Km
#   latBreak = c(-25, -12)
#   sstBreaks = c(0, 14, 18, 20, 23, 26, 28, 40)
#   sssBreaks = c(25, 34, 34.8, 35.05, 40)
#   recatFactor = c(0, 8, 9, 12, 16, 20, 24, 28)
#   
#   # add validation for NAs in categories
#   latCat  = numeric(nrow(data))
#   sstCat  = as.numeric(cut(data[, sst], sstBreaks, labels = seq(length(sstBreaks) - 1)))
#   sssCat  = as.numeric(cut(data[, sss], sssBreaks, labels = seq(length(sssBreaks) - 1)))
#   cdiCat  = as.numeric(cut(data[, dc], c(0, dcBreak, Inf), labels = seq(length(dcBreak) + 1)))
#   latCat1 = as.numeric(cut(data[, lat], c(-60.1, latBreak[2], 10.1), labels = seq(length(latBreak))))
#   latCat2 = as.numeric(cut(data[, lat], c(-60.1, latBreak[1], 10.1), labels = seq(length(latBreak))))
#   latCat[data[, lon] > -140 & data[, lon] < -60] = latCat1[data[, lon] > -140 & data[, lon] < -60]
#   latCat[data[, lon] <= -140 | data[, lon] >= -60] = latCat2[data[, lon] <= -140 | data[, lon] >= -60]
#   monthCat = data[, month] < 4
#   
#   # Create vectors to first categorization (only works for length(cat)<10)
#   catLabels = as.numeric(apply(cbind(sstCat, sssCat), 1, paste, collapse = ""))
#   patLabels = as.numeric(apply(expand.grid(seq(length(sstBreaks) - 1), seq(length(sssBreaks) - 1)), 1,
#                                 paste, collapse = ""))
#   
#   # Make first categorization and recategorization
#   
#   waterMass = as.numeric(cut(match(catLabels, patLabels), breaks=recatFactor, 
#                              labels = seq_len(length(recatFactor)-1)))
#   
#   # Create vector to second categorization
#   
#   catLabels = as.numeric(apply(data.frame(cdiCat, latCat, waterMass), 1, paste, collapse = ""))  
#   patLabels = as.numeric(apply(expand.grid(seq_len(length(dcBreak) + 1), 
#                                            seq_len(length(latBreak)), 
#                                            seq_along(recatFactor)), 1, paste, collapse = ""))
#   
#   label = c(1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 4, 4, 7, 3, 
#             8, 4, 6, 6, 6, 6, 5, 6, 5, 6, 6, 6, 6, 6)
#   
#   # Make second categorization
#   waterMass = label[match(catLabels, patLabels)]
#   
#   # Make third categorization
#   ind = (waterMass == 7)
#   waterMass[ind] = 5 + as.numeric(monthCat[ind])
#   ind = (waterMass == 8)
#   waterMass[ind] = 5 - as.numeric(monthCat[ind])
#   
#   output[noNA] = waterMass
#   if(isTRUE(asFactors)) output = factor(x=labels[output], levels=labels)
#   
#   return(output)
# }
# 
