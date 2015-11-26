
# Swartzman method -----------------------------------------------------------------------

.watermass_swartzman <- function(data, sst, sss, lat, month, dc, asFactors){
  
  # Select just useable variables
  data <- data[,c(sst, sss, dc, lat, month)]
  
  # If there is not any complete row, out NA vector
  if(sum(complete.cases(data)) == 0) return(output)
  
  # Define breaks values
  # dcBreak <- 350*1.852 # 350 nautical miles to Km
  latBreaks <- c(-Inf, -13, -8, Inf)
  sstBreaks <- c(-Inf, 13, 13.5, 14, 17, 18, 19, 20, 21, 22, 24, 25, 26, Inf)
  sssBreaks <- c(-Inf, 34.0, 34.8, 35.05, 35.1, Inf)
  seasonBreaks <- c(0, 3, 6, 9, 12)
  dcBreaks <- c(-Inf, 5, Inf)
  
  # Making categorizations
  latCat <- cut(x = data$lat, breaks = latBreaks)
  sstCat <- cut(x = data$sst, breaks = sstBreaks)
  sssCat <- cut(x = data$sss, breaks = sssBreaks)
  seasonCat <- cut(x = data$month, breaks = seasonBreaks, 
                   labels = c("summer", "autumn", "winter", "spring"))
  dcCat <- cut(x = data$dc, breaks = dcBreaks)
  
  # Show levels and values of each category
  categories <- c("sstCat", "sssCat", "latCat", "seasonCat", "dcCat")
  #   for(i in categories){
  #     tempVector <- levels(get(i))
  #     names(tempVector) <- seq_along(tempVector)
  #     
  #     cat(paste0("\n\n Category: ", i, "\n"))
  #     print(tempVector)
  #   }
  
  # Set Conditions in paper:            #SST  #SSS  #Lat  #Season  #DC
  wmDefinitions <- list(ccw = list(list(4:5, 3, 1:2, 4, 1:2),
                                   list(4:6, 3, 1:2, 1, 1:2),
                                   list(4:5, 3, 1:2, 2:4, 1:2),
                                   list(3:4, 3, 1:2, 4, 1:2),
                                   list(6:13, 3, 3, 1, 1:2),
                                   list(6:13, 3, 3, 2:3, 1:2),
                                   list(5:13, 3, 3, 2:3, 1:2)),
                        
                        ssw = list(list(5:13, 5, 3, 2:3, 1:2)),
                        
                        sew = list(list(8:12, 2, 3, 2:3, 1:2),
                                   list(9:12, 2, 3, 3:4, 1:2)),
                        
                        stw = list(list(9:13, 1, 3, 4, 1:2),
                                   list(11:13, 1, 3, 1, 1:2),
                                   list(10:13, 1, 3, 2:3, 1:2),
                                   list(8:13, 1, 3, 2:3, 1:2)),
                        
                        mcs = list(list(4:11, 4, 3, c(2, 4), 1:2),
                                   list(3:11, 4, 3, 3, 1:2),
                                   list(4:5, 4:5, 3, 4, 1:2),
                                   list(4:6, 4:5, 3, 1, 1:2),
                                   list(4:5, 4:5, 3, 2:3, 1:2),
                                   list(3:4, 4:5, 3, 2:3, 1:2)),
                        
                        mesc = list(list(6:13, 3, 3, c(2, 4), 1:2),
                                    list(7:13, 3, 3, 1, 1:2),
                                    list(5:13, 3, 3, 3, 1:2),
                                    list(4:5, 2, 3, c(2, 4), 1:2),
                                    list(4, 2, 3, 1, 1:2)),
                        
                        mrw = list(list(4:7, 1:2, 1:2, c(2, 4), 1),
                                   list(3:7, 1:2, 1:2, 3, 1)),
                        
                        caw = list(list(2:5, 1:2, 1, 3, 2)))
  
  # Get combinations on big data frame
  allDefinitions <- NULL
  for(i in seq_along(wmDefinitions)){
    tempDefinitions_1 <- wmDefinitions[[i]]
    
    for(j in seq_along(tempDefinitions_1)){
      tempDefinitions_2 <- tempDefinitions_1[[j]]
      
      tempDefinitions_3 <- do.call("expand.grid", tempDefinitions_2)
      colnames(tempDefinitions_3) <- categories
      
      tempDefinitions_3$wm <- names(wmDefinitions)[i]
      
      allDefinitions <- rbind(allDefinitions, tempDefinitions_3)
    }
  }
  
  # Set water mass definition in format sst-sss-lat-season-dc
  allDefinitions <- data.frame(wm = allDefinitions$wm,
                               definition = with(allDefinitions, 
                                                 paste(sstCat, sssCat, latCat, 
                                                       seasonCat, dcCat, sep = "-")),
                               stringsAsFactors = FALSE)
  
  # Set categorized values in format sst-sss-lat-season-dc
  data <- data.frame(sstCat = as.numeric(sstCat),
                     sssCat = as.numeric(sssCat),
                     latCat = as.numeric(latCat),
                     seasonCat = as.numeric(seasonCat),
                     dcCat = as.numeric(dcCat))
  # data$wm <- apply(data[,categories], 1, paste, collapse = "-")
  data$wm <- with(data, paste(sstCat, sssCat, latCat, seasonCat, dcCat, sep = "-"))
  
  # Get water mass for each row
  output <- allDefinitions$wm[match(data$wm, allDefinitions$definition)]
  
  # convert output as factor
  if(isTRUE(asFactors))
    output <- factor(x = toupper(output), 
                     levels = c("CCW", "SSW", "SEW", "STW", 
                                "MCS", "MESC", "MRW", "CAW"))
  
  return(output)
}



# Oliveros method ------------------------------------------------------------------------

.watermass_oliveros <- function(data, sst, sss, lat, lon, month, dc, asFactors){
  labels <- c("ASA", "ATS", "AMSA", "AES", "ACF", "ASS")
  
  data = data[, c(sst, sss, dc, lat, lon, month)]
  output = rep(NA, nrow(data))
  noNA = complete.cases(data)
  if(sum(noNA)==0) return(output)
  
  data = data[noNA, ]
  
  # Define breaks values
  dcBreak = 350*1.852 # 350 nautical miles to Km
  latBreak = c(-25, -12)
  sstBreaks = c(0, 14, 18, 20, 23, 26, 28, 40)
  sssBreaks = c(25, 34, 34.8, 35.05, 40)
  recatFactor = c(0, 8, 9, 12, 16, 20, 24, 28)
  
  # add validation for NAs in categories
  latCat  = numeric(nrow(data))
  sstCat  = as.numeric(cut(data[, sst], sstBreaks, labels = seq(length(sstBreaks) - 1)))
  sssCat  = as.numeric(cut(data[, sss], sssBreaks, labels = seq(length(sssBreaks) - 1)))
  cdiCat  = as.numeric(cut(data[, dc], c(0, dcBreak, Inf), labels = seq(length(dcBreak) + 1)))
  latCat1 = as.numeric(cut(data[, lat], c(-60.1, latBreak[2], 10.1), labels = seq(length(latBreak))))
  latCat2 = as.numeric(cut(data[, lat], c(-60.1, latBreak[1], 10.1), labels = seq(length(latBreak))))
  latCat[data[, lon] > -140 & data[, lon] < -60] = latCat1[data[, lon] > -140 & data[, lon] < -60]
  latCat[data[, lon] <= -140 | data[, lon] >= -60] = latCat2[data[, lon] <= -140 | data[, lon] >= -60]
  monthCat = data[, month] < 4
  
  # Create vectors to first categorization (only works for length(cat)<10)
  catLabels = as.numeric(apply(cbind(sstCat, sssCat), 1, paste, collapse = ""))
  patLabels = as.numeric(apply(expand.grid(seq(length(sstBreaks) - 1), seq(length(sssBreaks) - 1)), 1,
                               paste, collapse = ""))
  
  # Make first categorization and recategorization
  
  waterMass = as.numeric(cut(match(catLabels, patLabels), breaks=recatFactor, 
                             labels = seq_len(length(recatFactor)-1)))
  
  # Create vector to second categorization
  
  catLabels = as.numeric(apply(data.frame(cdiCat, latCat, waterMass), 1, paste, collapse = ""))  
  patLabels = as.numeric(apply(expand.grid(seq_len(length(dcBreak) + 1), 
                                           seq_len(length(latBreak)), 
                                           seq_along(recatFactor)), 1, paste, collapse = ""))
  
  label = c(1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 4, 4, 7, 3, 
            8, 4, 6, 6, 6, 6, 5, 6, 5, 6, 6, 6, 6, 6)
  
  # Make second categorization
  waterMass = label[match(catLabels, patLabels)]
  
  # Make third categorization
  ind = (waterMass == 7)
  waterMass[ind] = 5 + as.numeric(monthCat[ind])
  ind = (waterMass == 8)
  waterMass[ind] = 5 - as.numeric(monthCat[ind])
  
  output[noNA] = waterMass
  if(isTRUE(asFactors)) 
    output = factor(x = labels[output], levels = labels)
  
  return(output)
}


# Zuta method -----------------------------------------------------------------------

.watermass_zuta <- function(data, sst, sss, lat, depth, asFactors){
  
  # If there is not value for depth
  if(is.null(data$depth)){
    depth <- "depth"
    data$depth <- 0.1
  }
  
  # Select just useable variables
  data <- data[,c(sst, sss, depth, lat)]
  
  # If there is not any complete row, out NA vector
  if(sum(complete.cases(data)) == 0) return(output)
  
  # Define breaks values
  # dcBreak <- 350*1.852 # 350 nautical miles to Km
  latBreaks <- c(-Inf, -8, Inf)
  sstBreaks <- c(-Inf, 14, 18, 20, 25, 27, Inf)
  sssBreaks <- c(-Inf, 33.8, 34.8, 35.1, 35.7, Inf)
  depthBreaks <- c(0, 20, 40, 50, 80, 100, 120, Inf)
  
  # Making categorizations
  latCat <- cut(x = data$lat, breaks = latBreaks)
  sstCat <- cut(x = data$sst, breaks = sstBreaks)
  sssCat <- cut(x = data$sss, breaks = sssBreaks)
  depthCat <- cut(x = data$depth, breaks = depthBreaks)
  
  # Show levels and values of each category
  categories <- c("sstCat", "sssCat", "latCat", "depthCat")
  #   for(i in categories){
  #     tempVector <- levels(get(i))
  #     names(tempVector) <- seq_along(tempVector)
  #     
  #     cat(paste0("\n\n Category: ", i, "\n"))
  #     print(tempVector)
  #   }
  
  # Set Conditions in paper:            #SST  #SSS  #Lat  #Season  #DC
  wmDefinitions <- list(ats = list(list(5:6, 1, 1:2, 1)),
                        
                        aes = list(list(4:6, 2, 1:2, 1:2)),
                        
                        ass = list(list(3:5, 4, 1:2, 1:5)),
                        
                        acf = list(list(2, 3, 2, 1:4),
                                   list(1:2, 3, 1, 1:6)),
                        
                        mcs = list(list(2, 4:5, 1:2, 1:5)),
                        
                        mesc = list(list(3:6, 3, 1:2, 1:3),
                                    list(2, 2, 1:2, 1:3)))
  
  # Get combinations on big data frame
  allDefinitions <- NULL
  for(i in seq_along(wmDefinitions)){
    tempDefinitions_1 <- wmDefinitions[[i]]
    
    for(j in seq_along(tempDefinitions_1)){
      tempDefinitions_2 <- tempDefinitions_1[[j]]
      
      tempDefinitions_3 <- do.call("expand.grid", tempDefinitions_2)
      colnames(tempDefinitions_3) <- categories
      
      tempDefinitions_3$wm <- names(wmDefinitions)[i]
      
      allDefinitions <- rbind(allDefinitions, tempDefinitions_3)
    }
  }
  
  # Set water mass definition in format sst-sss-lat-season-dc
  allDefinitions <- data.frame(wm = allDefinitions$wm,
                               definition = with(allDefinitions, 
                                                 paste(sstCat, sssCat, latCat, 
                                                       depthCat, sep = "-")),
                               stringsAsFactors = FALSE)
  
  # Set categorized values in format sst-sss-lat-season-dc
  data <- data.frame(sstCat = as.numeric(sstCat),
                     sssCat = as.numeric(sssCat),
                     latCat = as.numeric(latCat),
                     depthCat = as.numeric(depthCat))
  # data$wm <- apply(data[,categories], 1, paste, collapse = "-")
  data$wm <- with(data, paste(sstCat, sssCat, latCat, depthCat, sep = "-"))
  
  # Get water mass for each row
  output <- allDefinitions$wm[match(data$wm, allDefinitions$definition)]
  
  # convert output as factor
  if(isTRUE(asFactors))
    output <- factor(x = toupper(output), 
                     levels = c("ATS", "AES", "ASS", "ACF", "MCS", "MESC"))
  
  return(toupper(output))
}