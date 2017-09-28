
#' Get distribution records from on-line databases
#'
#' @param sp The scientific name of the species
#' @param limit The maximum number of records, NULL gets all records.
#' @param source The data source (currently only "gbif" and "obis" are supported.)
#' @param verbose boolean, do you want to know what's happening?
#' @param ... May be used later
#'
#' @return
#' @export
#'
#' @examples
#' getDistributionRecords(sp="Manta birostris", limit=300)
getDistributionRecords = function(sp, limit=NULL, source="gbif", verbose=FALSE, ...) {
  
  out = lapply(source, .getOccurrence, sp=sp, limit=limit, verbose=verbose)
  out = do.call(rbind, out)
  return(out)
  
}


# Internal ----------------------------------------------------------------

.getOccurrence = function(source, sp, limit=NULL, verbose=FALSE) {
  
  out = switch(source, 
               gbif = .getOccurrence_GBIF(sp=sp, limit=limit, verbose=verbose),
               obis = .getOccurrence_OBIS(sp=sp, limit=limit, verbose=verbose)
  )
  return(out)
  
}

.getOccurrence_OBIS = function(sp, limit=NULL, verbose=FALSE) {
  
  if(!requireNamespace("robis", quietly = TRUE)) 
    stop("You need to install the 'robis' package (from github).")
  
  if(length(sp)!=1) stop("You must provide only one species name.")
  sp = checkScientificName(sp)
  vars = c("decimalLongitude", "decimalLatitude", "year", "month", "day", "basisOfRecord")
  dat = occurrence(sp, verbose=verbose)
  dat$eventDate = as.Date(dat$eventDate)
  dat$year = year(dat$eventDate)
  dat$month = month(dat$eventDate)
  dat$day = day(dat$eventDate)
  out = dat[, vars]
  names(out)[1:2] = c("lon", "lat")
  out = out[complete.cases(out), ]
  nc = nchar(out$basisOfRecord)
  out$basisOfRecord[nc<2] = "Unknown"
  out = as_tibble(out)
  out = remove_rownames(out)
  out$basisOfRecord = as.factor(out$basisOfRecord)
  
  spName = names(which.max(table(dat$scientificName)))
  spAuthor = names(which.max(table(dat$scientificNameAuthorship)))
  
  attr(out, "sp") = spName
  attr(out, "author") = spAuthor
  attr(out, "sciName") = sprintf("%s %s", spName, spAuthor)
  attr(out, "accessed") = lubridate::today()
  attr(out, "source") = "OBIS"
  attr(out, "org") = "Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO"
  attr(out, "web") = "www.iobis.org"
  att = attributes(out)
  attr(out, "reference") = sprintf("%s (%d) Distribution records of %s %s [Dataset]. Available: %s. %s. Accessed: %s)\n", 
                                   att$source, year(att$accessed), att$sp, att$author, att$org, att$web, att$accessed)
  class(out) = c("occ_df", class(out))
  
  return(out)
  
}

.getOccurrence_GBIF = function(sp, limit=NULL, verbose=FALSE) {
  
  if(!requireNamespace("rgbif", quietly = TRUE)) 
    stop("You need to install the 'rgbif' package.")
  
  if(length(sp)!=1) stop("You must provide only one species name.")
  sp = checkScientificName(sp)
  vars = c("decimalLongitude", "decimalLatitude", "year", "month", "day", "basisOfRecord")
  tmp = rgbif::occ_search(scientificName = sp, limit=20)
  nrec = if(is.null(limit)) tmp$meta$count else limit
  dat = rgbif::occ_search(scientificName = sp, limit=nrec)$data
  out = dat[, vars]
  names(out)[1:2] = c("lon", "lat")
  out = out[complete.cases(out), ]
  out$basisOfRecord = toCamel(out$basisOfRecord, split="_")
  out$basisOfRecord[out$basisOfRecord=="Observation"] = "Occurrence"
  out = remove_rownames(out)
  out$basisOfRecord = as.factor(out$basisOfRecord)
  
  tmp = unlist(strsplit(names(which.max(table(dat$scientificName))), split=" "))
  spName = paste(tmp[1:2], collapse=" ")
  spAuthor = paste(tail(tmp, -2), collapse=" ")
  
  attr(out, "sp") = spName
  attr(out, "author") = spAuthor
  attr(out, "sciName") = sprintf("%s %s", spName, spAuthor)
  attr(out, "accessed") = lubridate::today()
  attr(out, "source") = "GBIF"
  attr(out, "org") = "GBIF. The Global Biodiversity Information Facility"
  attr(out, "web") = "www.gbif.org"
  att = attributes(out)
  attr(out, "reference") = sprintf("%s (%d) Distribution records of %s %s [Dataset]. Available: %s. %s. Accessed: %s)\n", 
                                   att$source, year(att$accessed), att$sp, att$author, att$org, att$web, att$accessed)
  
  
  class(out) = c("occ_df", class(out))
  
  return(out)
} 


# S3 methods --------------------------------------------------------------

print.occ_df = function(x, n=6, ...) {
  
  att = attributes(x)
  spHead = sprintf("\nDistribution records of %s\n\n", unique(att$sciName))
  
  cat(spHead)
  
  tibble:::print.tbl(x, n=n, ...)
  
  cat("\nCite this dataset as:\n")
  cat(paste(att$reference, collapse=""))
  
  return(invisible())
  
}

rbind.occ_df = function(..., deparse.level=1) {
  
  out = rbind.data.frame(..., deparse.level=deparse.level, 
                         stringsAsFactors=FALSE)
  
  attr(out, "sp")        = sapply(list(...), attr, which="sp")
  attr(out, "author")    = sapply(list(...), attr, which="author")
  attr(out, "sciName")   = sapply(list(...), attr, which="sciName")
  attr(out, "accessed")  = lapply(list(...), attr, which="accessed")
  attr(out, "source")    = lapply(list(...), attr, which="source")
  attr(out, "org")       = lapply(list(...), attr, which="org")
  attr(out, "web")       = lapply(list(...), attr, which="web")
  attr(out, "reference") = sapply(list(...), attr, which="reference")
  
  return(out)
  
}





