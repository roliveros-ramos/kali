
#' Get distribution records from on-line databases
#'
#' @param sp The scientific name of the species
#' @param limit The maximum number of records, NULL gets all records.
#' @param source The data source (currently only "gbif" and "obis" are supported.)
#' @param year Range of years to be included in the search.
#' @param verbose boolean, do you want to know what's happening?
#' @param ... May be used later
#'
#' @return
#' @export
#'
#' @examples
#' getDistributionRecords(sp="Manta birostris", limit=300)
getDistributionRecords = function(sp, limit=NULL, source="gbif", year=NULL, 
                                  verbose=FALSE, ...) {
  
  out = list()
  for(isource in source) {
    message(sprintf("\nDownloading occurrence data from %s.\n", isource))
    out[[isource]] = .getOccurrence(sp=sp, source=isource, limit=limit, verbose=verbose, 
                 year=year)
  }
  out = do.call(rbind, out)
  check = setdiff(names(out), "provider") # check for duplicated up to provider.
  ind = duplicated(as.data.frame(out[, check]))
  ndup = sum(ind, na.rm=TRUE)
  if(ndup>0) message(sprintf("\nDetecting up to %d duplicated records, check the dataset before use.", ndup))
  return(out)
  
}


# S3 methods --------------------------------------------------------------

#' @export
print.occ_df = function(x, n=NULL, ...) {
  
  attList = attributes(x)
  spHead = sprintf("\nDistribution records of %s\n\n", unique(attList$sciName))
  
  cat(spHead)
  
  tibble:::print.tbl_df(x, n=n, ...)
  
  cat("\nCite this dataset as:\n")
  cat(paste(attList$reference, collapse=""))
  
  return(invisible(NULL))
  
}

#' @export
rbind.occ_df = function(..., deparse.level=1) {
  
  out = rbind.data.frame(..., deparse.level=deparse.level, 
                         stringsAsFactors=FALSE)
  
  attr(out, "sp")        = sapply(list(...), attr, which="sp")
  attr(out, "author")    = sapply(list(...), attr, which="author")
  attr(out, "sciName")   = sapply(list(...), attr, which="sciName")
  attr(out, "accessed")  = sapply(list(...), attr, which="accessed")
  attr(out, "source")    = sapply(list(...), attr, which="source")
  attr(out, "org")       = sapply(list(...), attr, which="org")
  attr(out, "web")       = sapply(list(...), attr, which="web")
  attr(out, "reference") = sapply(list(...), attr, which="reference")
  
  return(out)
  
}

#' @export
'[.occ_df' = function(x, i, j, drop = FALSE) {
  
  attList = attributes(x)
  # out = NextMethod("[")
  class(x) =  tail(class(x), -1)
  out = tibble:::'[.tbl_df'(x, i, j, drop)
  attList$row.names = attr(out, "row.names")
  attList$names = attr(out, "names")
  attributes(out) = attList
  return(out)
}

# head.occ_df = function (x, n = 6L, ...) {
#     # class(x) =  tail(class(x), -1)
#     stopifnot(length(n) == 1L)
#     n <- if (n < 0L) 
#       max(nrow(x) + n, 0L)
#     else min(n, nrow(x))
#     out = x[seq_len(n), , drop = FALSE]
#     # class(out) = c("occ_df", class(x))
#     return(out)
# }
# 
# tail.occ_df = function (x, n = 6L, ...) {
#   class(x) =  tail(class(x), -1)
#   stopifnot(length(n) == 1L)
#   nrx <- nrow(x)
#   n <- if (n < 0L) 
#     max(nrx + n, 0L)
#   else min(n, nrx)
#   x[seq.int(to = nrx, length.out = n), , drop = FALSE]
# }

# Internal ----------------------------------------------------------------

.getOccurrence = function(source, sp, limit=NULL, year=NULL, verbose=FALSE) {
  
  out = switch(source, 
               gbif = .getOccurrence_GBIF(sp=sp, limit=limit, year=year, 
                                          verbose=verbose),
               obis = .getOccurrence_OBIS(sp=sp, limit=limit, year=year, 
                                          verbose=verbose)
  )
  return(out)
  
}

.getOccurrence_OBIS = function(sp, limit=NULL, year=NULL, verbose=FALSE) {
  
  if(!requireNamespace("robis", quietly = TRUE)) 
    stop("You need to install the 'robis' package (from github).")
  
  if(length(sp)!=1) stop("You must provide only one species name.")
  sp = check_taxon(sp)
  vars = c("decimalLongitude", "decimalLatitude", "year", "month", "day", 
           "depth", "lifeStage", "basisOfRecord")
  
  startdate = NULL
  enddate   = NULL
  
  if(!is.null(year)) {
    if(!is.numeric(year)) stop("Argument 'year' must be numeric")
    year = na.omit(year)
    year = range(c(year, floor(min(year)), ceiling(max(year))), na.rm=TRUE)
    startdate = as.Date(sprintf("%d-01-01", year[1]))
    enddate   = as.Date(sprintf("%d-12-31", year[2]))
  }
  
  
  dat = robis::occurrence(sp, verbose=verbose, limit=limit,
                          startdate=startdate, enddate=enddate)
  if(nrow(dat)<1) {
  msg = sprintf("\nRetrieved %d records of %d (%0.2f%%)\n", 0, 0, 100)
  cat(msg)
    return(NULL)
    }
  dat$eventDate = as.Date(dat$eventDate)
  dat$year = year(dat$eventDate)
  dat$month = month(dat$eventDate)
  dat$day = day(dat$eventDate)
  out = dat[, vars]
  names(out)[1:2] = c("lon", "lat")
  n = nrow(out)
  out = out[complete.cases(out[, c("year", "month", "day")]), ]
  n = n - nrow(out)
  if(n>0) message(sprintf("\nRemoved %d records without date information.", n))
  nc = nchar(out$basisOfRecord)
  out$basisOfRecord[nc<2] = "Unknown"
  out$basisOfRecord = toCamel(out$basisOfRecord, split="_")
  out$basisOfRecord[out$basisOfRecord=="Observation"] = "Occurrence" # why?
  out$lifeStage[out$lifeStage=="Unknown"] = NA
  out = as_tibble(out)
  out = tibble::remove_rownames(out)
  out$basisOfRecord = as.factor(out$basisOfRecord)
  out$provider = "OBIS"
  spName = names(which.max(table(dat$scientificName)))
  spAuthor = names(which.max(table(dat$scientificNameAuthorship)))
  
  attr(out, "sp") = spName
  attr(out, "author") = spAuthor
  attr(out, "sciName") = sprintf("%s %s", spName, spAuthor)
  attr(out, "accessed") = lubridate::today()
  attr(out, "source") = "OBIS"
  attr(out, "org") = "Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO"
  attr(out, "web") = "www.iobis.org"
  attr(out, "all_vars") = names(dat)
  att = attributes(out)
  attr(out, "reference") = sprintf("%s (%d) Distribution records of %s %s [Dataset]. Available: %s. %s. Accessed: %s)\n", 
                                   att$source, year(att$accessed), att$sp, att$author, att$org, att$web, att$accessed)
  class(out) = c("occ_df", class(out))
  
  return(out)
  
}

.getOccurrence_GBIF = function(sp, limit=NULL, year=NULL, verbose=FALSE) {
  
  if(!requireNamespace("rgbif", quietly = TRUE)) 
    stop("You need to install the 'rgbif' package.")
  
  if(length(sp)!=1) stop("You must provide only one species name.")
  sp = check_taxon(sp)
  vars = c("decimalLongitude", "decimalLatitude", "year", "month", "day", 
           "depth", "lifeStage", "basisOfRecord")
  
  oyear = year
  
  if(!is.null(year)) {
    if(!is.numeric(year)) stop("Argument 'year' must be numeric")
    year = na.omit(year)
    year = range(c(year, floor(min(year)), ceiling(max(year))), na.rm=TRUE)
    year = sprintf("%d,%d", year[1], year[2])
  }
  
  tmp = rgbif::occ_search(scientificName = sp, limit=0, year=year, hasCoordinate=TRUE)
  ntot = tmp$meta$count
  if(ntot==0) {
    msg = sprintf("\nRetrieved %d records of %d (%0.2f%%)\n", 0, 0, 100)
    cat(msg)
    return(NULL)
  }

  nrec = if(is.null(limit)) ntot else min(limit, ntot)  
  
  if(nrec >= 1e5) {
    
    ncount = .getOccNumber_GBIF(sp, year=oyear)
    split = attr(ncount, "split")
    split[1] = split[1] - 1 
    nquery = length(split) - 1
    
    output = list()
    for(i in seq_len(nquery)) {
      
      iyear = split[i + c(0, 1)] + c(1, 0)
      msg = sprintf("Retrieving records from %d to %d...", iyear[1], iyear[2])
      message(msg)
      iyear = sprintf("%d,%d", iyear[1], iyear[2])
      tmp = rgbif::occ_search(scientificName = sp, limit=1e5, year=iyear, hasCoordinate=TRUE)$data
      msg = sprintf("\tRetrieved %d records (%0.2f%%).\n", nrow(tmp), 100*nrow(tmp)/ntot)
      message(msg)
     
      tmp = tmp[, c(vars, "scientificName")]
      output[[i]] = tmp
      
    }
    
    dat = do.call(rbind, output)
    
  } else {
    
    dat = rgbif::occ_search(scientificName = sp, limit=nrec)$data
    msg = sprintf("\nRetrieved %d records of %d (%0.2f%%)\n", nrec, ntot, 100*nrec/ntot)
    message(msg)
    
  }

  out = dat[, vars]
  names(out)[1:2] = c("lon", "lat")
  n = nrow(out)
  out = out[complete.cases(out[, c("year", "month", "day")]), ]
  n = n - nrow(out)
  if(n>0) message(sprintf("Removed %d records without date information.", n))
  out$basisOfRecord = toCamel(out$basisOfRecord, split="_")
  out$basisOfRecord[out$basisOfRecord=="Observation"] = "Occurrence" # why?
  out$lifeStage[out$lifeStage=="Unknown"] = NA # why?
  out = tibble::remove_rownames(out)
  out$basisOfRecord = as.factor(out$basisOfRecord)
  out$provider = "GBIF"
  
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
  attr(out, "all_vars") = names(dat)
  att = attributes(out)
  attr(out, "reference") = sprintf("%s (%d) Distribution records of %s %s [Dataset]. Available: %s. %s. Accessed: %s)\n", 
                                   att$source, year(att$accessed), att$sp, att$author, att$org, att$web, att$accessed)
  
  
  class(out) = c("occ_df", class(out))
  
  return(out)
} 


# Count occurrences -------------------------------------------------------

.getOccNumber_GBIF = function(sp, year=NULL, hasCoordinate=TRUE) {
  
  if(!requireNamespace("rgbif", quietly = TRUE)) 
    stop("You need to install the 'rgbif' package.")
  
  if(length(sp)!=1) stop("You must provide only one species name.")
  sp = check_taxon(sp)
  
  if(is.null(year)) year = c(1900, lubridate::year(lubridate::today()))
  if(!is.numeric(year)) stop("Argument 'year' must be numeric")
  year = na.omit(year)
  year = range(c(year, floor(min(year)), ceiling(max(year))), na.rm=TRUE)
  year = seq(from=year[1], to=year[2])
  
  out = rep(NA_integer_, length=length(year))
  
  for(i in seq_along(year)) {
    out[i] = rgbif::occ_search(scientificName = sp, limit=0, year=year[i], hasCoordinate=TRUE)$meta$count
  }
  names(out) = year
  
  tot = sum(out, na.rm=TRUE)
  
  nquery = ceiling(tot/1e5) + 1
  crit = tot/nquery
  ctot = cumsum(out)
  
  split = numeric(nquery+1)
  split[1] = 1
  split[nquery+1] = length(year)
  for(i in seq_len(nquery-1)) split[i+1] = which.max(ctot > (crit*i)) - 1
  split = year[split]
  
  attr(out, "split") = split
  
  return(out)
  
}


