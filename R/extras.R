getTime = function(data) {
  names(data) = tolower(names(data))
  if(all(c("year", "month", "day") %in% names(data)))
    return(with(data, round(year + (month-1)/12 + (day-1)/365, 5)))
  if(all(c("year", "month") %in% names(data)))
    return(with(data, round(year + (month-1)/12 + (15-1)/365, 5)))
  if("date" %in% names(data)) {
    dates = ymd(data$date)
    year  = year(dates)
    month = month(dates)
    day   = day(dates)
    return(round(year + (month-1)/12 + (day-1)/365, 5))
  }
  stop("Cannot create time from data.")
}


# trimGAM = function(x) {
#   
#   x$model = NULL
#   x$prior.weights = NULL
#   x$offset = NULL
#   x$weights = NULL
#   x$linear.predictors = NULL #
#   x$db.drho = NULL
#   x$hat = NULL             
#   x$mgcv.conv = NULL       
#   x$rank = NULL            
#   x$Ve = NULL   #           
#   x$Vp = NULL              
#   x$Vc = NULL   #           
#   x$scale.estimated = NULL              
#   x$outer.info = NULL              
#   x$optimizer = NULL              
#   x$scale = NULL              
#   x$sig2 = NULL              
#   x$sp = NULL              
#   x$gamma = NULL
#   x$G = NULL    #
#   x$qrx = NULL  #
#   x$iter =  NULL
#   x$assign = NULL
#   x$boundary = NULL
#   x$control = NULL
#   x$converged = NULL
#   
#   return(x)
# }

trimGAM = function(object) {
  object$Sl <- object$qrx <- object$R <- object$F <- NULL
  object$Ve <- object$Vc <- object$G <- object$residuals <- NULL
  object$fitted.values <- object$linear.predictors <- NULL
  return(object)
}

num2date = function(x) {
  # requires(lubridate)
  year = floor(x)
  days = x%%1*(365+leap_year(year)) + 1
  dates = as.Date(paste(year, "01", "01", sep="-"))
  yday(dates) = floor(round(days, 6))
  return(dates)
}


load_object = function(file, object="env") {
  name = paste("tmp", basename(tempfile()), sep="_")
  attach(what=file, name=name)
  on.exit(detach(name, character.only = TRUE))
  if(!exists(object, where=2, inherits=FALSE)) 
    stop("object not found.")
  return(get(object, pos=2, inherits=FALSE))
}

toCamel = function(x, split=" ", lower=FALSE) {
  .capitalize = function(x) {
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    return(x)
  }
  x = strsplit(x, split=split)
  x = lapply(x, tolower)
  x = lapply(x, .capitalize)
  x = sapply(x, paste, collapse="")
  if(isTRUE(lower)) substr(x, 1, 1) = tolower(substr(x, 1, 1))
  return(x)
}

check_taxon = function(x, na.return=FALSE, verbose=TRUE) {
  
  if(!requireNamespace("taxize", quietly = TRUE)) 
    stop("You need to install the 'taxize' package.")
  
  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")
  
  .checkScientificName = function(x) {
    tmp = taxize::gnr_resolve(names = x, canonical = TRUE)$matched_name2
    tmp = names(which.max(table(tmp)))
    isNA = FALSE
    if(is.null(tmp)) {
      if(isTRUE(na.return)) {
        if(isTRUE(verbose)) message(sprintf("Name '%s' not found, returning NA.", x))
        tmp = NA_character_
        isNA = TRUE
      } else {
        if(isTRUE(verbose)) message(sprintf("Name '%s' not found, returning original name.", x))
        tmp = x 
      }
    }
    id = identical(x, tmp)
    if(!id & !isNA & isTRUE(verbose)) message(sprintf("Species '%s' was corrected to '%s'.", x, tmp))
    return(tmp)
  }
  out = unlist(sapply(x, .checkScientificName))
  if(is.null(out)) {
    
  }
  return(out)
}

get_taxon = function(x, rank, db="itis") {
  
  if(length(rank)!=1) stop("Only one taxon is allowed.")
  
  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")
  
  db = match.arg(db, c("itis", "ncbi", "both"))
  isna = is.na(x)
  tmp = tax_name(query=x[!isna], db=db, get=rank, messages=FALSE, ask=FALSE)
  out = character(length(x))
  out[which(!isna)] = tmp[, rank]
  return(out)
  
}

#' Weighted random sampling with a reservoir
#'
#' Implementation of the Weighted random sampling with a reservoir (without replacement)
#' (Efraimidis & Spirakis, 2006) algorithm.
#' @param x a vector of one or more elements from which to choose.
#' @param size a non-negative integer giving the number of items to choose.
#' @param prob a vector of weights for obtaining the elements of the vector being sampled.
#'
#' @references Efraimidis & Spirakis (2006). Weighted random sampling with a reservoir
#' @return A vector of length \code{size} with elements drawn from \code{x}
#' @export
#'
#' @examples
#' N = 1000
#' x = seq_len(N)
#' prob = c(rep(0.1, N/2), rep(1, N/2))
#' x_sample = sample.weighted(x=x, prob=prob, size=N/2)
#' hist(x_sample)
sample.weighted = function(x, size, prob) {
  if(length(x) == 1L && is.numeric(x) && is.finite(x) && x >= 
      1) x = seq_len(x)
  if (missing(size)) size = length(x)
  if(length(prob)!=length(x)) stop("'prob' must match the length of 'x'.")
  if(size>length(x)) stop("cannot take a sample larger than the population.")
  nr = runif(length(x))
  ind = order(log(nr)/prob, decreasing=TRUE)[seq_len(size)]
  return(x[ind])
}
