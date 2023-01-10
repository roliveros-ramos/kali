
#' @export
check_taxon = function(x, na.return=FALSE, verbose=TRUE) {
  
  if(!requireNamespace("taxize", quietly = TRUE)) 
    stop("You need to install the 'taxize' package.")
  
  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")
  
  .checkScientificName = function(x, na.return, verbose) {
    if(is.na(x)) return(NA)
    tmp = taxize::gnr_resolve(sci = x, canonical = TRUE)$matched_name2
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
    if(!id & !isNA & isTRUE(verbose)) message(sprintf("Taxon '%s' was corrected to '%s'.", x, tmp))
    return(tmp)
  }
  
  sx = na.omit(unique(x))
  
  out = unlist(sapply(sx, .checkScientificName, 
                      na.return=na.return, verbose=verbose))
  
  out = out[x]
  names(out) = NULL
  
  if(is.null(out)) {
    # case?
  }
  
  return(out)
}

#' @export
get_taxon = function(x, rank, db="itis") {
  
  if(!requireNamespace("taxize", quietly = TRUE)) 
    stop("You need to install the 'taxize' package.")
  
  if(length(rank)!=1) stop("Only one taxon is allowed.")
  
  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")
  
  db = match.arg(db, c("itis", "ncbi", "both"))
  isna = is.na(x)
  xna = unique(x[!isna])
  tmp = taxize::tax_name(sci=xna, db=db, get=rank, messages=FALSE, ask=FALSE)[, rank]
  xind = which(is.na(tmp) & .is_binomial(xna))
  if(length(xind)>0) {
    last_x = .get_first(xna[xind])
    tmp[xind] = taxize::tax_name(query=last_x, db=db, get=rank, messages=FALSE, ask=FALSE)[, rank]    
  }
  names(tmp) = xna
  out = tmp[x]
  names(out) = NULL
  return(out)
  
}

#' @export
get_classification = function(x) {
  x = check_taxon(x, verbose=FALSE)
  taxa = c('kingdom', 'subkingdom', 'infrakingdom', 'phylum', 'subphylum', 
           'infraphylum', 'superclass', 'class', 'superorder', 'order', 'suborder', 
           'family', 'subfamily', 'genus', 'species')
  out = setNames(rep(NA, length(taxa)), nm=taxa)
  x = taxize::classification(x, db="itis", messages=FALSE, ask=FALSE)[[1]]
  out[x$rank] = x$name
  out = as.list(out)
  return(out)
}

# Internal ----------------------------------------------------------------

.is_binomial = function(x) {
  out = sapply(strsplit(x, split=" "), length) == 2
  return(out)
}

.get_first = function(x) {
  out = sapply(strsplit(x, split=" "), "[", i=1)
  return(out)
}

