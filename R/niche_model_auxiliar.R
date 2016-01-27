# convert matrix or arrays to presence/absence data
toPA= function(x, thr=0, prob=FALSE) {
  if(!isTRUE(prob)) {
    if(all(thr==0)) return(.toPA(x))
    thr = sort(unique(c(-Inf, thr, 1)))
    probs = seq(from=0, to=1, len=length(thr)-1)
    xcut = cut(x, breaks=thr, labels=probs)
    x[] = as.numeric(as.character(xcut))
    return(x)
  } else {
    return(.toPA(x)*x)
  }
}

# to be deprecated
# toPAm = function(x, thr=c(0,0.5)) toPA(x, thr=thr) 

# remove absences using teseract exclusion
clearAbsences = function(data, control=NULL, species, vars, alpha=0.05, lowerOnly=FALSE) {
  
  dataN    = data[data[,species]==0,]
  dataPx   = data[data[, species]==1, ]
  dataP    = if(is.null(control)) dataPx else control[control[, species]==1, ]
  
  dataN_PA = cleanZeros(data=dataN, control=dataP, vars=vars, 
                        alpha=alpha, lowerOnly=lowerOnly)
  data_PA  = rbind(dataPx, dataN_PA)
  
  return(data_PA)
}

# remove zeros using teseract exclusion
cleanZeros = function(data, control, vars, alpha=0.05, lowerOnly=FALSE) {
  
  if(length(alpha)==1) alpha = rep(alpha, len=length(vars))
  if(length(alpha)!=length(vars)) stop("You must provide as many confidence thresholds as control vars.")
  if(length(lowerOnly)==1) lowerOnly = rep(lowerOnly, len=length(vars))
  if(length(lowerOnly)!=length(vars)) stop("You must provide as many confidence thresholds type as control vars.")
  ind = rep(TRUE, len=nrow(data))
  for(i in seq_along(vars)) {
    ind = ind & .insider(data, control, vars[i], alpha[i], lowerOnly[i])
  }
  ind = !ind & !is.na(ind)
  output = data[ind, ]
  
  return(output)
}

# balance prevalence
balancePA = function(data, PA, sp, rpa=0.1, seed=771104) {
  
  if(rpa<0 | rpa>1) error("ratio must be in [0,1].")
  ind.pa = which(PA[,sp]==0 & !is.na(PA[,sp]))
  ind.a  = which(data[,sp]==0 & !is.na(data[,sp]))
  ind.p  = which(data[,sp]==1 & !is.na(data[,sp]))
  
  np  = length(ind.p)
  nn  = length(ind.a)
  npa = max(np - nn, 0)
  ppa = npa/np
  
  if(ppa >= rpa) {
    set.seed(seed)
    ind.pa = sample(ind.pa, npa)
  } else {
    npa = floor(rpa*np)
    set.seed(seed)
    ind.pa = sample(ind.pa, npa)
    ind.a  = sample(ind.a, np-npa)
  }
  cat("Balancing prevalence in data for", sQuote(sp), ":\n")
  cat("Pseudo-Absence ratio = ", round(npa/np,1), " (minimum ratio= ", rpa,")\n", sep="")
  print( c(Presences=np, Absences=np-npa, 'Pseudo-Absences'=npa) )
  output = rbind(data[ind.p,], 
                 data[ind.a,],
                 PA[ind.pa,])
  
  cat("\nFinal data summary:\n")
  print(table.summary(output, sp=sp))
  cat("\n")
  
  return(output)
}

# table.summary (rename!!!)

