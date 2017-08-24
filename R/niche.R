pretty.factor = function(x, ...) return(factor(levels(x)))

xrange = function(x, n) {
  UseMethod("xrange")
}

xrange.default = function(x, n) {
  xr = range(pretty(x, n=100), na.rm=TRUE)
  out  = seq(from=xr[1], to=xr[2], length.out=n)
  return(out)
} 

xrange.factor = function(x, n) return(factor(levels(droplevels(x))))
quantile.factor = function(x, n) return(factor(levels(droplevels(x))))

calculateNiche = function(object, model=NULL, nmax=1e6, doIt=FALSE, 
                          alpha=0.99, cluster=NULL, FUN=mean, 
                          verbose=FALSE, ...) {
    
  on.exit(invisible(gc()))
  
  if(!inherits(object, "niche.models")) 
    stop("object must be of class 'niche.models'")

  delta  = (1-alpha)/2 # for environmental range calculation
  
  modelNames = model
  if(is.null(modelNames)) modelNames = names(object$model)
  
  indM = modelNames %in% names(object$model)
  modelNames = modelNames[indM]
  if(length(modelNames)==0) stop("No models match model fitted, check 'model' names.")
  
  if(sum(!indM)>0) {
    discard = paste(modelNames[!indM], collapse=", ")
    msg = if(sum(!indM)==1) sprintf("One model have been discarded: %s.", discard) else 
      sprintf("Some (%d) models have been discarded: %s.", sum(!indM), discard)
    warning(msg)
  }  
  
  models = object$model[modelNames]
    
  usedVars = lapply(object$formulas[modelNames], FUN=.getUsedVars, x=names(object$train))
  usedVars = unique(unlist(usedVars))
  niche2Dvars  = combn(x = usedVars, 2)
  
  nvars = length(usedVars)
  # nvars = length(model$var.summary)
  
  if(nvars>=7 & !doIt) 
    stop("Calculating the niche for more than 6 variables
         may take a long time, set doIt=TRUE to proceed.")
  DateStamp("Starting...")

  cat(sprintf("Calculating niche for %s.\n", paste(usedVars, collapse=", ")))
  
  nmax = nvars*nmax
  n  = floor(max(min(300, floor((nmax)^(1/nvars))), 5))

  values = lapply(object$train[, usedVars], xrange, n=n)
  
  for(iModel in modelNames) values[[iModel]] = NA_real_

  newdata = do.call(expand.grid, values)
  insideNiche = matrix(FALSE, nrow=nrow(newdata), ncol=length(modelNames))
  colnames(insideNiche) = modelNames

  ranges = list()
  niches2D = list()
  
  pb = txtProgressBar(style=3)
  setTxtProgressBar(pb, 0)
  
  for(iModel in modelNames) {
    
    newdata[[iModel]] = predict(models[[iModel]], newdata = newdata, 
                                   type="response", cluster=cluster)

    # after model predictions
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (2*which(modelNames==iModel)-1)/(2*length(modelNames)+1))
    
    iThr = object$thr$validation["ReqSens", iModel]
    if(isTRUE(verbose)) cat("\nUsing iThr=", iThr, "for model", iModel, ".\n")

    insideNiche[, iModel] = (newdata[[iModel]] > iThr)
    
    # check ranges
    ranges[[iModel]] = lapply(newdata[insideNiche[, iModel], usedVars], 
                              quantile, probs=c(delta, 1-delta))
   
    insideNiche[, iModel] = insideTesseract(x = newdata[, usedVars], 
                                            index = insideNiche[, iModel])
     
    niches2D[[iModel]] = calculateNiche2D(data=newdata, model=iModel, 
                                          index=insideNiche[, iModel], 
                                          usedVars=usedVars, FUN=FUN)
    
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (2*which(modelNames==iModel))/(2*length(modelNames)+1))
    
  }

  DateStamp("\nEnding at")
  
  species = as.character(object$formulas[[1]][2]) #check for species in fitGAMs
  train = object$train[, c(species, usedVars)]
   
  info = list(env=usedVars, model=modelNames, pairs=niche2Dvars)
  
  pb = txtProgressBar(style=3)
  setTxtProgressBar(pb, 1)
  
  output = list(data=newdata, model=train, var=values, species=species, 
                formula=object$formulas, thr=object$thr, tolerance=ranges,
                niches2D=niches2D, info=info)
  
  class(output) = c("niche")
  return(output)  
  
}

# data = newdata # data frame with all the data
# model # the name of the model
# index # points inside the tesseract
# pair # combination of variables
# usedVars # environmental variables used
calculateNiche2D = function(data, model, index, usedVars, FUN=FUN) {
  
  niche2Dvars  = combn(x = usedVars, 2)
  out = list()
  for(i in seq_len(ncol(niche2Dvars))) {
    out[[i]] = .calculateNiche2D(data=data, model=model, index=index, 
                                 pair=niche2Dvars[,i], FUN=FUN)
  }
  return(out)
}

.calculateNiche2D = function(data, model, index, pair, FUN=mean) {
  
  nicheData = data[[model]]
  envData = data[pair]
  xranges = lapply(envData, FUN = function(x) sort(unique(x)))
  
  old = list()
  old$x = xranges[[1]]
  old$y = xranges[[2]]
  z1 = tapply(X=nicheData, INDEX=envData, FUN=FUN, na.rm=TRUE)
  nicheData[!index] = NA
  z2 = tapply(X=nicheData, INDEX=envData, FUN=FUN, na.rm=TRUE)
  z2[is.na(z2)] = z1[is.na(z2)] 
  old$z = z2 
  old$labs$x = pair[1]
  old$labs$y = pair[2]
  
  return(old)
}


# plot for niche class ----------------------------------------------------


plot.niche = function(x, vars, what, plot=TRUE, n=1, thr=NULL, 
                      criteria="MinROCdist", type="prob", 
                      col=tim.colors(64), hull=FALSE, 
                      col.hull="black", lwd=2, ...) {
  
  if(!(what %in% x$info$model)) stop("Model names ('what') not found.")
  if(length(vars)!=2) stop("The argument 'vars' must have length 2.")
  
  npair = .getPair(x$info$pairs, vars)
  reverse = attr(npair, "reverse")
  
  xthr = x$thr$validation[, what, drop=FALSE]
  
  x = x$niches2D[[what]][[npair]]
  if(isTRUE(reverse)) {
    x0 = x$x; y0 = x$y
    x$x = y0; x$y = x0
    x$z = t(x$z)
  }
 
  x$labs = list(x=vars[1], y=vars[2]) 
  x$thr = xthr
  
  out = switch(type, 
               prob = .plotNicheProb2(x=x, vars=vars, plot=plot, n=n, thr=thr, 
                                      criteria=criteria, col=col,
                                      hull=hull, col.hull=col.hull, lwd=lwd, ...),
               hull = .plotNicheHull(x=x, vars=vars, plot=plot, n=n, thr=thr, 
                                     criteria=criteria, col=col.hull, 
                                     lwd=lwd, ...),
               pa   = .plotNicheProb2(x=x, vars=vars, plot=plot, n=n, thr=thr, 
                                      criteria=criteria, col=col,
                                      hull=hull, col.hull=col.hull, lwd=lwd, toPA=TRUE, ...),
               stop("Invalid plot type."))
  return(invisible(out))
  
}

points.niche = function(x, vars, pch=".", col=c("black", "grey"), alpha=0.9, n=1, 
                        what="presence", ...) {

  what = match.arg(arg=what, choices=c("presence", "absence", "both"))

  col = makeTransparent(col, alpha=alpha)
  
  if(what=="both" & length(col)!=2) stop("Argument 'col' has to be length 2.") 
  usePositive = if(what=="presence") TRUE else FALSE
  useNegative = if(what=="absence") TRUE else FALSE
  if(what=="both") {
    usePositive = TRUE 
    useNegative = TRUE
    colP = col[1]
    colN = col[2]
  } else colP = colN = col[1]
  
  if(n<0 | n>1) stop("n must be in [0,1].")
  dat = x$model[,vars]
  
  indN = which(x$model[, x$species] == 0)
  indN = sample(indN, size=floor(n*length(indN)), replace=FALSE)
  
  indP = which(x$model[, x$species] == 1)
  indP = sample(indP, size=floor(n*length(indP)), replace=FALSE)
  
  if(useNegative) points(dat[indN, ], pch=pch, col=colN, ...)
  if(usePositive) points(dat[indP, ], pch=pch, col=colP, ...)
  
  return(invisible())
}

.plotNicheProb2 = function(x, vars, FUN=mean, plot=TRUE, n=200, thr=NULL, 
                          criteria="MinROCdist", toPA=FALSE,
                          col=tim.colors(64), hull=FALSE, col.hull="black", 
                          lwd=2, ...) {
  
  if(is.null(thr)) thr = x$thr[criteria, ]
  
  old = x[c("x", "y", "z", "labs")]
    
  nold = mean(length(old$x), length(old$y))
  
  if(nold<n) {
    new = list()
    new$x = xrange(old$x, n=n)
    new$y = xrange(old$y, n=n)
    new$z = interp.surface.grid(obj=old, grid.list=new, ...)$z
    new$labs = old$labs    
    out = new
  } else out = old
 
  if(isTRUE(toPA)) out$z = toPA.default(x=out$z, thr=thr)
 
  if(isTRUE(hull)) {
    xout = out
    xout$z = toPA.default(x=xout$z, thr=thr)
    out$hull = .calculateNicheHull(out=xout)
  } 
  
  if(isTRUE(plot)) {
    image(out, xlab=out$labs$x, ylab=out$labs$y, zlim=c(0,1), col=col)    
    if(isTRUE(hull)) lines(out$hull, col=col.hull, lwd=lwd, ...)
  }
  return(invisible(out))
}

.plotNicheProb3 = function(x, vars, FUN=mean, plot=TRUE, n=200, thr=NULL, 
                           criteria="MinROCdist", toPA=FALSE,
                           col=tim.colors(64), hull=FALSE, col.hull="black", 
                           lwd=2, ...) {
  
  if(is.null(thr)) thr = x$thr[criteria, ]
  mthr = x$thr["ReqSens",]
  
  old = list()
  old$x = x$var[[vars[1]]]
  old$y = x$var[[vars[2]]]
  # x$data$niche = toPA.default(x=x$data$niche, thr=thr)
  x$data$niche[x$data$niche<mthr] = NA
  old$z = tapply(X=x$data$niche, INDEX=x$data[vars], FUN=FUN, na.rm=TRUE)
  old$z[is.na(old$z)] = mthr
  old$labs$x = vars[1]
  old$labs$y = vars[2]
  
  nold = mean(length(old$x), length(old$y))
  
  if(nold<n) {
    new = list()
    new$x = pretty(old$x, n=n)
    new$y = pretty(old$y, n=n)
    new$z = interp.surface.grid(obj=old, grid.list=new, ...)$z
    new$labs = old$labs    
    out = new
  } else out = old
  
  # out$z[out$z<=mthr] = 0
  
  if(isTRUE(toPA)) out$z = toPA.default(x=out$z, thr=thr)
  
  if(isTRUE(hull)) {
    xout = out
    xout$z = toPA.default(x=xout$z, thr=thr)
    out$hull = .calculateNicheHull(out=xout)
  } 
  
  if(isTRUE(plot)) {
    image(out, xlab=out$labs$x, ylab=out$labs$y, zlim=c(0,1), col=col)    
    if(isTRUE(hull)) lines(out$hull, col=col.hull, lwd=lwd, ...)
  }
  return(invisible(out))
}

.calculateNicheHull <- function(out) {
  
  coords = list()
  inside = which(out$z==1)
  coords$x = matrix(out$x, nrow=nrow(out$z), ncol=ncol(out$z))[inside]
  coords$y = matrix(out$y, nrow=nrow(out$z), ncol=ncol(out$z), byrow=TRUE)[inside]
  
  hull = .getConvexHull(coords)
  hull = .doSmooth(hull)
  
  return(hull)
  
}

.getConvexHull = function(x) {
  hpts = chull(cbind(x$x,x$y))
  hpts = c(hpts, hpts[1])
  out = list(x=x$x[hpts], y=x$y[hpts])
  return(out)
}

.doSmooth = function(x) {

  if(length(x$x)<4) return(list(x=NULL, y=NULL))
  out = list()
  xout1 = seq(0, 1, len=length(x$x))
  xout2 = seq(0, 1, len=30*length(x$x))
  xout = sort(unique(c(xout1, xout2)))
  out$x = approx(x=seq(0, 1, len=length(x$x)), y=x$x, xout=xout)$y
  out$y = approx(x=seq(0, 1, len=length(x$y)), y=x$y, xout=xout)$y
  out = bezier::bezier(seq(0, 1, len=length(out$x)), out)
  return(out)
}

.plotNicheHull = function(x, vars, FUN=median, plot=TRUE, n=200, criteria="MinROCdist",
                          thr=NULL, add=FALSE, col="black", lwd=2, ...) {
  
  ithr=x$thr[criteria,]
  
  out = .plotNicheProb2(x=x, vars=vars, FUN=FUN, plot=FALSE, n=n, thr=thr, 
                       criteria=criteria, toPA=TRUE)
 
   
  hull = .calculateNicheHull(out=out)
  
  if(isTRUE(plot)) {
    if(!isTRUE(add)) {
      plot.new()
      plot.window(xlim=range(out$x), ylim=range(out$y))
      # title(xlab=out$labs$x, ylab=out$labs$y)
      axis(1); axis(2); box()
    }
    lines(hull, col=col, lwd=lwd, ...)    
  }
  return(invisible(hull))
}

density.niche = function(x, var, plot=TRUE, vertical=FALSE,
                         axes=TRUE, col=c("black", "grey"),
                         lty=c(1,1), lwd=c(1,1), allData=TRUE,
                         ...) {
  dat = x$model[, var]
  indP = which(x$model[, species]==1)
  d1 = density(dat[indP])
  d2 = density(dat)
  if(isTRUE(plot)) {
    plot.new()
    if(!vertical) {
      plot.window(xlim=range(pretty(dat)), 
                  ylim=range(c(d1$y, d2$y)))
      lines(d1$x, d1$y, col=col[1], lty=lty[1], lwd=lwd[1])
      if(isTRUE(allData)) lines(d2$x, d2$y, col=col[2], lty=lty[2], lwd=lwd[2])
      
    } else {
      plot.window(xlim=range(c(d1$y, d2$y)), 
                  ylim=range(pretty(dat)))
      lines(d1$y, d1$x, col=col[1], lty=lty[1], lwd=lwd[1])
      if(isTRUE(allData)) lines(d2$y, d2$x, col=col[2], lty=lty[2], lwd=lwd[2])
    }
    if(isTRUE(axes)) {
      axis(1)
      axis(2)
      box()
    }
  }
  return(invisible(list(presence=d1, all=d2)))
}


# OLD ---------------------------------------------------------------------

calculateNiche2 = function(model, nmax=1e5, doIt=FALSE, req.sens=0.90, 
                          cost=list(FPC=1, FNC=10), alpha=0.95, ...) {
  nvars = length(model$var.summary)
  if(nvars>=7 & !doIt) 
    stop("Calculating the niche for more than 6 variables
         may take a long time, set doIt=TRUE to proceed.")
  if(nvars>=7 & !doIt) DateStamp("Starting...")
  
  n  = floor(max(min(100, floor((nmax)^(1/nvars))), 10))
  
  # first exploration to set limits
  values = lapply(model$var.summary, xrange, n=n/3)
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  
  fmla = .fmla2txt(model$formula)
  species = as.character(model$formula[2])
  model$model$fitted = predict(model, newdata = model$model, type="response")
  thr = try(calculateThresholds(data=model$model, coordNames=names(model$var.summary),
                                models="fitted", obs=species, req.sens=req.sens, FPC=cost$FPC, FNC=cost$FPC))
  
  # final calculation within limits
  values = lapply(newdata[newdata$niche>thr["ReqSpec", 1], ], xrange, n=n)
  delta  = (1-alpha)/2
  ranges = lapply(newdata[newdata$niche>thr["ReqSpec", 1], ], quantile, 
                  probs=c(delta, 1-delta))
  
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  
  if(nvars>=7 & !doIt) DateStamp("Ending at")
  
  performance = try(PredictivePerformance(data=model$model, coordNames=names(model$var.summary),
                                          models="fitted", obs=species, st.dev=FALSE))
  output = list(data=newdata, var=values, model=model$model, species=species, formula=fmla, thr=thr,
                performance=performance, tolerance=ranges)
  class(output) = c("niche")
  return(output)  
}


calculateNiche3 = function(object, model=NULL, nmax=1e5, doIt=FALSE, req.sens=0.90, 
                          cost=list(FPC=1, FNC=10), alpha=0.95, ...) {
  
  if(!inherits(object, "niche.models")) 
    stop("object must be of class 'niche.models'")
  
  if(is.null(model)) model =  .getBestModel(object)
  
  model = object$model[[model]]
  
  nvars = length(model$var.summary)
  if(nvars>=7 & !doIt) 
    stop("Calculating the niche for more than 6 variables
         may take a long time, set doIt=TRUE to proceed.")
  if(nvars>=7 & !doIt) DateStamp("Starting...")
  
  n  = floor(max(min(100, floor((nmax)^(1/nvars))), 10))
  
  # first exploration to set limits
  values = lapply(model$var.summary, xrange, n=n/3)
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  
  fmla = .fmla2txt(model$formula)
  species = as.character(model$formula[2])
  model$model$fitted = model$fitted.values
  thr = try(calculateThresholds(data=model$model, coordNames=names(model$var.summary),
                                models="fitted", obs=species, req.sens=req.sens, FPC=cost$FPC, FNC=cost$FPC))
  
  # final calculation within limits
  values = lapply(newdata[newdata$niche>thr["ReqSpec", 1], ], xrange, n=n)
  
  # check range calculation
  delta  = (1-alpha)/2
  ranges = lapply(newdata[newdata$niche>thr["ReqSpec", 1], ], quantile, 
                  probs=c(delta, 1-delta))
  
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  
  if(nvars>=7 & !doIt) DateStamp("Ending at")
  
  performance = try(PredictivePerformance(data=model$model, coordNames=names(model$var.summary),
                                          models="fitted", obs=species, st.dev=FALSE))
  
  missVar = names(model$var.summary)[!(names(model$var.summary) 
                                       %in% names(model$model))]
  
  mVars = if(!is.null(model$na.action)) {
    object$train[-model$na.action, missVar]
  } else object$train[, missVar]
  
  model$model = cbind(model$model, mVars)
  
  output = list(data=newdata, var=values, model=model$model, species=species, formula=fmla, thr=thr,
                performance=performance, tolerance=ranges)
  
  class(output) = c("niche")
  return(output)  
}



# Auxiliar ----------------------------------------------------------------

.getPair = function(x, pair) {
  reverse = FALSE
  out = which(apply(x, 2, identical, y=pair))
  if(length(out)==0) {
    reverse = TRUE
    out = which(apply(x, 2, identical, y=rev(pair)))
  }
  if(length(out)==0) stop("The pair doesn't match with any combination.")
  attr(x=out, which="reverse") = reverse
  return(out)
}
