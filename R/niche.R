pretty.factor = function(x, ...) return(factor(levels(x)))

calculateNiche = function(model, nmax=1e5, doIt=FALSE, req.sens=0.95, cost=list(FPC=1, FNC=10), ...) {
  nvars = length(model$var.summary)
  if(nvars>=7 & !doIt) 
    stop("Calculating the niche for more than 6 variables
         may take a long time, set doIt=TRUE to proceed.")
  if(nvars>=7 & !doIt) DateStamp("Starting...")
  n = max(min(100, floor((nmax)^(1/nvars))), 10)
  values = lapply(model$var.summary, pretty, n=n)
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  if(nvars>=7 & !doIt) DateStamp("Ending at")
  fmla = .fmla2txt(model$formula)
  species = as.character(model$formula[2])
  model$model$fitted = predict(model, newdata = model$model, type="response")
  thr = calculateThresholds(data=model$model, coordNames=names(model$var.summary),
                            models="fitted", obs=species, req.sens=req.sens, FPC=cost$FPC, FNC=cost$FPC)
  performance = PredictivePerformance(data=model$model, coordNames=names(model$var.summary),
                                      models="fitted", obs=species, st.dev=FALSE)
  output = list(data=newdata, var=values, model=model$model, species=species, formula=fmla, thr=thr,
                performance=performance)
  class(output) = c("niche")
  return(output)  
}

plot.niche = function(x, vars, FUN=median, plot=TRUE, n=200, thr=NULL, 
                      type="prob", add=FALSE, col="black", lwd=2, bezier=TRUE, ...) {
  out = switch(type, 
               prob = .plotNicheProb(x=x, vars=vars, FUN=FUN, plot=plot, n=n, thr=thr, ...),
               hull = .plotNicheHull(x=x, vars=vars, FUN=FUN, plot=plot, n=n, thr=thr, 
                                     add=add, col=col, lwd=lwd, bezier=bezier, ...),
               stop("Invalid plot type."))
  return(invisible(out))
}

points.niche = function(x, vars, pch=".", col="black", alpha=0.9, n=1, ...) {
  if(n<0 | n>1) stop("n must be in [0,1].")
  dat = x$model[,vars]
  indN = which(x$model[, x$species] == 0)
  indP = which(x$model[, x$species] == 1)
  col = makeTransparent(col, alpha=alpha)
  ind = sample(indP, size=floor(n*length(indP)), replace=FALSE)
  points(dat[ind, ], pch=pch, col=col, ...)
  return(invisible())
}

.plotNicheProb = function(x, vars, FUN=median, plot=TRUE, n=200, thr=NULL, ...) {
  
  old = list()
  old$x = x$var[[vars[1]]]
  old$y = x$var[[vars[2]]]
  old$z = tapply(X=x$data$niche, INDEX=x$data[vars], FUN=FUN, ...)
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
  
  if(!is.null(thr)) out$z = toPA(out$z, thr=thr)
  
  if(isTRUE(plot)) {
    image(out, xlab=out$labs$x, ylab=out$labs$y, zlim=c(0,1), col=tim.colors(64))    
  }
  return(invisible(out))
}

.getConvexHull = function(x) {
  hpts = chull(cbind(x$x,x$y))
  hpts = c(hpts, hpts[1])
  out = list(x=x$x[hpts], y=x$y[hpts])
  return(out)
}

.doSmooth = function(x, n) {
  out = list()
  xout = seq(0, 1, len=0.25*n)
  if(length(x$x)<4) return(list(x=NULL, y=NULL))
  out$x = approx(x=seq(0, 1, len=length(x$x)), y=x$x, xout=xout)$y
  out$y = approx(x=seq(0, 1, len=length(x$y)), y=x$y, xout=xout)$y
  out = bezier::bezier(seq(0, 1, len=5*n), out)
  return(out)
}

.plotNicheHull = function(x, vars, FUN=median, plot=TRUE, n=200, 
                          thr=NULL, add=FALSE, col="black", lwd=2, bezier=TRUE, ...) {
  
  if(is.null(thr)) thr=x$thr["ReqSens",]
  
  out = .plotNicheProb(x=x, vars=vars, FUN=FUN, plot=FALSE, n=n, thr=thr)
  
  coords = list()
  inside = which(out$z==1)
  coords$x = matrix(out$x, nrow=nrow(out$z), ncol=ncol(out$z))[inside]
  coords$y = matrix(out$y, nrow=nrow(out$z), ncol=ncol(out$z), byrow=TRUE)[inside]
  
  hull = .getConvexHull(coords)
  if(isTRUE(bezier)) hull = .doSmooth(hull, n=n)
  
  if(isTRUE(plot)) {
    if(!isTRUE(add)) {
      plot.new()
      plot.window(xlim=range(out$x), ylim=range(out$y))
      title(xlab=out$labs$x, ylab=out$labs$y)
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
