
calculateNiche = function(model, nmax=1e5, doIt=FALSE, ...) {
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
  output = list(data=newdata, var=values)
  class(output) = c("niche")
  return(output)  
}

plot.niche = function(x, vars, FUN=median, plot=TRUE, n=200, thr=NULL, 
                      type="prob", add=FALSE, col="black", lwd=2, bezier=TRUE,...) {
  out = switch(type, 
               prob = .plotNicheProb(x=x, vars=vars, FUN=FUN, plot=plot, n=n, thr=thr, ...),
               hull = .plotNicheHull(x=x, vars=vars, FUN=FUN, plot=plot, n=n, thr=thr, 
                                     add=add, col=col, lwd=lwd, bezier=bezier, ...),
               stop("Invalid plot type."))
  return(invisible(out))
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
    image.plot(out, xlab=out$labs$x, ylab=out$labs$y, zlim=c(0,1))    
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
  out$x = approx(x=seq(0, 1, len=length(x$x)), y=x$x, xout=xout)$y
  out$y = approx(x=seq(0, 1, len=length(x$y)), y=x$y, xout=xout)$y
  out = bezier::bezier(seq(0, 1, len=5*n), out)
  return(out)
}

.plotNicheHull = function(x, vars, FUN=median, plot=TRUE, n=200, 
                          thr=NULL, add=FALSE, col="black", lwd=2, bezier=TRUE, ...) {
  
  if(is.null(thr)) thr=0.5
  
  old = list()
  old$x = x$var[[vars[1]]]
  old$y = x$var[[vars[2]]]
  old$z = tapply(X=x$data$niche, INDEX=x$data[vars], FUN=FUN)
  old$labs$x = vars[1]
  old$labs$y = vars[2]
  
  nold = mean(length(old$x), length(old$y))
  
  if(nold<n) {
    new = list()
    new$x = pretty(old$x, n=n)
    new$y = pretty(old$y, n=n)
    new$z = interp.surface.grid(obj=old, grid.list=new)$z
    new$labs = old$labs    
    out = new
  } else out = old
  
  out$z = toPA(out$z, thr=thr)
  
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
