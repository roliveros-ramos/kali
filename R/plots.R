# plot.niche.models alias for image.niche.model

# image.map ---------------------------------------------------------------

#' @title Nice image plots for maps
#' @description This function plots an image for a georeferenced matrix, 
#' assuming x and y axes are longitude and latitude.  
#' @param lon vector or matrix of longitude values, following the conventions
#' of \code{image}.
#' @param lat vector or matrix of latitude values, following the conventions
#' of \code{image}.
#' @param z matrix of values, following the conventions of \code{image}.
#' @param \dots Additional parameters to be passed to \code{image.plot}.
#' @author Ricardo Oliveros-Ramos, modified from xxx's \code{image.plot} from 
#' the fields package
#' @examples
#' load(anchovy)
#' image.map(lon=anchovy$lon, lat=anchovy$lat, z=anchovy$z)
#' @export image.map
image.map = function (lon, lat, z, center=0, legend=TRUE, hires=FALSE, add = FALSE, nlevel = 1000, horizontal = FALSE, 
                      legend.shrink = 0.9, legend.width = 1.2, slice=NULL,
                      legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
                      bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
                      col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)), 
                      lab.breaks = NULL, axis.args = NULL, legend.args = NULL, axes=TRUE,
                      midpoint = FALSE, border = TRUE, lwd = 1, land.col="black", labels = TRUE,
                      sea.col="white", boundaries.col="grey", grid=FALSE, grid.col="white", ...) {
  
  if(!is.null(attr(z, "longitude"))) lon = attr(z, "lon")
  if(!is.null(attr(z, "latitude"))) lat = attr(z, "lat")
    
  if(length(dim(z))==3) {
    if(!is.null(slice)) z = z[, , slice]
    if(is.null(slice)) stop("Must specify an slice for plotting arrays.")
  }
  
  if(length(dim(z))>3) stop("Only arrays of 3-dimensions are supported.")
  
  pm = .findPrimeMeridian(lon)
  
  if(!isTRUE(legend)) {
    .image.mapnl(lon=lon, lat=lat, z=z, center=center, hires=hires, add=add, nlevel=nlevel, 
                 col=col, land.col=land.col, sea.col=sea.col, boundaries.col=boundaries.col, 
                 grid.col=grid.col, grid=grid, axes=axes, border=border, ...)
    return(invisible())
  }
  
  old.par = par(no.readonly = TRUE)
  info = imageplot.info(x=lon, y=lat, z=z, ...)
  if (add) {
    big.plot = old.par$plt
  }
  if (legend.only) {
    graphics.reset = TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar = ifelse(horizontal, 3.1, 5.1)
  }
  temp = imageplot.setup(add = add, legend.shrink = legend.shrink, 
                         legend.width = legend.width, legend.mar = legend.mar, 
                         horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot = temp$smallplot
  bigplot = temp$bigplot
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(x=lon, y=lat, z=z, add = add, col = col, axes=FALSE, 
            xlab="", ylab="", ...)
      mapDetails(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE, 
                 axes=axes, border=border, boundaries.col=boundaries.col,
                 grid=grid, grid.col=grid.col, water=sea.col, labels=labels)    
    }
    else {
      poly.image(x=lon, y=lat, z=z, add = add, col = col, midpoint = midpoint, 
                 border = border, lwd.poly = lwd, axes=FALSE, 
                 xlab="", ylab="",...)
      mapDetails(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE, 
                 axes=axes, border=border, boundaries.col=boundaries.col,
                 grid=grid, grid.col=grid.col, water=sea.col, labels=labels)  
    }
    big.par = par(no.readonly = TRUE)
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  ix = 1
  minz = info$zlim[1]
  maxz = info$zlim[2]
  binwidth = (maxz - minz)/nlevel
  midpoints = seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy = midpoints
  iz = matrix(iy, nrow = 1, ncol = length(iy))
  breaks = list(...)$breaks
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    axis.args = c(list(side = ifelse(horizontal, 1, 4), 
                       mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                       at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args = c(list(side = ifelse(horizontal, 1, 4), 
                       mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                  axis.args)
  }
  if (!horizontal) {
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    }
    else {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
  }
  else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    }
    else {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
  }
  do.call("axis", axis.args)
  box()
  if (!is.null(legend.lab)) {
    legend.args = list(text = legend.lab, side = ifelse(horizontal, 
                                                        1, 4), line = legend.mar - 2)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  mfg.save = par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    par(mar = par("mar"))
    invisible()
  }
  
  return(invisible())
}


# plot.map ----------------------------------------------------------------

#' @title Nice scatter plots for georeferenced data
#' @description This function makes and scatter plot given the longitude and latitude
#' of the data points.  
#'
#' @param x vector of longitude values, or a data.frame containing variables named 'lon'
#' and 'lat'
#' @param y vector of latitude values.
#' @param xlim 
#' @param ylim 
#' @param center 
#' @param hires 
#' @param land.col 
#' @param sea.col 
#' @param boundaries.col 
#' @param grid.col 
#' @param grid 
#' @param cex 
#' @param pch 
#' @param main 
#' @param add 
#' @param axes 
#' @param border 
#' @param asp 
#' @param axs 
#' @param xaxs 
#' @param yaxs 
#' @param cex.axis 
#' @param interior 
#' @param fill 
#' @param countries 
#' @param ... 
#' @param domain
#'
#' @author Ricardo Oliveros-Ramos
#' @examples
#' load(mackerel)
#' plot.map(mackerel)
#' plot.map(domain="peru")
#' @export plot.map
plot.map = function(x, y=NULL, xlim=NULL, ylim=NULL, domain=NULL, center=0, 
                    hires=FALSE, land.col="darkolivegreen4", sea.col="aliceblue", 
                    boundaries.col = "black", grid.col="white", grid=TRUE,
                    cex=0.5, pch=19, main=NULL, add=FALSE, axes=TRUE, 
                    border=!axes, asp=NA, axs="i", xaxs=axs, yaxs=axs, cex.axis=0.75, 
                    interior=FALSE, fill=TRUE, countries=FALSE, nx=NULL, ny=nx, ...) {
  
  if(missing(x)) x = NA
  
  if(is.data.frame(x) & is.null(y)) {
    names(x) = tolower(names(x))
    y = x$lat
    x = x$lon
  }
  
  if(is.character(x)) {
    domain = x
    x      = NA
  }
  
  xy = xy.coords(x, y)
  xy$xlab = ""
  xy$ylab = ""
  
  if(!is.null(domain)) {
    if(!is.null(xlim)) warning("Domain specified, ignoring 'xlim'.")
    xlim = .getDomain(domain, "x")
    if(!is.null(ylim)) warning("Domain specified, ignoring 'ylim'.")
    ylim = .getDomain(domain, "y")
  }
  
  if(is.null(xlim)) xlim = findXlim(xy$x)
  if(is.null(ylim)) ylim = pmax(-90, pmin(range(pretty(xy$y), na.rm=TRUE), 90))
  
  xlim = addPM(xlim)
  pm = attr(xlim, "pm")
  
  xy$x = checkLongitude(xy$x, primeMeridian = pm)

  if(!add) {
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs, asp=asp)
    .plotSea(col=sea.col)
    mapDetails(primeMeridian=pm, hires=hires, col=land.col, interior=interior, 
               axes=axes, border=border, boundaries.col=boundaries.col,
               grid=grid, grid.col=grid.col, cex.axis=cex.axis, fill=fill, 
               water=sea.col, countries=countries, nx=nx, ny=ny)    
    title(main=main)
  }
  points(xy, cex=cex, pch=pch, ...)
  
  return(invisible(pm))
}


# mapDetails --------------------------------------------------------------

mapDetails = function(primeMeridian="center", hires=FALSE, col="black", interior=FALSE, 
                      axes=TRUE, border=TRUE, boundaries.col="black",
                      grid=TRUE, grid.col="white", cex.axis=0.75, fill=TRUE, 
                      boundary = TRUE, water=NULL, countries=FALSE, nx, ny, labels=TRUE, ...) {
  
  primeMeridian = match.arg(primeMeridian, choices=c("center", "left"))
  
  if(hires) {
    if(!requireNamespace("mapdata", quietly = TRUE)) {
      warning("You need to install the 'mapdata' package, using hires=FALSE.")
      hires = FALSE
    }
  }
  
  mapa =  if(hires) {
    if(primeMeridian=="center") "mapdata::worldHires" else "mapdata::world2Hires"
  } else {
    if(primeMeridian=="center") "world" else "world2"
  }

  # mapa = if(primeMeridian=="center") "world" else "world2"
  
  wrap = ifelse(primeMeridian=="center", c(-180,180), c(0,360))
  
  if(isTRUE(grid)) grid(nx=nx, ny=ny, col=grid.col, lty=1)
  
  map(database=mapa, fill = fill, col = col, add = TRUE, interior=FALSE, 
      border=col, boundary = boundary, ...)
  
  if(!is.null(water)) {
    if(primeMeridian=="center") {
      map("lakes", fill = TRUE, col = water, add = TRUE, border=water, ...)
    } else {
      lakes = map("lakes", plot=FALSE)
      lakes$x = checkLongitude(lakes$x, primeMeridian = "left")
      map(lakes, fill = TRUE, col = water, add = TRUE, border=water, ...)
    }
  }
  
  if(isTRUE(countries)) 
    map(database=mapa, fill = FALSE, col = boundaries.col, add = TRUE)
  
  
  if(axes) {
    map.axes2(cex.axis=cex.axis)
    if(isTRUE(labels)) {
      mtext("LONGITUDE", 1, line = 1.8*cex.axis/0.75, cex = 0.9*par("cex"))
      mtext("LATITUDE", 2, line = 2.4*cex.axis/0.75, cex = 0.9*par("cex"))    
    }
  } else {
    if(border) box()
  }
  
  return(invisible())
}

#' @export plot.domain
plot.domain = function(domain=NA, xlim=NULL, ylim=NULL, col="red", lwd=1, fill=NA, ...) {
  
  if(is.null(xlim)) xlim = .getDomain(domain, "x")
  if(is.null(xlim)) stop("latitudinal domain limits has not been specified.")
  
  if(is.null(ylim)) ylim = .getDomain(domain, "y")
  if(is.null(ylim)) stop("longitudinal domain limits has not been specified.")
  
  xpol = c(xlim, rev(xlim))
  ypol = c(ylim[1],ylim[1],ylim[2],ylim[2])
  polygon(x=xpol, y=ypol, border=col, lwd=lwd, col=fill, ...)
  
  return(invisible())
}


countMap = function(data, var, lat, lon, dx=1, dy=dx, ...) {
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  
  latAsFactor = cut(data[,"lat"], coords$psi$lat, labels=FALSE)
  lonAsFactor = cut(data[,"lon"], coords$psi$lon, labels=FALSE)
  
  map = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor),
               FUN=length)
  
  rows = seq_along(coords$rho$lat)
  cols = seq_along(coords$rho$lon)
  
  complete.rows = as.numeric(rownames(map))
  complete.cols = as.numeric(colnames(map))
  
  missing.rows = rows[!(rows %in% complete.rows)]
  missing.cols = cols[!(cols %in% complete.cols)]
  
  map = cbind(map, matrix(ncol=length(missing.cols), nrow=nrow(map)))
  map = rbind(map, matrix(ncol=ncol(map), nrow=length(missing.rows)))
  
  irows = sort(c(complete.rows, missing.rows), index.return=TRUE)$ix
  icols = sort(c(complete.cols, missing.cols), index.return=TRUE)$ix
  
  map = map[irows,icols]
  map = map[nrow(map):1, ] # correct lat direction
  map = .rotate(map) # write map for ncdf file, origin in down left corner
  
  rownames(map) = coords$rho$lon
  colnames(map) = coords$rho$lat
  
  image.map(coords$lon, coords$lat, map, ...)
  
  return(map)
}

barplot2 = function(hist, horiz=TRUE, ...) {
  
  xax = if(isTRUE(horiz)) 2 else 1
  yax = if(isTRUE(horiz)) 1 else 2
  
  hist2 = data.frame(counts=hist$counts,
                     mids=hist$mids,
                     density=hist$density)
  
  hist2$bp = barplot(hist$counts, horiz=horiz, axes=FALSE, ...)
  
  axis(yax)
  md = lm(bp ~ mids, data=hist2)
  posx = predict(md, newdata=data.frame(mids=hist$breaks))
  axis(xax, at=posx, labels=coord2text(hist$breaks, type="lat"),
       las=xax)
  box()
  
  hist$axis = posx
  
  return(invisible())
}



# keep checking

plotMAPt = function(lon, lat, map, thr, ...) {
  
  image.map(lon=lon, lat=lat, 0 + map>thr, ...)
}

printMaps = function(object, dir, prefix, width=800, height=800, res=NA, ...) {
  
  filename = file.path(dir, prefix)
  
  plot(object, png=TRUE, width=width, height=height, res=res, filename=filename, ...)
  
  plot(object, type="seasonal", png=TRUE, width=width, height=height, res=res, filename=filename, ...)
  
  plot(object, type="climatology", png=TRUE, width=2*width, height=1.5*height, 
       res=res, filename=filename, ...)
  
  return(invisible())
}

