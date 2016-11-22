# General use functions ---------------------------------------------------

# shortcut to rename objects (currently)
#' @export
rename = function(x, old, new, ...) {
  UseMethod("rename")
}

# data.frame/grid transformations -----------------------------------------

# the next three functions have to be documented in the same help
# page (map2data, table2grid, table2array), think renaming
# map (grid) to data.frame
map2data = function(z, x, y) {
  
  # Transform a matrix-like object z to a data.frame
  # object
  # z is a matrix of dimensions length(x) times
  # length(y) 
  
  if(nrow(z) !=length(x) || ncol(z)!=length(y)) {
    stop("dimension of matrix z must be length(x) times
         # length(y)")
  }
  
  X = matrix(x, ncol=length(y), nrow=length(x))
  Y = matrix(y, ncol=length(y), nrow=length(x), byrow=T)
  
  X = .array2vector(X)
  Y = .array2vector(Y)
  Z = .array2vector(Z)
  
  out = data.frame(x=X, y=Y, z=z)
  
  return(out)
}

 
# convert a data.frame to matrix (grid)
table2grid = function(data, var, lat, lon, dx=dy, dy=dx, FUN=sum, ...) {
  
  FUN = match.fun(FUN)
  
  if(length(lat)< 2 | length(lon) < 2) stop("At lest a pair of lat and lon values must be provided.")
  
  latL = length(lat) 
  lonL = length(lon)
  
  if(latL == 2) {  
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
    latCut = coords$psi$lat
    cnames = coords$rho$lat
  } else {
    latCut = sort(lat)
    cnames = latCut[-1] - diff(latCut)
  }
  
  if(lonL == 2) {  
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
    lonCut = coords$psi$lon
    rnames = coords$rho$lon
  } else {
    lonCut = sort(lon)
    rnames = lonCut[-1] - diff(lonCut)
  }
  
  latAsFactor = cut(data[,"lat"], latCut, labels=FALSE)
  lonAsFactor = cut(data[,"lon"], lonCut, labels=FALSE)
  
  map = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor),
               FUN=FUN, ...)
  
  rows = seq_along(cnames)
  cols = seq_along(rnames)
  
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
  
  rownames(map) = rnames
  colnames(map) = cnames
  
  return(map)
}

# convert a data.frame to an array
table2array = function(data, var, lat, lon, start, end, 
                       dx, dy=dx, frequency=12, FUN=sum, toPA=FALSE) {
  
  FUN = match.fun(FUN)
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  latAsFactor  = cut(data[,"lat"], coords$psi$lat, labels=FALSE)
  lonAsFactor  = cut(data[,"lon"], coords$psi$lon, labels=FALSE)
  timeAsFactor = cut(data[,"time"], times$bounds, labels=FALSE)
  
  map0 = tapply(data[,var], INDEX=list(latAsFactor, lonAsFactor, timeAsFactor),
                FUN=FUN, na.rm=TRUE)
  
  rows    = seq_along(coords$rho$lat)
  cols    = seq_along(coords$rho$lon)
  slices  = seq_along(times$center)
  
  complete.rows   = as.numeric(rownames(map0))
  complete.cols   = as.numeric(colnames(map0))
  complete.slices = as.numeric(dimnames(map0)[[3]])
  
  missing.rows   = rows[!(rows %in% complete.rows)]
  missing.cols   = cols[!(cols %in% complete.cols)]
  missing.slices = slices[!(slices %in% complete.slices)]
  
  map = array(dim=c(length(rows), length(cols), length(slices)))
  map[seq_along(complete.rows), seq_along(complete.cols), 
      seq_along(complete.slices)] = map0 
  
  irows = sort(c(complete.rows, missing.rows), index.return=TRUE)$ix
  icols = sort(c(complete.cols, missing.cols), index.return=TRUE)$ix
  itime = sort(c(complete.slices, missing.slices), index.return=TRUE)$ix
  
  map = map[irows,icols, itime]
  
  map = map[nrow(map):1, , ] # correct lat direction
  map = .rotate(map) # write map for ncdf file, origin in down left corner
  
  dimnames(map) = list(coords$rho$lon, coords$rho$lat, round(times$center,2))
  
  if(toPA) map=.toPA(map)
  
  return(map)
}


# Auxiliar functions for maps ---------------------------------------------

getIsoline = function(lon, lat, z, level, cutoff=0) {
  
  .getLine = function(x) rbind(cbind(x$x, x$y), NA)
  out = grDevices::contourLines(x=lon, y=lat, z=z, nlevels=1, levels=level)
  ll = unlist(lapply(out, function(x) length(x$x)))
  ns = length(ll)
  ll2 = sort(ll)
  cat("Creating ", ns, " isolines of lengths ", paste(ll2[-ns], collapse=", "),
      " and ", ll2[ns], ".\n", sep="")
  cat("Removing ", sum(ll<=cutoff)," isolines of length lower than cutoff=", cutoff, ".\n", sep="")
  out = out[ll>cutoff]
  out = do.call(rbind, lapply(out, FUN=.getLine))
  out = out[-nrow(out),]
  out = as.data.frame(out)
  names(out) = c("lon", "lat")
  
  return(out)
}








extractValidData = function(files, control, var="traj", output=NULL) {
  
  # files is a set of ncdf files paths
  # extract data from variables with the same dimension of control, and
  # matrix with the same lat and lon dimensions. 
  
  cat("Reading control variable...\n")
  control = open.ncdf(control)
  control.var = get.var.ncdf(control, var)
  control.dim = dim(control.var)
  ix = which(!is.na(control.var))
  close.ncdf(control)
  gc(verbose=FALSE)
  base = NULL
  for(file in files) {
    cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      isValidArray  = identical(data1$var[[var]]$varsize, control.dim)
      isValidMatrix = identical(data1$var[[var]]$varsize, control.dim[1:2])
      if( !isValidArray & !isValidMatrix ) {
        cat("Skipping variable", var,"\n")  
        next
      } else {
        
        cat("Adding variable", var,"\n")
        base.names = colnames(base)
        if(isValidMatrix) {        
          base = cbind(base, 
                       array(get.var.ncdf(data1, var), dim=control.dim)[ix])
        }
        if(isValidArray) {
          base = cbind(base, get.var.ncdf(data1, var)[ix])
        }
        colnames(base) = c(base.names, var)
        
      }
      
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  
  cat("Writing final data base.\n")
  base = data.frame(base)
  if(!is.null(output)) write.csv(base, output)
  
  return(base)
}

ncdf2data = function(files, slices, control, var, control.dim=NULL) {
  
  # files is a set of ncdf files paths
  # extract data from variables with the same dimension of control, and
  # matrix with the same lat and lon dimensions. 
  
  cat("Reading control variable...\n")
  if(is.null(control.dim)) {
    control = open.ncdf(control)
    control.dim = control$var[[var]]$varsize
    cat("Extracting data with dimensions", control.dim, "\n")
    close.ncdf(control)
    gc(verbose=FALSE)    
  }
  base = NULL
  for(file in files) {
    cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      isValidArray  = identical(data1$var[[var]]$varsize, control.dim)
      isValidMatrix = identical(data1$var[[var]]$varsize, control.dim[1:2])
      if( !isValidArray & !isValidMatrix ) {
        cat("Skipping variable", var,"\n")  
        next
      } else {
        cat("Adding variable", var,"\n")
        base.names = colnames(base)
        if(isValidMatrix) {
          newvar = array(get.var.ncdf(data1, var), dim=control.dim)[,,slices]
        }
        if(isValidArray) {
          newvar = get.var.ncdf(data1, var)[,,slices]
        }
        dim(newvar) = NULL
        base = cbind(base, newvar)
        colnames(base) = c(base.names, var)        
      }
      
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  
  cat("Writing final data base.\n")
  base = data.frame(base)
  
  return(base)
}

.makeNameNcdfVar = function(var, data) {
  
  longname = data$var[[var]]$longname
  units    = data$var[[var]]$units
  output = paste0(sQuote(var), " (", longname, ", ", units, ")")
  
  return(output)
}

extractData = function(files, data, lat, lon, start, end, 
                       dx, dy=dx, frequency=12, 
                       dim.names = c("lat", "lon", "time"),
                       verbose=TRUE) {
  
  if(dim.names[1] %in% names(data))
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  latAsFactor  = cut(data[, dim.names[1]], coords$psi$lat, labels=FALSE)
  lonAsFactor  = cut(data[, dim.names[2]], coords$psi$lon, labels=FALSE)
  timeAsFactor = cut(data[, dim.names[3]], times$bounds, labels=FALSE)
  
  ix  = cbind(lonAsFactor, latAsFactor, timeAsFactor)
  ix2 = cbind(lonAsFactor, latAsFactor)
  
  control.dim = c(length(coords$rho$lon),
                  length(coords$rho$lat),
                  length(times$center))                  
  
  for(file in files) {
    if(isTRUE(verbose)) cat("Reading", file,"...\n")
    data1 = open.ncdf(file)
    for(var in names(data1$var)) {
      if(!identical(data1$var[[var]]$varsize, control.dim)) {
        if(identical(data1$var[[var]]$varsize, control.dim[1:2])) {
          if(isTRUE(verbose)) cat("\tAdding variable", .makeNameNcdfVar(var, data1),"\n")
          data[, var] = get.var.ncdf(data1, var)[ix2]  
        } else {
          if(isTRUE(verbose)) cat("Skipping variable", sQuote(var),"(Dimensions do not match.)\n")  
          next  
        }
      } else {
        if(isTRUE(verbose)) cat("\tAdding variable", .makeNameNcdfVar(var, data1),"\n")
        data[, var] = get.var.ncdf(data1, var)[ix]
      }
    }
    close.ncdf(data1)
    gc(verbose=FALSE)
  }
  if(isTRUE(verbose)) cat("Writing final data base.\n\n")
  
  return(data)
}

extractEnvData = extractData

extractData2 = function(files, data, lat, lon, start, end, 
                       dx, dy=dx, frequency=12, 
                       dim.names = c("lat", "lon", "time"),
                       verbose=TRUE) {
  
  if(dim.names[1] %in% names(data))
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  latAsFactor  = cut(data[, dim.names[1]], coords$psi$lat, labels=FALSE)
  lonAsFactor  = cut(data[, dim.names[2]], coords$psi$lon, labels=FALSE)
  timeAsFactor = cut(data[, dim.names[3]], times$bounds, labels=FALSE)
  
  ix  = cbind(lonAsFactor, latAsFactor, timeAsFactor)
  ix2 = cbind(lonAsFactor, latAsFactor)
  
  control.dim = c(length(coords$rho$lon),
                  length(coords$rho$lat),
                  length(times$center))                  
  
  for(file in files) {
    if(isTRUE(verbose)) cat("Reading", file,"...\n")
    data1 = nc_open(file)
    for(var in names(data1$var)) {
      if(!identical(data1$var[[var]]$varsize, control.dim)) {
        if(identical(data1$var[[var]]$varsize, control.dim[1:2])) {
          if(isTRUE(verbose)) cat("\tAdding variable", .makeNameNcdfVar(var, data1),"\n")
          data[, var] = ncvar_get(data1, var)[ix2]  
        } else {
          if(isTRUE(verbose)) cat("Skipping variable", sQuote(var),"(Dimensions do not match.)\n")  
          next  
        }
      } else {
        if(isTRUE(verbose)) cat("\tAdding variable", .makeNameNcdfVar(var, data1),"\n")
        data[, var] = ncvar_get(data1, var)[ix]
      }
    }
    nc_close(data1)
    gc(verbose=FALSE)
  }
  if(isTRUE(verbose)) cat("Writing final data base.\n\n")
  
  return(data)
}


table2ncdf = function(file, species, 
                      lat, lon, dx, dy,
                      start, end, frequency=12,
                      toPA=FALSE, longNames=NULL,
                      FUN=mean, colClasses=NA) {
  
  
  cat("Reading data from file", paste0(file,"...", "\n"))
  
  precision = if(toPA) "integer" else "float"
  units     = if(toPA) "1/0" else "abundance"
  unitsName = if(toPA) "Presence/Absence" else "Abundance"
  
  output = sub(pattern=".csv",".nc", file)
  if(is.null(longNames)) longNames = paste(species, unitsName)
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  times  = createTimeAxis(start=start, end=end, frequency=frequency, center=TRUE)
  
  dim1 = dim.def.ncdf("lon",  "degrees", coords$rho$lon)
  dim2 = dim.def.ncdf("lat", "degrees", coords$rho$lat)
  dim3 = dim.def.ncdf("time", "months", times$center)
  
  biology = read.csv(file=file, colClasses=colClasses)
  
  checkNames = all(c("year", "month") %in% names(biology))
  if(!checkNames) stop(paste("Input file has no", 
                             sQuote(year), "or", sQuote(month),
                             "fields."))
  
  biology$time = biology$year + biology$month/12 - 1/24
  biology$traj = 1
  
  species = c(species, "traj")
  longNames = c(longNames, "Survey trajectory")
  
  cat("Creating maps from data...\n")
  
  species.data = list()
  ncdf.vars = list()
  
  for(var in seq_along(species)) {
    cat("Species:", paste0(species[var],"... "))
    species.data[[var]] = table2array(biology, 
                                      var=species[var], 
                                      lat=lat, lon=lon, 
                                      start=start, end=end, 
                                      dx=dx, FUN=FUN, 
                                      toPA=toPA)
    
    ncdf.vars[[var]]    = var.def.ncdf(name=species[var], units=units, 
                                       dim=list(dim1,dim2,dim3), missval=-9999, 
                                       longname=longNames[var], 
                                       prec=precision)
    
    cat("DONE.\n")
  }
  
  cat("Creating ncdf file...\n")
  
  nc.new = create.ncdf(output, ncdf.vars)
  
  for(var in seq_along(species)) {
    # put values to vars
    cat("\t Writing", paste0(species[var],"...\n"))
    put.var.ncdf(nc.new, ncdf.vars[[var]], species.data[[var]])
  }
  
  cat("Saving data.\n")
  
  close.ncdf(nc.new)
  
  cat("File", output, "has been written.","\n")     
  
  return(invisible(output))
}

trainingDatabase = function(files, control, var="traj", 
                            output=NULL) {
  
  extractValidData(files=files, control=control, var=var, 
                   output=output)  
}

writePredictionData = function(files, control, var="traj", path,
                               start, end, frequency=12) {
  
  times  = createTimeAxis(start=start, end=end, frequency=frequency, 
                          center=TRUE)$center
  
  control = open.ncdf(control)
  control.dim = control$var[[var]]$varsize
  cat("Extracting data with dimensions", control.dim, "\n")
  close.ncdf(control)
  gc(verbose=FALSE)
  
  n = ceiling(log(length(times)))
  
  for(t in seq(along=times)) {
    
    cat("----------------------\n")
    cat("----- time step =", t, "\n")
    cat("----------------------\n")
    base = ncdf2data(files = files, slice=t,
                     control = control, 
                     var = var,
                     control.dim=control.dim)
    
    file = file.path(path, sprintf(paste0("base_env_%0",n,"d.csv"),t))
    write.csv(base, file)
    
  }
  
  return(invisible(control.dim[1:2]))
}

setDataBase = function(base, factors) {
  
  base = setFactors(base, factors)
  base = base[complete.cases(base),]
  
  return(base)
}

setFactors = function(base, factors) {
  
  for(f in factors) {
    base[, f] = as.factor(base[, f])
  }
  
  return(base)
}

makePredictions = function(input, output, run, factors) {
  
  root = getwd()
  
  pred.files = dir(path=input)
  
  for(t in seq_along(pred.files)) {
    cat("Running map ", t, ": ", as.character(Sys.time()),"\n", sep="")
    data.file = file.path(input, pred.files[t])
    newdata = read.csv(data.file, row.names=1)
    newdata = setFactors(newdata, factors)
    complete = complete.cases(newdata)
    new.data = newdata[complete,]
    
    setwd(run)
    Projection(Proj=new.data, Proj.name=t, 
               BinRoc=TRUE, BinKappa=TRUE, BinTSS=TRUE)
    
    dir.complete = file.path(getwd(),paste0("proj.", t))
    file.complete = paste("Proj_", t, "_complete", sep="")
    save(complete, file=file.path(dir.complete, file.complete))
    
    setwd(root)
    
    file.copy(from=file.path(dir.complete, dir(path=dir.complete)), 
              to=file.path(output,dir(path=dir.complete)), overwrite=TRUE)
    
    unlink(dir.complete, recursive=TRUE, force=TRUE)
    
  }
  
}

fitModels = function(run, model.session,
                     GLM = FALSE, TypeGLM = "simple", Test = "AIC", GBM = FALSE, 
                     No.trees = 5000, GAM = FALSE, Spline = 3, CTA = FALSE, CV.tree = 50, 
                     ANN = FALSE, CV.ann = 5, SRE = FALSE, quant = 0.025, FDA = FALSE, 
                     MARS = FALSE, RF = FALSE, NbRunEval = 1, DataSplit = 100, 
                     NbRepPA = 0, strategy = "sre", coor = NULL, distance = 0, 
                     nb.absences = NULL, Yweights = NULL, VarImport = 0, Roc = FALSE, 
                     Optimized.Threshold.Roc = FALSE, Kappa = FALSE, TSS = FALSE, 
                     KeepPredIndependent = FALSE) {
  
  cat("Running BIOMOD.\n")
  cat("Starting at", as.character((Sys.time())), "\n")
  
  root = getwd()
  on.exit(setwd(root))
  setwd(run)
  
  Models(GLM = GLM, TypeGLM = TypeGLM, Test = Test, GBM = GBM, 
         No.trees = No.trees, GAM = GAM, Spline = Spline, CTA = CTA, CV.tree = CV.tree, 
         ANN = ANN, CV.ann = CV.ann, SRE = SRE, quant = quant, FDA = FDA, 
         MARS = MARS, RF = RF, NbRunEval = NbRunEval, DataSplit = DataSplit, 
         NbRepPA = NbRepPA, strategy = strategy, coor = coor, distance = distance, 
         nb.absences = nb.absences, Yweights = Yweights, VarImport = VarImport, 
         Roc = Roc, Optimized.Threshold.Roc = Optimized.Threshold.Roc, 
         Kappa = Kappa, TSS = TSS, KeepPredIndependent = KeepPredIndependent)
  
  cat("Ended at", as.character((Sys.time())), "\n")
  save.image(model.session)
  
  setwd(root)
  
  return(invisible())
}

getPredictionMap = function(output, time, species) {
  
  root = getwd()
  setwd(file.path(root,output))
  i = time
  proy = paste0("Proj_", i, "_",species)
  comp = paste0("Proj_", i, "_complete")
  load(proy)
  load(comp)
  proy = get(proy)
  pred = array(dim=c(pred.dim, 9))
  pred[complete] = as.vector(proy[,,1,1])/1000 
  
  return(pred)
}

.getDistance = function (xy, ref, ref2abs, index.return = FALSE) {
  
  .getDist = function(xy, ref, ref2abs, index.return) {
    grados = acos(ref[, 4] * xy[4] + ref[, 3] * xy[3] * cos(ref[,2] - xy[2]))
    ind = which.min(grados)
    dist = c(grados[ind], ref2abs[ind])
    if(isTRUE(index.return)) 
      dist = c(dist, ind)
    return(dist)
  }
  out = t(apply(xy, 1, .getDist, ref = ref, ref2abs = ref2abs, 
                index.return = index.return))
  return(out)
}

.getDistance.old = function(data, ref, ref2abs=NULL, lon="lon", lat="lat") {
  
  x    = as.matrix(data[,c(lat, lon)])
  cx   = as.matrix(ref[,c(lat, lon)])
  
  grados2km  = 1.852*60*180/pi
  grados2rad = pi/180
  
  if(!is.null(ref2abs)) ref2abs = ref2abs/grados2km
  
  ref = cx*grados2rad
  ref = cbind(ref, cos(ref[,1]), sin(ref[,1])) # c(lat, lon, cos(lat), sin(lat))
  xy  = x*grados2rad  
  xy  = cbind(xy, cos(xy[,1]), sin(xy[,1]))
  
  output = array(dim=c(nrow(data), 1 + !is.null(ref2abs)))
  ind = complete.cases(xy)
  output[ind,] = .getDistance(xy=xy[ind, ], ref=ref, ref2abs=ref2abs)*grados2km
  
  return(output)
}

getDistance = function (data, ref, ref2abs = NULL, lon = "lon", lat = "lat", 
                        index.return=FALSE) {
  
  data     = data[,c(lat,lon)]
  dat      = data[!duplicated(data),]
  data$ind = 1:nrow(data)
  
  x = as.matrix(dat[, c(lat, lon)])
  cx = as.matrix(ref[, c(lat, lon)])
  grados2km = 1.852 * 60 * 180/pi
  grados2rad = pi/180
  if (!is.null(ref2abs)) ref2abs = ref2abs/grados2km
  ref = cx * grados2rad
  ref = cbind(ref, cos(ref[, 1]), sin(ref[, 1]))
  xy = x * grados2rad
  xy = cbind(xy, cos(xy[, 1]), sin(xy[, 1]))
  ncols = 1 + isTRUE(index.return) + !is.null(ref2abs)
  output = array(dim = c(nrow(dat), ncols))
  ind = which(complete.cases(xy))
  output[ind, ] = .getDistance(xy = xy[ind, ], ref = ref, 
                               ref2abs = ref2abs,
                               index.return=index.return) 
  output[, 1] = output[, 1] * grados2km
  if(!is.null(ref2abs)) output[, 2] = output[, 2] * grados2km
  
  xnames = if(!is.null(ref2abs)) c("dist", "ref2abs") else "dist"
  xnames = if(isTRUE(index.return)) c(xnames, "ix") else xnames
  colnames(output) = xnames
  output = data.frame(dat, output)
  output = merge(data, output, all=TRUE, sort=FALSE)
  ind    = sort(output$ind, index=TRUE)$ix
  output = drop(as.matrix(output[ind, xnames]))
  
  return(output)
}

getSignedDistance = function(data, ref, abs, lon="lon", lat="lat", digits=3) {
  
  grados2km  = 1.852*60*180/pi
  
  DateStamp("Calculating distance from data to abs.\t")
  data2abs = getDistance(data=data, ref=abs, lon=lon, lat=lat)
  DateStamp("Calculating distance from ref to abs.\t")
  ref2abs  = getDistance(data=ref, ref=abs, lon=lon, lat=lat)
  DateStamp("Calculating distance from data to ref.\t")
  data2ref = getDistance(data=data, ref=ref, ref2abs=ref2abs, lon=lon, lat=lat)
  
  sign = -2*(data2abs < data2ref[,2]) + 1 # inside is negative
  
  dist = sign*data2ref[,1]
  
  output = list(dist=dist, abs=data2abs)
  
  return(output)
}

read.Herve = function(file, skip=9, units=TRUE, k=1000, output=".") {
  
  .parseLine = function(line) {
    x = unlist(strsplit(line," "))
    x = x[x!=""]
    return(x)
  }
  
  .parseHeader = function(line) {
    cat("\t Parsing header...\n")
    var = paste(strsplit(file, "\\.")[[1]][c(5:4)],collapse=".")
    var = paste0(var,".")
    header = .parseLine(dat[skip+1])
    header = sub("\\.", "",tolower(header))
    header = sub("longitude", "lon", header)
    header = sub("latitude", "lat", header)
    header = sub("val", "aver", header)
    header = sub("stddev", "sd", header)
    ind = header=="aver"
    rmax = sum(ind)
    header[ind] = paste0(var, "avg_r", seq_len(rmax)-1)
    ind = header %in% c("min", "max", "sd", "mode")
    header[ind] = paste0(var, header[ind], "_r", rmax-1)
    return(header)
  }
  
  .readHerve = function(line, ncols) {
    x = .parseLine(line)
    if(length(x)!=ncols) {
      out = c(x[2:6], rep(NA, ncols-7))
    } else {
      x[x=="********"] = NA
      x[is.nan(x)] = NA
      out = x[-c(1,ncols)]
    }
    out = as.numeric(out)
    names(out) = NULL
    return(out)
  }
  
  cat("\t Reading file", file, "\n")
  dat = readLines(file)
  header = .parseHeader(dat[skip+1])
  ncols = length(header)
  cat("\t\t",header,"\n")
  if(units) unit = .parseLine(dat[skip+1+units])
  dat = dat[-seq_len(skip+1+units)]
  
  cat("\t Parsing data...\n")
  x4 = array(dim=c(length(dat), ncols-2))
  pp = as.integer(seq(from=0, to=nrow(x4), len=k+1))
  
  for(j in seq_len(k)) {
    ind = seq.int(pp[j]+1, pp[j+1])
    x4a = lapply(dat[ind], .readHerve, ncols=ncols)
    x4[ind,] = t(matrix(unlist(x4a), nrow=ncols-2))
  }
  
  cat("\t Creating data frame...\n")
  x4 = as.data.frame(x4)
  names(x4) = header[-c(1,ncols)]
  cat("\t Done.\n")
  
  file.out = paste(strsplit(file, "\\.")[[1]][c(1,4:6,8,2)],collapse=".")
  file.out = file.path(output, file.out)
  cat("\t Saving", sQuote(file.out), "...\n")
  write.csv(x4, file.out, row.names=FALSE)
  
  return(invisible(file.out))
}

DateStamp = function(...) cat(..., "\t\t | ", date(), "\n\n")

isInside = function(x, range, value=FALSE, prob=FALSE) {
  ind = (x >= range[1]) & (x <= range[2])
  if(isTRUE(prob)) {
    if(all(is.na(x))) return(NA)
    return(sum(ind, na.rm=TRUE)/length(ind)) 
  }
  if(isTRUE(value)) return(x[which(ind)]) 
  return(ind)
}


limitingFactor = function(data, ranges) {
  vars = names(data)[names(data) %in% names(ranges)]
  ranges = ranges[vars]
  .getProp = function(var, data, ranges) 
    isInside(x=data[, var], range=ranges[[var]], prob=TRUE)
  out = sapply(vars, FUN=.getProp, data=data, ranges=ranges)
  return(out)
}

geq = function(x, thr) {
  
  ind = which(x>=thr & !is.na(x))
  out = x[ind]
  names(out) = ind
  
  return(out)  
}

leq = function(x, thr) {
  
  ind = which(x<=thr & !is.na(x))
  out = x[ind]
  names(out) = ind
  
  return(out)  
}

.insider = function(data, control, var, alpha, lowerOnly) {
  
  # TO_DO: handle factors
  delta = alpha/2
  alpha = if(!lowerOnly) c(delta, 1 - delta) else c(1-alpha, 1)
  alpha = quantile(control[,var], prob=alpha, na.rm=TRUE)
  ind = (data[,var] > alpha[1]) & (data[,var] < alpha[2])
  
  return(ind)
}




setPredictionFiles = function(prefix, start, end, aux=NULL, dir=".") {
  
  times = createTimeAxis(start=start, end=end, frequency=12, center=TRUE)
  files = paste0("Y",times$year, "M", sprintf("%02d", times$month))
  files = paste(prefix, files, sep="-")
  files = paste0(files, ".csv")
  files = file.path(dir, files)
  output = list(files=files, aux=aux, time=times)
  files.check = file.exists(files)
  aux.check = if(!is.null(aux)) file.exists(aux) else TRUE
  if(!all(files.check)) warning("Some prediction files don't exist")
  if(!all(aux.check)) warning("Auxiliar prediction files don't exist")
  
  return(output)
}

.getPngFileName = function(filename, replacement) {
  
  if(is.null(filename)) filename = "map"
  filename = gsub(pattern="\\.png", "", filename)
  filename = paste0(paste(filename, replacement, sep="-"), ".png")
  
  return(filename)
}

saveAnimation = function(object, ...) {
  
  UseMethod("saveAnimation")
}

saveAnimation.default = function(object, file, dir=getwd(),
                                 interval=0.5, ...) {
  
  n = dim(object)[3]
  DateStamp("\nCreating animation (",n," time steps).", sep="")
  
  try(suppressMessages(saveGIF(
    {
      for(i in seq_len(n)) {
        cat(i,"")
        image.plot(object[,,i], ...)
      }
      },
    movie.name="temp.gif",
    img.name="slice", clean=TRUE, verbose=FALSE, 
    interval=interval, loop=1, check=TRUE, autobrowse=FALSE)),
    silent=TRUE)
  tmp = file.path("temp.gif")
  x = file.copy(from=tmp, to=file.path(dir, file), overwrite=TRUE)
  DateStamp("DONE.")
  
  return(invisible(x))

}

addLHT = function(object, sp, lifespan=1, 
                  ages=seq_len(ceiling(lifespan))-1, ...) {
  
  pars = list(sp=sp, lifespan=lifespan, ages=ages, ...)
  object$info$LHT = pars
  
  return(invisible(object))
}

normalize = function(x) {
  
  if(length(dim(x))> 3) stop("'x' cannot be an array of dimension greater than 3.")
  
  if(length(dim(x))< 3) out = x/sum(x,na.rm=TRUE) 
  if(length(dim(x))==3) out = apply(x, 3, normalize)
  dim(out) = dim(x)
  return(out)
  
}

.inArea = function(x, lat, lon) {
  
  cnull = c(-Inf, Inf)
  if(is.null(lat)) lat = cnull
  if(is.null(lon)) lon = cnull
  ilat = x$lat >= lat[1] & x$lat <= lat[2]
  ilon = x$lon >= lon[1] & x$lon <= lon[2]
  ind = ilat & ilon
  
  return(ind)

}

.ponderate = function(w, x) {
  
  ind = complete.cases(cbind(as.numeric(w), as.numeric(x)))
  w = w[ind]
  x = x[ind]
  out = weighted.mean(x=x, w=w, na.rm=TRUE)
  
  return(out)
  
}

ponderate = function(x, aux, vars=names(aux), lat=NULL, lon=NULL) {
  
  DateStamp("Starting calculations...")
  vv = vars %in% names(aux)
  if(sum(vv)==0) stop("No valid variables 'vars' in auxiliar file.")
  vars = vars[vv]
  
  dw = dim(x)
  if(length(dw)<2) stop("'w' must be an array or matrix")
  ncol = if(length(dw)==2) 1 else dw[3]
  w = matrix(normalize(x), ncol=ncol)
  
  if(!is.null(lat) | !is.null(lon)) {
    if(!is.null(lat)) lat = sort(lat) 
    if(!is.null(lon)) lon = sort(lon) 
    ind = .inArea(aux, lat=lat, lon=lon)
    aux = aux[ind,]
    w   = w[ind,]    
  }
  
  output = NULL
  for(var in vars) {
    DateStamp("Computing statistics for ", sQuote(var), "...", sep="")
    output$center[[var]] = apply(w*aux[,var], 2, .na.sum)
    output$density[[var]] = getDensity(x=aux[,var], w=w)
    
  }
  
  output$center = as.data.frame(output$center)
  DateStamp("DONE.")
  
  return(output)
  
}

.na.sum = function(x) if(!all(is.na(x))) sum(x, na.rm=TRUE) else NA

.getDensity = function(w, x, n=256, adjust=1.5, range=NULL, nLim=10) {
  
  if(all(is.na(w)) || all(is.na(x))) return(rep(NA, length=n))
  if(sum(w!=0)<= nLim) return(rep(NA, length=n))
  if(sum(w)==0) return(rep(0, length=n))
  ran = if(is.null(range)) range(pretty(x), na.rm=TRUE) else range
  ind = complete.cases(data.frame(w,x))
  x = x[ind]
  w = normalize(w[ind])
  out = density(x=x, weights=w, n=n, adjust=adjust,
                from=ran[1], to=ran[2])
  out = out$y
  
  return(out)
}

getDensity = function(x, w, n=256, adjust=1.5, range=NULL, nLim=10) {
  
  out = apply(w, 2, .getDensity, x=x, n=n, adjust=adjust, range=range, nLim=nLim)
  ale = if(is.null(range)) range(pretty(x)) else range
  ale = seq(from=ale[1], to=ale[2],len=n)
  out = list(xaxis=ale, densities=out)
  
  return(out)
}

getProfiles = function(var, w, by, data, ref=NULL, n=256, nLim=10, 
                       adjust=1.5, alpha=NULL) {
  
  if(!is.null(ref)) {
    if(!(ref %in% w)) stop("Reference model ", sQuote(ref), " not found in model list 'w'.") 
  }
  if(!is.null(alpha)) {
    if(length(alpha)==1) alpha = rep(alpha, len=length(w))
    if(length(alpha)!=length(w)) 
      stop("Parameter alpha doesn't match the length of w")
  } else {
    alpha = rep(0.5, len=length(w))
  }
  validModels = w %in% names(data)
  if(length(validModels)!=length(w)) {
    msg = paste("Models", w[!validModels], "are ignored.")
    warning(msg)
  }
  w     = w[validModels]
  alpha = alpha[validModels]
  ind   = !duplicated(w)
  w     = w[ind]
  alpha = alpha[ind]
  w     = unique(c(ref, w))
  
  sortBy = sort(unique(data[, by]))
  ale = data[, var]
  ali = range(pretty(ale), na.rm=TRUE)
  ale = seq(from=ali[1], to=ali[2],len=n)
  
  output = list()
  output$x = ale
  output$by = sortBy
  output$n = numeric()
  output$npos = array(dim=c(length(sortBy), length(w)+1))
  patterns = array(NA, dim=c(length(sortBy),length(w), 7))
  profiles = list()
  colnames(output$npos) = c("nobs", w)
  rownames(output$npos) = sortBy
  
  for(i in seq_along(w)) {
    model = w[i]
    profiles[[model]] = list()
    for(j in seq_along(sortBy)) {
      dat = data[data[, by]==sortBy[j], ]
      output$n[j] = nrow(dat)
      xx = dat[, var]
      ww = dat[, model]
      npos = sum(ww>0, na.rm=TRUE)
      output$npos[j, i+1] = npos
      nx = if(!is.null(ref)) sum(dat[, ref]>0, na.rm=TRUE) else npos  
      thisProfile = if( (nrow(dat) >= 5*nLim) & (nx > nLim) ) {
        .getDensity(x = xx, w = ww, n=n, range=ali, adjust=adjust, nLim=0)
      } else rep(NA, length=n)
      
      patterns[j, i, ] = .sevennum(x=ale, w=thisProfile)
      profiles[[model]][[j]] = thisProfile
      
    }
    
    profiles[[model]] = as.data.frame(profiles[[model]])
    names(profiles[[model]]) = sortBy
    
  }
  
  patNames = c("mean", "median", "mode", "P05", "P25", "P75", "P95")
  dimnames(patterns) = list(sortBy, w, patNames)
  
  if(!is.null(ref)) {
    output$patterns$observed = patterns[, 1,]
    output$patterns$models   = patterns[,-1,]
    output$profiles$observed = profiles[[ref]]
    output$profiles$models   = profiles[-1]
    
  } else {
    output$patterns$observed = NULL
    output$patterns$models   = patterns
    output$profiles$observed = NULL
    output$profiles$models   = profiles
  }
  
  output$npos[, 1] = output$n
  output$npos = as.data.frame(output$npos)
  
  return(output)
  
}

.sevennum = function(x, w) {
  if(all(is.na(w))) return(rep(NA, 7))
  w = w/sum(w, na.rm=TRUE)
  cw = round(cumsum(w), 3)
  ind = seq_len(sum(cw<1)+1)
  ww = cw[ind]
  xx = x[ind]
  prob = function(p) {
    fun = splinefun(x=ww, y = xx)
    out = fun(p) 
    return(out)
  }
  out = NULL
  out["mean"] = weighted.mean(x, w)
  out["median"] = prob(0.50)
  out["mode"] = x[which.max(w)]
  out["p05"] = prob(0.05)
  out["p25"] = prob(0.25)
  out["p75"] = prob(0.75)
  out["p95"] = prob(0.95)
  
  return(out)
  
}

getRleIndex = function(x) {
  
  x = rle(x)
  nslice = length(x$length)
  nrow   = cumsum(x$lengths)
  nrow   = cbind(begin=c(1, nrow[-nslice] + 1), end=nrow)
  
  return(nrow)
  
}

getValuesinRange = function(index, x, FUN, ...) {
  
  FUN = match.fun(FUN)
  .FUN = function(index, x, ...) FUN(x[index[1]:index[2]], ...)
  output = apply(index, 1, .FUN, x=x, ...)
  
  return(output)
}

getYearMin = function(object, index) {
  
  ind = getRleIndex(object$info$time[[index]])
  out = getValuesinRange(index=ind, x=object$info$time$year, FUN=min)
  out = out - min(out, na.rm=TRUE)
  
  return(out)
}

getYearMax = function(object, index) {
  
  ind = getRleIndex(object$info$time[[index]])
  out = getValuesinRange(index=ind, x=object$info$time$year, FUN=max)
  out = out - min(out, na.rm=TRUE) + 1
  
  return(out)
}

getSeasonalMap = function(object) {
  
  nrow = getRleIndex(object$info$time$season)
  map  = object$prediction
  nslice = nrow(nrow)
  output = array(dim=c(dim(map)[1:2], nslice))
  
  for(i in seq_len(nslice)) {
    ali = nrow[i,1]:nrow[i,2]
    ale = map[,, ali]
    output[,,i] = apply(ale, 1:2, mean, na.rm=TRUE)*object$info$mask
  }
  
  return(output)
}

year2month = function(data, year=seq_along(data), month) {
  
  var.name = rev(strsplit(deparse(substitute(data)), "\\$")[[1]])[1]
  years = rep(year, each=12)
  months = rep(1:12, len=length(years))
  new.data = rep(NA, len=length(years))
  new.data[months==month] = data
  output = data.frame(year=years, month=months, new.data)
  names(output)[3] = var.name
  
  return(output)
}

.removeSector3 = function(x, coords, lon, lat) {
  
  xlat = isInside(coords$lat, lat)
  xlon = isInside(coords$lon, lon)
  xpred = x
  
  if(!is.null(lat) & !is.null(lon)) {
    xpred[xlon, xlat, ] = 0
  } else {
    if(is.null(lon)&!is.null(lat)) {
      xpred[,xlat,] = 0
    }
    if(!is.null(lon)&is.null(lat)) {
      xpred[xlon,,] = 0      
    }
  }
  xpred[is.na(x)] = NA
  return(xpred)
  
}

.removeSector2 = function(x, coords, lon, lat) {
  
  xlat = isInside(coords$lat, lat)
  xlon = isInside(coords$lon, lon)
  xpred = x
  
  if(!is.null(lat) & !is.null(lon)) {
    xpred[xlon, xlat] = 0
  } else {
    if(is.null(lon)&!is.null(lat)) {
      xpred[,xlat] = 0
    }
    if(!is.null(lon)&is.null(lat)) {
      xpred[xlon,] = 0      
    }
  }
  xpred[is.na(x)] = NA
  
  return(xpred)
}

removeSector = function(x, coords, lon, lat) {
  
  ndim = length(dim(x))
  if(ndim==2) out = .removeSector2(x, coords, lon, lat)
  if(ndim==3) out = .removeSector3(x, coords, lon, lat)
  
  return(out)
}

makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
}

.removeVars = function(x, vars) {
  
  vars = intersect(vars, names(x))
  if(length(vars)==0) return(x)
  
  for(var in vars) {
    x[, var] = NULL
  }
  
  return(x)
}

createPredictionFiles = function(files, aux.files=NULL, lat, lon, start, end, 
                                 dx, dy=dx, frequency=12, output=".",
                                 FUN=NULL, FUN2=NULL, 
                                 verbose=FALSE, prefix="base-env", ...) {
  
  if(!is.null(FUN)) FUN = match.fun(FUN)
  
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  DateStamp("Starting at")
  
  if(!is.null(aux.files)) {
    DateStamp("Creating auxiliar prediction file...")
    xy = createAuxPredictionFiles(files=aux.files, lat=lat, lon=lon, 
                                  start=start, end=end, dx, dy=dx, 
                                  frequency=frequency, output=output,
                                  FUN=FUN2, verbose=TRUE, 
                                  prefix=prefix, ...)
    auxVars = setdiff(names(xy), c("lat", "lon"))
  } else {
    coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
    xy = data.frame(lat=as.numeric(coords$LAT), lon=as.numeric(coords$LON))
    xy$lon = round(xy$lon, 3)
    xy$lat = round(xy$lat, 3)
  }
  
  day = 15 # to do
  times  = createTimeAxis(start=start, end=end, center=TRUE, frequency=frequency)
  dates = data.frame(year=times$year, month=times$month)
  dates$time = dates$year + (dates$month-1)/12 + (day-1)/365
  
  for(t in seq_len(nrow(dates))) {
    DateStamp("Time step", t, ":", as.numeric(dates[t,1:2]), "\t\t")
    xy$time =  dates$time[t]
    xy$month = dates$month[t]
    out = extractData2(files=files, data=xy, 
                      lat=lat, lon=lon, start=start, end=end, dx=dx, dy=dy,
                      verbose=verbose)
    
    if(!is.null(FUN)) out = FUN(out, ...)
    
    out$time = NULL
    if(!is.null(aux.files)) out = .removeVars(out, vars=auxVars)
    
    file.out = paste0(prefix, "-Y",dates$year[t], "M",
                      sprintf("%02d",dates$month[t]),".csv")
    
    write.csv(out, file.path(output, file.out), row.names=FALSE)
  }
  
  DateStamp("Finishing at") 
  
  return(invisible())
}

createAuxPredictionFiles = function(files, lat, lon, start, end, 
                                    dx, dy=dx, frequency=12, output=".",
                                    FUN=NULL, verbose=FALSE, 
                                    prefix="base-env-auxiliar", ...) {
  
  if(!is.null(FUN)) FUN = match.fun(FUN)
  if(!file.exists(output)) dir.create(output, recursive=TRUE)
  
  
  coords = createGridAxes(lat=lat, lon=lon, dx=dx, dy=dy)
  
  out      = data.frame(lat=as.numeric(coords$LAT), lon=as.numeric(coords$LON))
  out$lon  = round(out$lon, 3)
  out$lat  = round(out$lat, 3)
  out$time =  start[1] + 14/365
  
  out = extractData2(files=files, data=out, 
                    lat=lat, lon=lon, start=start, end=end, dx=dx, dy=dy,
                    verbose=verbose)
  
  if(!is.null(FUN)) {
    DateStamp("Post-processing auxiliar file.")
    out = FUN(out, ...)
  }
  out$time = NULL
  
  file.out = paste0(prefix,"-auxiliar.csv")
  
  write.csv(out, file.path(output, file.out), row.names=FALSE)
  
  DateStamp("Finishing at") 
  
  return(invisible(out))
}



