
AUC = function(data, coordNames = c("lat", "lon"), obs="observed",
               models=NULL, st.dev=TRUE, na.rm=TRUE) {

  modelNames = models
  if(is.null(models)) {
    models = which(!(names(data) %in% c(coordNames, obs)))
    modelNames = names(data)[models]
  }
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  output = NULL
  for(i in seq_along(models)) {
    output = rbind(output, auc(DATA,st.dev=st.dev, which.model=i,na.rm=na.rm))
  }
  rownames(output) = modelNames
  
  return(output)
}

kappa = function(data, coordNames = c("lat", "lon"), obs="observed",
                 models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  modelNames = models
  if(is.null(models)) {
    models = which(!(names(data) %in% c(coordNames, obs)))
    modelNames = names(data)[models]
  }
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  
  thr = as.numeric(optimal.thresholds(DATA=DATA, opt.methods="MaxKappa")[-1])
  output = NULL
  for(i in seq_along(models)) {
    output = rbind(output, Kappa(cmx(DATA,threshold=thr[i], which.model=i,
                                     na.rm=na.rm), st.dev=st.dev))
  }
  rownames(output) = modelNames
  
  return(output)
}

TSS = function(data, coordNames = c("lat", "lon"), obs="observed",
               models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  modelNames = models
  if(is.null(models)) {
    models = which(!(names(data) %in% c(coordNames, obs)))
    modelNames = names(data)[models]
  }
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  
  thr = as.numeric(optimal.thresholds(DATA=DATA, opt.methods="MaxSens+Spec")[-1])
  output = NULL
  for(i in seq_along(models)) {
    CMX = cmx(DATA,threshold=thr[i], which.model=i, na.rm=na.rm)
    output = rbind(output, cbind(specificity(CMX, st.dev=st.dev), 
                                 sensitivity(CMX, st.dev=st.dev)))
  }
  if(st.dev) {
    tss = cbind(output[,1]+output[,3]-1, sqrt(output[,2]^2+output[,4]^2))
    colnames(tss) = c("TSS", "TSS.sd")
  } else {
    tss = cbind(output[,1]+output[,2]-1)
    colnames(tss) = c("TSS")
  }
  output = cbind(tss, output)
  rownames(output) = modelNames
  
  return(output)
}

PredictivePerformance = function(data, coordNames = c("lat", "lon"), obs="observed",
                                 models=NULL, st.dev=TRUE, na.rm=TRUE) {
  
  data = data[complete.cases(data), ]
  out1 = AUC(data=data, coordNames = coordNames, obs=obs, 
             models=models, st.dev=TRUE, na.rm=na.rm)
  out2 = kappa(data=data, coordNames = coordNames, obs=obs, 
               models=models, st.dev=TRUE, na.rm=na.rm)
  out3 = TSS(data=data, coordNames = coordNames, obs=obs, 
             models=models, st.dev=TRUE, na.rm=na.rm)
  
  output1 = cbind(AUC=out1[,1], Kappa=out2[,1], out3[,c(1,3,5)])
  output2 = cbind(AUC.sd=out1[,2], Kappa.sd=out2[,2], out3[,c(2,4,6)])
  
  output = if(st.dev) {
    list(statistic=output1, sd=output2)
  } else {
    output1
  }
  
  return(output)
}

plotROC = function(data, coordNames = c("lat", "lon"), obs="observed",
                   models=NULL, opt.thresholds=TRUE, opt.methods=c(4,10),
                   req.sens=0.95, ...) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  auc.roc.plot(DATA=DATA, opt.thresholds=opt.thresholds, 
               opt.methods=opt.methods, req.sens=req.sens, ...)
  
  return(invisible())
}

plotThreshold = function(data, coordNames = c("lat", "lon"), obs="observed",
                         models=NULL, opt.thresholds=TRUE, opt.methods=c(4,10),
                         req.sens=0.95, which.model=1, ...) {
  
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  error.threshold.plot(DATA=DATA, which.model=which.model, opt.thresholds=opt.thresholds, 
                       opt.methods=opt.methods, req.sens=req.sens, ...)
  
  return(invisible())
}

calculateThresholds = function(data, coordNames = c("lat", "lon"), 
                               obs="observed", 
                               models=NULL, opt.methods=2:12,
                               req.sens=0.95, req.spec=0.5,
                               FPC=1, FNC=10, ...) {
  
  data = data[complete.cases(data), ]
  if(is.null(models)) models = !(names(data) %in% c(coordNames, obs))
  observed = data[,obs]
  if(is.factor(observed)) observed = as.numeric(as.character(observed))  
  DATA = data.frame(1, observed, data[,models])
  out = optimal.thresholds(DATA=DATA, opt.methods=opt.methods, 
                           req.sens=req.sens, req.spec=req.spec,
                           FPC=FPC, FNC=FNC, ...)
  output = as.matrix(out[,-1])
  rownames(output) = out[,1]
  
  return(output)
}

cleanBIOMOD = function() {
  
  if(exists("Biomod.material", env=globalenv())) rm(list="Biomod.material", envir=globalenv())
  if(exists("DataBIOMOD", env=globalenv())) rm(list="DataBIOMOD", envir=globalenv())
  if(exists("DataEvalBIOMOD", env=globalenv())) rm(list="DataEvalBIOMOD", envir=globalenv())
  
  return(invisible())
}

table.summary = function(..., sp) {
  warning("table.summary function is deprecated, '' instead")
  return(invisible(NULL))
}

xxx = function(..., sp) {
  
  names  = as.character(substitute(list(...))[-1L])
  object = list(...)
  output = NULL
  tnames = NULL
  for(i in seq_len(length(object))) {
    tnames = c(tnames, levels(as.factor(object[[i]][,sp])))
  }
  tnames = unique(tnames)
  output = array(0, dim=c(length(object), length(tnames)))
  colnames(output) = tnames
  rownames(output) = names
  for(i in seq_len(length(object))) {
    tnames = levels(as.factor(object[[i]][,sp]))
    output[i,tnames] = table(object[[i]][,sp])
  }
  total = apply(output, 1, sum, na.rm=TRUE)
  percent = round(output/total,3)
  output = cbind(output, total, percent)
  
  return(output)
}


splitDataSet = function(data, var, factor=0.15, seed=771104) {
  
  if(class(data)=="niche.models") {
    output = list(train=data$train, val=data$val)
  }
  if(is.data.frame(data)) {
    data = data[!is.na(data[,var]),]
    dataN = data[data[,var]==0,]
    dataP = data[data[,var]==1,]
    nn = nrow(dataN)
    np = nrow(dataP)
    set.seed(seed)
    indn = sample(nn, factor*nn)
    indp = sample(np, factor*np)
    training = rbind(dataN[-indn,], dataP[-indp,])
    validation = rbind(dataN[indn,], dataP[indp,])
    output = list(train=training, val=validation)
  }
  class(output) = c("niche.models", class(output))
  
  return(output)    
}

.gam.fmla = function(y, var, spline=NULL, loess=NULL, factors=NULL) {
  
  if(is.null(loess)) {
    loess = rep(FALSE, length(var))
    names(loess) = var
  }
  
  #   ss = c(sst=2, sss=2, chl=2)
  #   lo = c(sst=FALSE, sss=FALSE, chl=TRUE)
  #   
  #   ft = rep("s", length(var))
  #   ft[lo] = "lo"
  #   lt = rep(",%02d", length(var))
  #   lt[lo] = ""
  #   x = paste0(ft, "(",var, lt, ")")
  #   x = paste(x, collapse=" + ")
  
  if(!is.null(spline)) {
    fmla = as.list(spline[var])
    fmla = append(fmla, paste0("s(", var, ",%d)", collapse=" + "), after=0)
    fmla = do.call(sprintf, fmla)
  } else {
    fmla = paste0(var, collapse=" + ")
  }
  
  if(!is.null(factors)) {
    
    fmla = paste0(c(fmla, factors), collapse=" + ") 
    
  }
  
  fmla = paste0(y, " ~ ",fmla)
  
  return(fmla)
}

gam.fmla = function(y, var, spline=NULL, factors=NULL) {
  
  fmla = .gam.fmla(y=y, var=var, spline=spline, factors=factors)
  fmla = as.formula(fmla, env=.GlobalEnv)
  
  return(fmla)
}

.fmla2txt = function(fmla) {
  
  out = paste(as.character(fmla)[c(2,1,3)], collapse=" ")
  return(out)
}

.getObjectName = function(x) {
  out = deparse(substitute(x))
  out = deparse(eval(substitute(quote(x))))
  
  return(out)
}

fakeGAM = function(object) {
  
  obj = deparse(eval(substitute(quote(object), env=as.environment(-1))))
  DateStamp("Fitting models for", sQuote(obj), "dataset.")
  
  return(obj)
}

fakeGAM2 = function(object) {
  
  obj = deparse(substitute(object, env=as.environment(-1)))
  DateStamp("Fitting models for", sQuote(obj), "dataset.")
  
  return(obj)
}



fitGAMs = function(object, formulas, FUN=identity, 
                   name=NULL, link="logit", bigData=FALSE) {
  
  name = deparse(substitute(object))
  
  FUN = match.fun(FUN)
  object$transform = FUN
  
  gc(verbose=FALSE)
  
  DateStamp("Fitting models for", sQuote(name), "dataset.")
  
  train = FUN(object$train)
  val   = FUN(object$val)
  
  if(!is.list(object$predicted)) object$predicted = NULL
  
  models = names(formulas)
  
  aic = object$fit[,"AIC"]
  bic = object$fit[,"BIC"]
  var = as.character(formulas[[1]])[2]
  
  object$predicted$lon = object$train[,"lon"]
  object$predicted$lat = object$train[,"lat"]
  object$predicted$observed = as.numeric(as.character(object$train[,var]))
  object$predicted = as.data.frame(object$predicted)
  
  object$validation$lon = object$val[,"lon"]
  object$validation$lat = object$val[,"lat"]
  object$validation$observed = as.numeric(as.character(object$val[,var]))
  object$validation = as.data.frame(object$validation)
  
  for(i in seq_along(formulas)) {
    model.name = models[i]
    model.formula = formulas[[i]]
    DateStamp("Training model ", model.name, ":\n", .fmla2txt(formulas[[i]]), sep="")
    object$formulas[[model.name]] = model.formula 
    #model.vars = .getModelVars2(model.formula, train)
    # TO_DO: filter complete cases
    if(isTRUE(bigData)) {
      model = mgcv::bam(model.formula, data = train, family = binomial(link=link))
    } else {
      model = mgcv::gam(model.formula, data = train, family = binomial(link=link))
    }
    
    gc(verbose=FALSE)
    model$anova = anova(model)
    model$call$family[2] = link
    object$models[[model.name]] = model
    # object$preplot[[models[i]]] = with(object$train, preplot(model))
    
    aic[model.name] = AIC(model)
    bic[model.name] = BIC(model)
    object$predicted[, model.name] = predict(model, newdata=train, type="response")
    object$validation[, model.name] = predict(model, newdata=val, type="response")
  }
  
  DateStamp("Computing Predictive Performance...")
  
  object$fit = cbind(AIC=aic, BIC=bic)
  
  object$performance$training   = PredictivePerformance(object$predicted, st.dev=FALSE)
  object$performance$validation = PredictivePerformance(object$validation, st.dev=FALSE)
  
  object$threshold$training  = calculateThresholds(object$predicted)
  object$threshold$validation = calculateThresholds(object$validation)
  
  DateStamp("DONE.")
  
  class(object$formulas) = c("niche.models.formulas", class(object$formulas))
  
  return(object)
}




getPredictions = function(data) {
  
  if(is.list(data$predicted)) return(do.call(cbind, data$predicted)) else return(invisible(NULL))
}

print.niche.models.formulas = function(x, ...) {
  
  out = matrix(sapply(x, .fmla2txt), ncol=1)
  colnames(out) = ""
  rownames(out) = names(x)
  print(out)
  
  return(out)
}

print.niche.models = function(x, ...) {
  
  model.group = deparse(substitute(x))
  best = .getBestModel(x)
  bfml = .fmla2txt(x$formulas[[best]])
  cat(model.group, "\n\tBest Model:", best, ":\n", bfml)
  
  return(invisible())
}

.getModelVars2 = function(formula, data) {
  
  ale  = as.character(formula)[3]
  vars = names(data)
  vars  = vars[unlist(lapply(vars, grepl, x=ale))]
  
  return(vars)  
}

.getModelVars = function(model) {
  
  ale  = as.character(model$formula)[3]
  vars = names(model$data)
  vars  = vars[unlist(lapply(vars, grepl, x=ale))]
  
  return(vars)  
}

predict.niche.models = function(object, model=NULL, pred.info, mask=NULL, radius=1,
                                precision=3, interp=!is.null(mask), FUN=identity,
                                req.sens=0.90, cost=list(FPC=1, FNC=10), ...) {
  
  if(!is.null(object$transform)) FUN0 = match.fun(object$transform)
  FUN  = match.fun(FUN)
  
  factor = 10^max(trunc(precision), 3)
  if( !is.null(model) && !is.na(model) ) {
    model.name = model
    model = getModel(object, model)
  } else {
    model.name = .getBestModel(object)
    model = getBestModel(object)
  }
  
  DateStamp("Calculating limiting variables")
  # limiting variables
  values = lapply(model$var.summary, xrange, n=10)
  newdata = do.call(expand.grid, values)
  newdata$niche = predict(model, newdata = newdata, type="response")
  species = as.character(model$formula[2])
  model$model$fitted = predict(model, newdata = model$model, type="response")
  thr = try(calculateThresholds(data=model$model, coordNames=names(model$var.summary),
                                models="fitted", obs=species, req.sens=req.sens, FPC=cost$FPC, FNC=cost$FPC))
  values = lapply(newdata[newdata$niche>thr["ReqSpec", 1], ], quantile, 
                  probs=c(0.025, 0.975))
  xvar = which(names(newdata) %in% names(values))
  limiting = matrix(nrow=length(pred.info$files), ncol=length(xvar)-1)
  
 
  # end limiting variables
  
  if(!is.null(pred.info$aux)) aux = read.csv(pred.info$aux)
  
  for(i in seq_along(pred.info$files)) {
    
    DateStamp("Prediction on ", pred.info$files[i], "...", sep="")
    
    newdata = read.csv(pred.info$files[i])
    
    if(i==1) {
      
      nlon = length(unique(newdata$lon))
      nlat = length(unique(newdata$lat))
      n = nlon*nlat
      
      if(nrow(newdata)!=n) stop("Predictions cannot be plotted in a grid") 
      pred = array(dim=c(nlon, nlat, length(pred.info$files)))
      
      lon = sort(unique(newdata$lon))
      lat = sort(unique(newdata$lat))
      
      coords = list(lon=lon, lat=lat)
      
    }
    
    if(!is.null(pred.info$aux)) {
      
      #       newdata$lat = round(newdata$lat, 2)
      #       newdata$lon = round(newdata$lon, 2)
      #       aux$lat = round(aux$lat, 2)
      #       aux$lon = round(aux$lon, 2)
      
      newdata = merge(newdata, aux, by=c("lon", "lat"), all=TRUE, sort=FALSE)
      
    }
    
    if(!is.null(object$transform)) newdata = FUN0(newdata)
    newdata = FUN(newdata)
    
    model.vars = .getModelVars(model) 
    ind = complete.cases(newdata[, model.vars])
    newdata = newdata[ind, ]
    
    if(nrow(newdata)==0) next
    
    temp = predict(model, type="response", newdata=newdata)
    temp = as.integer(round(factor*temp,3))
    pred[,,i][ind] = temp 
    
    # begin limiting variables
    limiting[i, ] = limitingFactor(data=newdata, ranges=values)
    # end limiting variables
  }
  
  xvar = which(names(newdata) %in% names(values))
  colnames(limiting) = names(newdata)[xvar]
  rownames(limiting) = round(pred.info$time$center, 2)
  
  
  if(!identical(dim(pred)[1:2], dim(mask)) & interp) {
    warning("'map' and 'mask' dimensions don't agree. Cannot perform interpolation")
    interp = FALSE
    mask   = NULL
  }
  
  if(interp) {
    DateStamp("Performing bilinear spatial interpolation over mask...")
    pred = fillMap(pred, mask=mask, radius=radius)
    DateStamp("Performing spline temporal interpolation...")
    pred = interpolateMap(pred, anomalies=FALSE)
  }
  
  thr = object$threshold$validation[,model.name]
  fit = object$fit[model.name,]
  per = rbind(training=object$performance$training[model.name,],
              validation=object$performance$validation[model.name,])
  fml = object$formulas[model.name]
  
  info = list(coords=coords, time=pred.info$time, factor=factor,
              threshold=thr, fit=fit, performance=per, formula=fml,
              model=model.name, transform=FUN, mask=mask)
  
  DateStamp("Computing climatologies...")
  pred.mean   = apply(pred, 1:2, median, na.rm=TRUE)
  pred.sd     = apply(pred, 1:2, sd, na.rm=TRUE)
  pred.clim   = climatology(pred, info$time$month)
  pred.season = climatology(pred, info$time$season)
  
  output = list(prediction=pred, info=info, mean=pred.mean, sd=pred.sd,
                climatology=pred.clim, season=pred.season, limiting=limiting)
  
  DateStamp("DONE.")
  
  class(output) = "prediction.niche.models"
  
  return(output)
}


window.prediction.niche.models = function(x, start=NULL, end=NULL, frequency=NULL, 
                                          deltat=NULL, extend=FALSE, ...) {
  
  if(is.null(start) & is.null(end)) return(x)
  
  year  = x$info$time$year
  month = x$info$time$month
  
  if(is.null(start)) start = c(year[1], month[1])
  if(is.null(end))   end   = c(tail(year, 1), tail(month, 1))
  
  oldtime = x$info$time$center
  newtime = createTimeAxis(start=start, end=end, frequency=12, center=TRUE)  
  date.min = min(newtime$center)
  date.max = max(newtime$center)
  
  time.ind = oldtime >= date.min & oldtime <= date.max 
  
  x$prediction = x$prediction[, , time.ind]
  x$info$time = newtime
  
  DateStamp("Computing climatologies for new time window...")
  x$mean        = apply(x$prediction, 1:2, median, na.rm=TRUE)
  x$sd          = apply(x$prediction, 1:2, sd, na.rm=TRUE)
  x$climatology = climatology(x$prediction, x$info$time$month)
  x$season      = climatology(x$prediction, x$info$time$season)
  
  return(x)
}

subset.prediction.niche.models = function(x, lat=NULL, lon=NULL) {
  
  if(is.null(lat) & is.null(lon)) return(x)
  
  xlat = x$info$coords$lat
  xlon = x$info$coords$lon
  
  if(is.null(lat)) lat = range(xlat)
  if(is.null(lon)) lon = range(xlon)
  
  lat = sort(lat)
  lon = sort(lon)
  
  indLat = xlat >= lat[1] & xlat <= lat[2]
  indLon = xlon >= lon[1] & xlon <= lon[2]
  
  x$prediction = x$prediction[indLon, indLat, ]
  x$info$coords$lon = xlon[indLon]
  x$info$coords$lat = xlat[indLat]
  
  DateStamp("Computing climatologies for spatial domain...")
  x$mean        = apply(x$prediction, 1:2, median, na.rm=TRUE)
  x$sd          = apply(x$prediction, 1:2, sd, na.rm=TRUE)
  x$climatology = climatology(x$prediction, x$info$time$month)
  x$season      = climatology(x$prediction, x$info$time$season)
  
  return(x)
}

getModel = function(object, model.name) {
  
  dataset = deparse(substitute(object))
  if(!(model.name %in% names(object$models))) stop("Model ", sQuote(model.name), " not found")  
  model = object$models[[model.name]]
  model$call$data    = eval(parse(text=paste0("get(\"", dataset,"\")$train")))
  model$call$formula = as.formula(.fmla2txt(model$formula))
  
  return(model)
}

getFormula = function(object, model.name) {
  
  if(!(model.name %in% names(object$models))) stop("Model ", sQuote(model.name), " not found")  
  fmla = object$formulas[[model.name]]
  
  return(fmla)
}

.getBestModel =function(object, criteria="AUC") {
  
  n = which.max(object$performance$validation[,criteria])
  model = rownames(object$performance$validation)[n]
  
  return(model)
}

getBestModel = function(object, criteria="AUC") {
  
  model = getModel(object=object, model=.getBestModel(object=object, 
                                                      criteria=criteria))
  return(model)
}

getThreshold = function(x, criteria="MinROCdist") {
  
  thr = x$info$threshold[criteria]
  return(thr)
}

getMap = function(object, ...) {
  
  UseMethod("getMap")
}

getMap.prediction.niche.models = function(object, date, toPA=FALSE, prob=FALSE, 
                                          criteria="MinROCdist", ...) {
  
  thr = getThreshold(x=object, criteria=criteria)
  dates = data.frame(year=object$info$time$year, month=object$info$time$month)
  if(!is.null(date)) {
    if(date[2]>12 | date[2] < 1 | date[2]%%1!=0) {
      stop("month has to be an integer between 1 and 12")
    }
    slice = which(dates$year==date[1] & dates$month==date[2])
    if(length(slice)==0) stop("date is out of range")
    z = object$prediction[,,slice]/object$info$factor
    
    if(toPA) z = toPA(z, thr, prob=prob)
    output = list(lon=object$info$coords$lon, lat=object$info$coords$lat, z=z)
    
    return(output)
  }
}

plot.prediction.niche.models = function(x, y=NULL, date=NULL, 
                                        slice=NULL, type=NULL, mfclim=c(3,4), 
                                        mar=c(4,4,1,1), oma=c(1,1,1,1),
                                        toPA=FALSE, prob=FALSE, criteria="MinROCdist", 
                                        png=FALSE, pdf=FALSE, filename=NULL, 
                                        width=800, height=800, res=NA, ...) {
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  thr = getThreshold(x=x, criteria=criteria)
  dates = data.frame(year=x$info$time$year, month=x$info$time$month)
  if(!is.null(date)) {
    if(date[2]>12 | date[2] < 1 | date[2]%%1!=0) {
      stop("month has to be an integer between 1 and 12")
    }
    slice = which(dates$year==date[1] & dates$month==date[2])
    if(length(slice)==0) stop("date is out of range")
    time.lab = paste(month.abb[date[2]],date[1], sep=" - ")
  } else {
    if(is.null(slice)) {
      time.lab = paste0("average ", min(dates$year), "-", max(dates$year))
    } else {
      time.lab = paste(month.abb[dates$month[slice]],
                       dates$year[slice], sep=" - ")
    }
  }
  if(is.null(date) & is.null(slice)) {
    z = x$mean/x$info$factor
    if(toPA) z = toPA(z, thr, prob=prob)
  } else {
    z = x$prediction[,,slice]/x$info$factor
    if(toPA) z = toPA(z, thr, prob=prob)
  }
  time.lab = toupper(time.lab)
  if(is.null(type)) {
    try(image.map(x$info$coords$lon, x$info$coords$lat, z, zlim=c(0,1), ...),
        silent=TRUE)
    points(-180,-180)
    mtext(time.lab, 3, line=1, adj=1, cex=0.75) 
    filename = .getPngFileName(filename, gsub(" ", "", time.lab))
  } else {
    if(type=="climatology") {
      par(mfrow=mfclim, mar=mar, oma=oma)
      z = x$climatology/x$info$factor
      if(toPA) z = toPA(z, thr, prob=prob)
      for(i in 1:12) {
        try(image.map(x$info$coords$lon, x$info$coords$lat, z[,,i], zlim=c(0,1),
                      ...),
            silent=TRUE)
        points(-180,-180)
        mtext(month.name[i], 3, line=1, adj=1, cex=0.75)            
      }
      filename = .getPngFileName(filename, "climatology")
    }
    
    if(type=="seasonal") {
      par(mfrow=c(2,2), mar=mar, oma=oma)
      z = x$season/x$info$factor
      if(toPA) z = toPA(z, thr, prob=prob)
      slab = toupper(dimnames(x$season)[[3]])
      for(i in 1:4) {
        try(image.map(x$info$coords$lon, x$info$coords$lat, z[,,i], zlim=c(0,1), 
                      ...),
            silent=TRUE)
        points(-180,-180)
        mtext(slab[i], 3, line=1, adj=1, cex=0.75)            
      }
      filename = .getPngFileName(filename, "seasonal")
    }
    
    if(!(type %in% c("seasonal", "climatology"))) {
      stop("Invalid 'type' argument, it can be 'seasonal' or 'climatology'.")
    }
  }
  
  if(isTRUE(png)) {
    png = match.fun("png")
    dev.copy(png, filename=filename, res=res, width=width, height=height)
    dev.off()
  }
  
  return(invisible())
}

saveAnimation.prediction.niche.models = function(object, file, dir=getwd(), interval=0.5, ...) {
  
  n = dim(object$prediction)[3]
  DateStamp("\nCreating animation (",n," time steps).", sep="")  
  
  try(suppressMessages(saveGIF(
    {
      for(i in seq_len(n)) {
        cat(i,"")
        plot(object, slice=i, ...)
      }
      },
    movie.name="temp.gif",
    img.name="niche_model", clean=TRUE, verbose=FALSE, 
    interval=interval, loop=1, check=TRUE, autobrowse=FALSE)),
    silent=TRUE)
  tmp = file.path("temp.gif")
  x = file.copy(from=tmp, to=file.path(dir, file), overwrite=TRUE)
  file.remove(tmp)
  DateStamp("DONE.")
  
  return(invisible(x))
}

climatology = function(object, index, FUN="median") {
  
  FUN = match.fun(FUN)
  dims = dim(object)
  ind  = unique(index)
  uind = length(ind)
  if(is.numeric(ind)) ind = sort(ind)
  out = array(0L, dim=c(dims[1:2], uind))
  for(i in seq_len(uind)) {
    out[,,i] = as.integer(apply(object[,, which(index == ind[i])], 1:2, FUN, na.rm=TRUE))
  }
  dimnames(out)[[3]] = ind
  
  return(out)
}

anova.niche.models = function(object, model=NULL, criteria="AUC", ...) {
  
  if(is.null(model)) model=.getBestModel(object=object, criteria=criteria)
  print(object$models[[model]]$anova)
  
  return(invisible(NULL))
}

plot.niche.models = function(x, model=NULL, 
                             mar=c(5,5,1,1), oma=c(0.5,1,4,1.5), ...) {
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  model.name = if(is.null(model)) .getBestModel(x) else model 
  dataset = deparse(substitute(x))
  if(model.name!="all") {
    
    if(!(model.name %in% names(x$models))) stop("Model ", sQuote(model.name), " not found")
    .plot.niche.models(x=x, model.name=model.name, dataset=dataset, mar=mar, oma=oma, ...)
    
  } else {
    
    models = names(x$models)
    
    for(model.name in models) {
      .plot.niche.models(x=x, model.name=model.name, dataset=dataset, mar=mar, oma=oma, ...)  
    }
    
  }
  
  return(invisible())
}

.plot.niche.models = function(x, model.name, dataset, 
                              mar=c(5,5,1,1), oma=c(0.5,1,4,1.5), ...) {
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  model = x$models[[model.name]]
  fml = .fmla2txt(model$formula)
  model$call$data    = eval(parse(text=paste0("get(\"", dataset,"\")$train")))
  model$call$formula = as.formula(fml)
  
  nvars = length(labels(model$terms))
  var = as.character(model$formula)[2]
  par(mfrow=getmfrow(nvars), mar=mar, oma=oma)
  
  plot(model, se=TRUE, ...)
  mtext(paste0(var, " (", dataset, ")"), 3, line=2, outer=TRUE)
  mtext(fml, 3, line=1, outer=TRUE)
  
}

getCoordinates = function(object, index) {
  
  if(is.logical(index)) index=which(index)
  if(!is.numeric(index)) stop("index must be logical or numerical")
  cols = col(object)[index]
  rows = row(object)[index]
  out = list(x=rows, y=cols)
  
  return(out)
}

map2coord = function(object) {
  
  out = data.frame(x = as.numeric(row(object)),
                   y = as.numeric(col(object)),
                   z = as.numeric(object)
                   )
  out = out[complete.cases(out),]
  
}

.getIndex = function(x, max, radius) {
  
  x = trunc(x)
  ind = seq(from=x-radius, to=x+radius, by=1)
  ind = ind[ind>0 & ind<=max]
  
  return(ind)
}

fillSquare = function(x, map, radius=1) {
  
  x = x[!complete.cases(x),]
  
  .fillSquare = function(x, map, radius) {
    ix = .getIndex(x[1], max=nrow(map), radius=radius)
    iy = .getIndex(x[2], max=ncol(map), radius=radius)
    z  = median(map[ix, iy], na.rm=TRUE)
    
    return(z)
  }
  
  out = apply(x, 1, .fillSquare, map=map, radius=radius)
  
  return(out)
}

fillMap = function(object, mask, radius=1, ...) {
  
  UseMethod("fillMap")
}

fillMap.matrix = function(object, mask, radius=1, fill.value=0) {
  
  allNA  = all(is.na(object))
  noMiss = !any(is.na(object)&!is.na(mask)) 
  if(!allNA & !noMiss) {
    if(!identical(dim(object), dim(mask))) stop("'object' and 'mask' dimension must agree.")
    coords = map2coord(object)
    ind    = is.na(object)&!is.na(mask)
    miss   = getCoordinates(object, ind)
    fill   = interpp(coords$x, coords$y, coords$z, miss$x, miss$y)
    map2   = object
    map2[ind] = fill$z
    fill   = as.data.frame(fill)
    fill$z[is.na(fill$z)] = fillSquare(x=fill, map=map2, radius=radius)
    fill$z[is.na(fill$z)] = fill.value # just in case!
    map2[ind] = fill$z
    map2[is.na(mask)] = NA
  } else {
    return(object)
  }
  
  return(map2)
}

fillMap.data.frame = function(object, mask=NULL, var, radius=1, dx=1/24, thr=10) {
  
  z = object[, var]
  latR = range(pretty(object$lat))
  lonR = range(pretty(object$lon))
  axes = createGridAxes(lat=latR, lon=lonR, dx=dx)
  fill = interpp(x=object$lat, y=object$lon, z=object[, var], 
                 xo=as.numeric(axes$LAT), yo=as.numeric(axes$LON))
  #   map = matrix(fill$z, nrow=nrow(axes$LON), ncol=ncol(axes$LON))
  #   mask = array(1, dim=dim(axes$LON))
  #   map = fillMap(map, mask, radius=2)
  #   fill$z = as.numeric(map)
  fill$z[fill$z<thr] = NA
  names(fill) = c("lat", "lon", var)
  fill = as.data.frame(fill)
  
  return(fill)
}

fillMap.array = function(object, mask, radius=1) {
  
  dmap = dim(object)
  out = apply(object, 3, fillMap, mask=mask, radius=radius)
  dim(out) = dmap
  
  return(out)
}

fillMap.prediction.niche.models = function(object, mask, radius=1) {
  
  if(!identical(dim(object$prediction)[1:2], dim(mask))) 
    stop("'object' and 'mask' dimension must agree.")
  
  DateStamp("Performing spatial interpolation over mask...")
  object$prediction = fillMap(object$prediction, mask=mask, radius=radius)
  object$info$mask = mask
  
  DateStamp("Performing temporal interpolation...")
  object$prediction = interpolateMap(object$prediction, anomalies=FALSE)
  
  DateStamp("Recalculating climatologies...")
  object$mean         = apply(object$prediction, 1:2, median, na.rm=TRUE)
  object$sd           = apply(object$prediction, 1:2, sd, na.rm=TRUE)
  object$climatology  = climatology(object$prediction, object$info$time$month)
  object$season       = climatology(object$prediction, object$info$time$season)
  
  DateStamp("DONE.")
  
  return(object)
}

.interpolateMap = function(x, anomalies=FALSE, ...) {
  
  if(all(is.na(x)) | all(!is.na(x))) return(x)
  
  out = x
  t = seq_along(x)
  
  if(anomalies) {
    xclim = tapply(x, INDEX=rep(1:12, len=length(x)), FUN=mean, na.rm=TRUE)
    xclim = rep(xclim, len=length(out))
    out = x - xclim
  }
  
  xs = splinefun(x=t, y=out)
  out[is.na(out)] = pmin(pmax(xs(t[is.na(out)]),min(out, na.rm=TRUE)),
                         max(out, na.rm=TRUE))
  
  if(anomalies) out = out + xclim
  
  return(out)
}

interpolateMap = function(object, anomalies=FALSE, ...) {
  
  pred = apply(object, 1:2, .interpolateMap, anomalies=anomalies, ...)
  pred = aperm(pred, c(2,3,1))
  
  return(pred)
}

plotDistributionMap <- function (object, var, dx, FUN=mean, log=FALSE, 
                                 hires=FALSE, thr=-Inf, ...) {
  # to be improved: interpolation
  FUN = match.fun(FUN)
  
  latR = range(pretty(object$lat))
  lonR = range(pretty(object$lon))
  
  axes = createGridAxes(lat=latR, lon=lonR, dx=dx)
  
  fill = fillMap(object, var=var, dx=dx/4)
  map  = table2grid(fill, var=var, 
                    lat=latR, lon=lonR, dx=dx, FUN=FUN)
  if(isTRUE(log)) map = log(map)
  map[map<thr] = NA
  image.map(axes$lon, axes$lat, map, hires=hires, ...)
  
  return(invisible())
}

