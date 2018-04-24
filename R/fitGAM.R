fitGAMs = function(object, formulas, FUN=identity, refit=FALSE, 
                   link="logit", select=NULL, bigData=FALSE, cluster=NULL,
                   useWeights=FALSE, method=NULL) {
  
  name = deparse(substitute(object))
  
  if(is.null(method) & isTRUE(bigData))  method = "fREML" 
  if(is.null(method) & !isTRUE(bigData)) method = "REML" 
    
  if(!is.null(cluster)) bigData = TRUE
  if(nrow(object$train)>1e5) bigData = TRUE
  
  FUN = match.fun(FUN)
  object$transform = FUN
  
  gc(verbose=FALSE)
  
  DateStamp("Fitting models for", sQuote(name), "dataset.")
  
  train = FUN(object$train)
  val   = FUN(object$val)
  
  if(!is.list(object$predicted)) object$predicted = NULL
  
  for(i in seq_along(formulas)) {
    formulas[[i]] = as.formula(formulas[[i]])
    environment(formulas[[i]]) = environment()
  }

  models = names(formulas)
  if(is.null(select)) select = rep(FALSE, length(models))
  if(length(select)==1) select = rep(select, length(models))
  if(length(select)!=length(models)) 
    stop("you have to specify one sQuote(select) per model.")
  names(select) = models
  
  aic = object$fit[,"AIC"]
  bic = object$fit[,"BIC"]
  
  var = as.character(formulas[[1]])[2]
  
  object$predicted$lon = object$train[, "lon"]
  object$predicted$lat = object$train[, "lat"]
  object$predicted$observed = as.numeric(as.character(object$train[,var]))
  object$predicted = as.data.frame(object$predicted)
  
  object$validation$lon = object$val[, "lon"]
  object$validation$lat = object$val[, "lat"]
  object$validation$observed = as.numeric(as.character(object$val[,var]))
  object$validation = as.data.frame(object$validation)
  
  for(i in seq_along(formulas)) {
    model.name = models[i]
    
    if((model.name %in% names(object$models)) & !refit) {
      DateStamp("Skipping ", model.name, ":already fitted.\n", sep="")
      next
    }
    
    DateStamp("Training model ", model.name, ":\n", .fmla2txt(formulas[[i]]), sep="")
    
    model.formula = formulas[[i]]
    
    ww = if(isTRUE(useWeights)) {
      calculatePointWeight(object$weights, model.formula)
    } else NULL
    
    object$formulas[[model.name]] = model.formula 
    #model.vars = .getModelVars2(model.formula, train)
    # TO_DO: filter complete cases
    if(isTRUE(bigData)) {
      model = mgcv::bam(model.formula, data = train, family = binomial(link=link),
                        select=select[model.name], method=method, cluster=cluster,
                        weights=ww)
    } else {
      model = mgcv::gam(model.formula, data = train, family = binomial(link=link),
                        select=select[model.name], method=method, weights=ww)
    }
    
    gc(verbose=FALSE)
    # model$anova = anova(model)
    model$call$family[2] = link
    object$models[[model.name]] = model
    # object$preplot[[models[i]]] = with(object$train, preplot(model))
    
    aic[model.name] = AIC(model)
    bic[model.name] = BIC(model)
    
    if(!is.null(model$na.action)) {
      fitted.values = rep(NA, length(model$na.action) + length(model$fitted.values))
      fitted.values[-model$na.action] = model$fitted.values 
    } else fitted.values = model$fitted.values
    
    # object$predicted[, model.name] = predict(model, newdata=train, type="response")
    object$predicted[, model.name] = fitted.values
    object$validation[, model.name] = predict(model, newdata=val, type="response", 
                                              cluster=cluster)
    
  }
  
  DateStamp("Computing Predictive Performance...")
  
  object$fit = cbind(AIC=aic, BIC=bic)
  
  object$performance$training   = PredictivePerformance(object$predicted, st.dev=FALSE)
  object$performance$validation = PredictivePerformance(object$validation, st.dev=FALSE)
  
  object$threshold$training  = calculateThresholds(object$predicted)
  object$threshold$validation = calculateThresholds(object$validation)
  
  DateStamp("DONE.")
 
  object$species = as.character(object$formulas[[1]][[2]])
   
  class(object$formulas) = c("niche.models.formulas", class(object$formulas))
  
  return(object)
}


# Data weighting ----------------------------------------------------------

.calculatePointWeight = function(x, n=100) {
  
  cuts = xrange(x, n=n+1)
  y = findInterval(x, cuts)
  tab  = table(y)
  tab2 = rep(0, n)
  tab2[as.numeric(names(tab))] = tab
  z = tab2[y]
  z = 1/z
  z = z/mean(z)
  return(z)
}

calculatePointWeight = function(w, formula) {
  wnames = names(w)
  if(is.null(wnames)) stop("Matrix of weights must have column names.")
  xvars = .getUsedVars(x=wnames, formula=formula)
  sw = w[, xvars, drop=FALSE]
  ww = apply(sw, 1, prod)
  ww = ww/mean(ww)
  return(ww)
}

.getUsedVars = function(formula, x) {
  xvars = paste(as.character(formula[[3]]), collapse=" + ")
  xvars = names(which(sapply(sapply(x, grep, x=xvars), length)>0))
  return(xvars)
}

# 
# rank = function(x) {
#   rk = sort(x, index.return=TRUE, decreasing = TRUE)$ix
#   out = (seq_along(x))[sort(rk, index.return=TRUE)$ix]
#   return(out)
# }
# 

performance = function(object, type="validation", rank=FALSE, 
                       aggregate=FALSE, sort=FALSE) {
  
  if(!inherits(object, "niche.models")) stop("object is not of class 'niche.models'.")
  if(isTRUE(aggregate)) rank = TRUE
  
  out = object$performance[[type]]
  
  if(isTRUE(rank)) {
    out = sapply(1-out, rank, ties="min")
    rownames(out) = names(object$models)
  }
  if(isTRUE(aggregate)) {
    out = rowSums(out)
    if(isTRUE(sort)) {
      out = sort(out)
      out = names(out) 
      }
    }
  
  return(out)
  
}

