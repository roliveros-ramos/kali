fitGAM2 = function(object, formulas, FUN=identity, refit=FALSE, 
                   link="logit", select=NULL, bigData=FALSE, cluster=NULL) {
  
  name = deparse(substitute(object))
  
  if(!is.null(cluster)) bigData = TRUE

  FUN = match.fun(FUN)
  object$transform = FUN
  
  gc(verbose=FALSE)
  
  DateStamp("Fitting models for", sQuote(name), "dataset.")
  
  train = FUN(object$train)
  val   = FUN(object$val)
  
  if(!is.list(object$predicted)) object$predicted = NULL
  
  formulas = lapply(formulas, as.formula)

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
    
    model.formula = formulas[[i]]
    DateStamp("Training model ", model.name, ":\n", .fmla2txt(formulas[[i]]), sep="")
    object$formulas[[model.name]] = model.formula 
    #model.vars = .getModelVars2(model.formula, train)
    # TO_DO: filter complete cases
    if(isTRUE(bigData)) {
      model = mgcv::bam(model.formula, data = train, family = binomial(link=link),
                        select=select[model.name], method="fREML", cluster=cluster)
    } else {
      model = mgcv::gam(model.formula, data = train, family = binomial(link=link),
                        select=select[model.name], method="REML")
    }
    
    gc(verbose=FALSE)
    # model$anova = anova(model)
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
