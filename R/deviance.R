
#' @export
explainedDeviance =  function(object, fixed, null, newdata, obs, ...) {
  UseMethod("explainedDeviance")  
}

#' @export
explainedDeviance.gam = function(object, fixed=NULL, null=NULL, newdata=NULL, obs=NULL, ...) {
  .explainedDeviance(object=object, fixed=fixed, null=null, newdata=newdata, obs=obs, isGAM=TRUE, ...)
}

#' @export
explainedDeviance.scam = function(object, fixed=NULL, null=NULL, newdata=NULL, obs=NULL, ...) {
  .explainedDeviance(object=object, fixed=fixed, null=null, newdata=newdata, obs=obs, isGAM=TRUE, ...)
}

#' @export
explainedDeviance.glm = function(object, fixed=NULL, null=NULL, newdata=NULL, obs=NULL, ...) {
  .explainedDeviance(object=object, fixed=fixed, null=null, newdata=newdata, obs=obs, isGAM=FALSE, ...)
}


.explainedDeviance = function(object, fixed=NULL, null=NULL, newdata=NULL, obs=NULL, isGAM=FALSE, ...) {
  
  cat(date(), "\n")
  
  if(!is.null(newdata) & is.null(obs)) stop("'obs' must be provided to use newdata.")
  
  if(is.null(null)) null = object
  if(inherits(null, "glm")) null = null$null.deviance
  
  term = stats::terms(object$formula, keep.order=TRUE, specials=c("s", "te", "ti", "t2"))
  sn  = attr(term,"term.labels")
  special = sort(unlist(attr(term, "specials")) - 1)
  if(length(special)==0) isGAM = FALSE
  
  if(isTRUE(isGAM)) {
    gp = mgcv::interpret.gam(gf = object$formula)
    sn0 = sapply(gp$smooth.spec, FUN="[[", i="label")
    sn[special] = sn0 # nice names
    # sp for smooths
    .getDim = function(i) length(grep(x=names(object$sp), pattern=sn0[i], fixed = TRUE))
    sl  = sapply(seq_along(sn0), FUN=.getDim) # work with bs="ad"
    # sl  = sapply(gp$smooth.spec, FUN="[[", i="dim")
    st  = cumsum(c(1, head(sl,-1)))
    en  = cumsum(sl)
    sps = setNames(mapply(seq, from=st, to=en), sn0)
  }
  
  ch  = setNames(lapply(sn, FUN=function(...) c(TRUE, FALSE)), sn)
  ch[fixed] = TRUE
  mm = do.call(expand.grid, ch) # mm: model matrix
  mm = mm[apply(mm, 1, any), ] # should we remove the constant model?
  
  .update = function(index, term) nfmla = stats::formula(term[index])
  .sp = function(index, sp, pos) sp[unlist(pos[unlist(index)])]
  .mexplDev = function(object, null) (null - object$deviance)/null
  
  fmlas = apply(mm, 1, FUN=.update, term=term)
  
  if(isTRUE(isGAM)) xsp = apply(mm[special], 1, FUN=.sp, sp=object$sp, pos=sps)
  
  xcall = object$call
  if(!is.null(xcall$weights)) xcall$weights = as.name('(weights)')
  
  # models = list()
  
  if(!is.null(newdata)) preds = matrix(NA_real_, ncol=length(fmlas), nrow=nrow(newdata))
  devs = numeric(length(fmlas))
  
  if(!is.null(newdata)) {
    preds[, 1] = predict(object, type="response", newdata=newdata)
  } else devs[1] = .mexplDev(object, null=null)
  
  for(i in seq_along(fmlas)[-1]) {
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, i/length(fmlas))
    icall = xcall
    icall$formula = stats::as.formula(fmlas[[i]], env=environment())
    icall$data  = quote(object$model)
    if(isTRUE(isGAM)) {
      isp = if(length(xsp[[i]])!=0) xsp[[i]] else NULL
      icall$sp = quote(isp)
    } 
    modelTmp = eval(icall)
    if(!is.null(newdata)) {
      preds[, i] = predict(modelTmp, type="response", newdata=newdata)
    } else devs[i] = .mexplDev(modelTmp, null=null)
  }
  
  mm2 = as.data.frame(sapply(mm, as.numeric))
  
  if(!is.null(newdata)) {
    .getDev = function(mu) object$family$dev.resids(y=obs, mu=mu, wt = 1)
    mm2$dev = colSums(apply(preds, 2, FUN=.getDev)) 
  } else mm2$dev = devs #sapply(models, FUN=.mexplDev, null=null)
  
  if(!is.null(fixed)) {
    mm2[, fixed] = NULL
    mm2 = cbind("_fixed_terms"=1, mm2)
  }
  
  edm = stats::glm(dev ~ . + 0 , data=mm2) 
  ed = stats::coef(edm)/sum(stats::coef(edm), na.rm=TRUE)  
  
  if(!is.null(fixed)) attr(ed, "_fixed_terms") = names(mm)[fixed]
  
  class(ed) = c("kali_explainedDeviance")
  
  pb = txtProgressBar(style=3)
  setTxtProgressBar(pb, 1)
  cat("\n")
  cat(date(), "\n")
  
  return(ed)
  
}

print.kali_explainedDeviance = function(x, digits=2, ...) {
  
  msg = sprintf("%%0.%df%%%%", digits)
  print(setNames(sprintf(msg, 100*x), nm = names(x)))
  
}

# Explained deviance ------------------------------------------------------

#' @export
explDev = function(..., null=NULL, full=FALSE, sort=TRUE) {
  
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  if(inherits(null, "glm")) null = null$deviance
  
  .explDev = function(object, null=NULL) {
    if(is.null(null)) null = object$null.deviance 
    (null - object$deviance)/null
  }
  
  .explDev2 = function(object, null=NULL) {
    if(is.null(null)) null = object$null.deviance 
    int = coef(object)["(Intercept)"]
    return(c(int, null, object$deviance, (null - object$deviance)/null))
  }
  
  if(!isTRUE(full)) {
    out = sapply(list(...), FUN=.explDev, null=null)
    names(out) = names
    if(isTRUE(sort)) out = sort(out, decreasing = TRUE)
  } else {
    out = t(sapply(list(...), FUN=.explDev2, null=null))
    rownames(out) = names
    colnames(out) = c("Intercept", "null", "deviance", "expl.deviance")
    out = as.data.frame(out)
    if(isTRUE(sort)) out = out[order(out$expl.deviance), ]
  }
  
  return(out)
  
}
