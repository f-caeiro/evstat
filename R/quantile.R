# Estimator of the Scale Second Order Parameter (beta)
xquantile <- function(x, log.x=FALSE, method="weissman", evi="hill", threshold=NA, p=0.001, shape2=NA, tau=0, ...){  
  est <- NULL
  stopifnot(is.numeric(x))
  #thr <- threshold
  #
  n0 <- length(x)
  #if (is.na(shape2)){shape2 <- adapt.shape2(x, log.x=log.x)}
  if (!log.x){
    lx <- log(x[x>0])
    log.x <- TRUE
  }
  else {lx <- x}
  #
  x <- sort(x)
  #rho <- shape2
  if (method == "weissman"){
    if (!log.x) lx <- log(x[x>0])
    x <- x[x>0]
    n <- length(x)
    evi <- evindex(x, method=evi, ...)
    est <- ((1:(n-1))/(n0*p))^evi*x[(n-1):1]
  }
  if (method == "plpwm"){
    if (log.x==FALSE) lx <- log(x[x>0])
    x <- x[x>0]
    n <- length(x)
    evi <- evindex(x, method="plpwm")
    est <- ((1:(n-1))/(n0*p))^evi*x[(n-1):1]
  }  
  ifelse(is.na(threshold), return(est), return(est[threshold]))
  #return(temp)
}


xprob <- function(x, log.x=FALSE, method="weissman", evi="hill", threshold=NA, level=200, nby=12, shape2=NA, tau=0){  
  est<-NULL
  stopifnot(is.numeric(x))
  #thr <- threshold
  #
  n0 <- length(x)
  if (is.na(shape2)){shape2 <- adapt.shape2(x, log.x=log.x)}
  if (log.x==FALSE){
    lx <- log(x[x>0])
    log.x<-T
  }
  else {lx <- x}
  #
  x<-sort(x)
  rho <- shape2
  if (method == "weissman"){
    if (log.x==FALSE) lx <- log(x[x>0])
    x <- x[x>0]
    n <- length(x)
    evi <- evindex(x, method=evi)
    est <- ((1:(n-1))/(n0))*(level/x[(n-1):1])^(-1/evi)
    est <- 1-(1-est)^nby
  }
  if (method == "plpwm"){
    if (log.x==FALSE) lx <- log(x[x>0])
    x <- x[x>0]
    n <- length(x)
    evi <- evindex(x, method="plpwm")
    est <- ((1:(n-1))/(n0))^evi*x[(n-1):1]
  }  
  ifelse(is.na(threshold), return(est), return(est[threshold]))
  #return(temp)
}


