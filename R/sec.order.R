#
# Estimator of the Shape Second Order Parameter (rho)
#
shape2 <- function(x, log.x=FALSE, method="fagh", threshold=NA, tau=0, adjust="neg.abs"){
  est <- NULL 
  stopifnot(is.numeric(x))
  #
  if (log.x==FALSE){x <- log(x[x>0])}
  n <- length(x)
  #
  # Estimator in Hall and Welsh (1985)
  if (method=="hw"){    
    xlog <- rev(sort(x))
    tmp <- cumsum(xlog)/(1:n)
    t <- 1/(tmp[-n]-xlog[-1])
    t0 <- floor(t[n^0.5])
    t1 <- floor(t[n^0.9])
    t2 <- floor(t[n^0.95])
    t <- (t1-t0)/(t2-t0)
    t <- abs(t)/log(n^0.9/n^0.95)
    est <- -abs(t)
  }  
  # Estimator in Fraga Alves et al. (2003)
  if (method=="fagh"){
    xlog1 <- rev(sort(x))
    xlog2 <- xlog1^2
    xlog3 <- xlog1^3
    t1 <- cumsum(xlog1)/(1:n)
    t2 <- cumsum(xlog2)/(1:n)
    t3 <- cumsum(xlog3)/(1:n)
    m1 <- t1[-n]-xlog1[-1]
    m2 <- (t2[-n]-(2*t1[-n]*xlog1[-1])+xlog2[-1])/2
    m3 <- (t3[-n]-(3*t2[-n]*xlog1[-1])+(3*t1[-n]*xlog2[-1])-(xlog3[-1]))/6
    
    ifelse(tau==0, temp<-(log(m1)-log(m2)/2)/(log(m2)/2-log(m3)/3), 
           temp<-(m1^tau-m2^(tau/2))/(m2^(tau/2)-m3^(tau/3)))
    est <- 3*(temp-1)/(temp-3)
  }
  #
  # Estimator in Caeiro and Gomes (2014)
  if (method=="cg"){  
    xlog <- rev(sort(x))
    ui<-(xlog[-n]-xlog[-1])
    m1 <- cumsum((1:(n-1))*ui)/(1:(n-1))
    m2 <- cumsum((1:(n-1))^(1.5)*ui)/((1:(n-1))^(1.5))*1.5
    m3 <- cumsum((1:(n-1))^2*ui)/((1:(n-1))^2)*2
    ifelse(tau==0, temp<-(log(m1)-log(m2))/(log(m2)-log(m3)),
           temp<-(m1^tau-m2^tau)/(m2^tau-m3^tau))
    est <- (2-temp)/(1-temp)
  }
  if (adjust=="neg.abs"){est<- -abs(est)}
  if (adjust=="min"){est<-ifelse(est>0, 0, est)}
  #if (est < -7) est <- -1
  ifelse(is.na(threshold), return(est), return(est[threshold]))
}


adapt.shape2 <- function(x, log.x=FALSE, shape2=c("fagh", "cg"), threshold=floor(length(x[x>0])^0.999), range=floor(length(x[x>0])^(c(0.995,0.999))), param=c(0,1), adjust="neg.abs"){
  shape2 = match.arg(shape2)
  for (i in 1:length(param)){
    x.r <- shape2(x, log.x=log.x, method=shape2, tau=param[i], adjust="neg.abs")
    x.n1 <- length(x.r)+1
    x1 <- range[1]
    x2 <- range[2]
    if (x1==x2){x1<-(x2-1)}    
    i1<-sum((x.r[x1:x2]-median(x.r[x1:x2]))^2)
    if (i==1){
      i0 <- i1
      est <- x.r[threshold]
    }
    if (i1<i0){
      i0 <- i1
      est <- x.r[threshold]
    }
  }
  return(est)
}


# Estimator of the Scale Second Order Parameter (beta)
scale2 <- function(x, log.x=FALSE, method="gm2002", threshold=floor(length(x[x>0])^(0.999)), shape2=NA, tau=0){  
  est<-NULL
  stopifnot(is.numeric(x))
  #
  if (is.na(shape2)){shape2 <- adapt.shape2(x, log.x=log.x)}
  if (log.x==FALSE){
    lx <- log(x[x>0])
    log.x<-T
  }
  else {lx <- x}
  #  
  rho <- shape2
  if (method == "gm2002"){
    lx <- sort(lx)
    n <- length(lx)
    xlog <- rev(lx)
    ui <- (xlog[-n]-xlog[-1])
    t0 <- cumsum((1:(n-1))*ui)/(1:(n-1))
    t1 <- cumsum((1:(n-1))^(1-rho)*ui)/((1:(n-1))^(1-rho))
    t2 <- cumsum((1:(n-1))^(1-2*rho)*ui)/((1:(n-1))^(1-2*rho))
    t3 <- (cumsum((1:(n-1))^(-rho))/((1:(n-1))^(-rho)))/(1:(n-1))    
    est <- (t3*t0-t1)/(t3*t1-t2)
    est <- est*((1:(n-1))/n)^(rho)    
  }
  ifelse(is.na(threshold), return(est), return(est[threshold]))
  #return(temp)
}


# Find the mode of a vector
sample.mode <- function(x) {
  un <- unique(x)
  m <- un[which.max(tabulate(match(x, un)))]
  c <- length(x[x==m])  
  return(list(mode=m,len.mode=c,max.mode=max(which(x==m))))
}

# Find the longest run of a vector rounded to dp decimal places
lrun <- function(x, dp=1){
  lgt.max<-0
  #for (i in 0:1){} x.runs<-rle(round(x+i*5*10^(-dp-1),dp))
  x.runs<-rle(round(x,dp))
  r.max<-which.max(x.runs$lengths)
  l.max<-x.runs$lengths[r.max]
  k1.max<-1
  if (r.max>1){k1.max<-k1.max+sum(x.runs$lengths[1:(r.max-1)])}
  s1<-k1.max
  s2<-k1.max+l.max-1  
  return(c(s1,s2))        
}  




mean.le <- function(x, p=1, log.x=TRUE){
  if (log.x == FALSE) x <- log(x[x>0])
  n <- length(x)
  xlog1 <- rev(sort(x))
  mm <- matrix(data=0, nrow=length(p), ncol=n-1)
  # Compute mean.le for every integer p
  if (length(which(p==round(p)))>0){
    idx <- which(p==round(p))
    a <- max(p[idx])
    m <- matrix(data=NA, nrow=a, ncol=n)
    for (i in 1:a) m[i,] <- cumsum(xlog1^i)/(1:n)
    for (i in 1:length(idx)){
      mm[idx[i],] <- ((-xlog1)^p[idx[i]])[-1]
      if (p[idx[i]] > 0) {
        for (j in 1:p[idx[i]]) mm[idx[i],] <- mm[idx[i],] + choose(p[idx[i]],j)*(m[j,-n])*((-xlog1)^(p[idx[i]]-j))[-1]
      }
    }
  }
  # Compute mean.le for every non integer p
  if (length(which(p!=round(p)))>0){
    idx <- which(p!=round(p))
    for (i in 1:length(idx)){
      for (j in 1:(n-1)){
        mm[idx[i],j] <- sum((xlog1[1:j]-xlog1[j+1])^(p[idx[i]]))/j
      }
    } 
  }
  return(mm)
}  

