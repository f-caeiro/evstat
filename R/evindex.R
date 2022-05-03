#
# Estimators of the Extreme Value Index
#
evindex <- function(x, log.x=FALSE, method="hill", threshold=NA, shape2=NA, scale2=NA, param=NA){
  x <- sort(x, decreasing=TRUE)
  # estimation of the second order parameters, if necessary
  rhobeta <- c("chill", "chill", "chill2", "wh", "NCH", "PRB", "PRBstar", "CHp", "CHstar", "ML", "MLb", "NML")
  if (is.element(method, rhobeta)){
  #if (method=="chill" | method=="chill1"| method=="chill2" | method=="PRB" | method=="PRBstar" | method=="CHp" | method=="CHstar" | method=="ML" | method=="MLb" | method=="NML"){
    if (is.na(shape2)){shape2 <- adapt.shape2(x, log.x=log.x)}
    if (is.na(scale2)){scale2 <- scale2(x, shape2=shape2, log.x=log.x)}
    x.rho <- shape2
    x.beta <- scale2
    if (is.nan(x.beta)) x.beta <- 1
    #if(!is.numeric(x.rho)) print(c(x.rho,x.beta))
    #print(c(x.rho,x.beta))
  }
  # Hill estimator
  if (is.element(method, c("hill", "chill", "chill2", "NCH"))){
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    tmp <- cumsum(x)/(1:n)
    est <- tmp[-n]-x[-1]
  }
  if (method=="hill-lg"){#para o log gamma - não funciona
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    t1 <- cumsum(x)/(1:n)
    #est <- (t1[-n]-x[-1])*(1-0.2/log(n/(1:(n-1))))
    est <- (t1[-n]-x[-1])*exp(-.2/log(n/(1:(n-1))))
  }
  # Weighted Hill estimator
  if (method=="whill"){
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    t1 <- cumsum(x)/(1:n)
    t2 <- cumsum((1:n)*x)/((1:n)*2:(n+1))
    m1 <- t1[-n]-x[-1]
    m2 <- t2[-n]-0.5*x[-1] 
    est <- (1+param[1])*m1-4*param[1]*m2
  }
  # Weighted Hill estimator (with a selected to reduce bias)
  if (method=="whillrb"){
    if (is.na(shape2)){shape2 <- adapt.shape2(x, log.x=log.x)}
    a <- 0
    if (shape2<0) param[1] <- -(2-shape2)/shape2
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    t1 <- cumsum(x)/(1:n)
    t2 <- cumsum((1:n)*x)/((1:n)*2:(n+1))
    m1 <- t1[-n]-x[-1]
    m2 <- t2[-n]-0.5*x[-1] 
    est <- (1+param[1])*m1-4*param[1]*m2
  }    
  # Weighted Hill estimator-notoptimized
  if (method=="whill2"){
    if (log.x==FALSE) x <- log(x[x>0])
    n<-length(x)
    #cs1 <- cumsum(xlog)/(1:n)
    #cs2 <- cumsum((1:n)*xlog)/(1:n)
    #est <- (1+param[1])*(cs1[-n]-xlog[-1])-2*param[1]*(cs2[-n]-(1+(1:(n-1)))*xlog[-1]/(2*(1:n)))
    est <- rep(0,n-1)
    for (k in 1:(n-1)){
      t1<-0
      for (i in 1:k) {
        t1 <- t1+(1+param[1]*(1-4*i/k))*(xlog[i]-xlog[k+1])
        #t1 <- t1+(xlog[i]-xlog[k+1])
      }
      est[k] <- t1/k
    }
  }
  if (method=="whill-ls"){
    if (log.x==FALSE){x<-log(x[x>0])}  
    n <- length(x)
    xlog <- sort(x, decreasing = T)
    ui <- (xlog[-n]-xlog[-1])
    m1 <- cumsum((1:(n-1))*ui)/(1:(n-1))
    m2 <- cumsum((1:(n-1))^2*ui)/((1:(n-1))^2)
    est <- (1+param[1])*m1-2*param[1]*m2
  }
  # Corrected Hill estimator
  if (method=="chill"){
    m1 <- (x.beta/(1-x.rho))*(n/(1:(n-1)))^(x.rho)  
    est <- est*(1-m1)
  }
  # Corrected Hill estimator
  if (method=="chill2"){
    m1 <- (x.beta/(1-x.rho))*(n/(1:(n-1)))^(x.rho)  
    est <- est*exp(-m1)
  }  
  # Corrected Hill estimator2
  if (method=="NCH"){
    m1 <- (x.beta/(1-x.rho))*(n/(1:(n-1)))^(x.rho)  
    est <- est*(2-exp(m1))
  }
  # Corrected WH Hill estimator (Gomes et al. 2008)
  if (method=="wh"){
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    est <- rep(NA,n-1)
    for (k in 1:(n-1)){
      if (k>1){
        ik <- (1:(k-1))/k
        wik <- (1-ik^(-x.rho))/(x.rho*log(ik))
        wik <- c(wik,1)
      } else {wik <- 1}
      est[k] <- sum(exp(-x.beta*(n/k)^(x.rho)*wik)*(x[1:k]-x[1+k]))/k
    } 
  }
  # Lehmer Hill estimator
  if (method=="lehmer"){
    if (log.x==FALSE) x<-log(x[x>0])
    p <- param[1]
    tt <- mean.le(x, p=c(p-1,p), log.x=T)
    est <- tt[2,]/(p*tt[1,])
    #n <- length(x)
    #xx <- rev(x)
    #est <- rep(NA,n-1)
    #for (i in 1:(n-1)){
    #  h1 <- sum((xx[1:i]-xx[i+1])^(p))/i
    #  h2 <- sum((xx[1:i]-xx[i+1])^(p-1))/i
    #  est[i] <- h1/(p*h2)
    #}  
  }
  # Generalized Hill estimator in Gomes and Martins (2001)
  if (method=="ghill1"){
    if (log.x==FALSE) x<-log(x[x>0])
    p <- param[1]
    tt <- mean.le(x, p=c(p-1,p), log.x=T)
    est <- tt[2,]/(p*tt[1,])  
  }    
  # Shifted Hill estimator
  if (method=="shill"){
    n <- length(x)
    #x <- rev(x)
    est <- rep(NA,n-1)
    for (i in 1:(n-1)){
      s <- (5/16)*(i/n)^(3/4)
      est[i] <- sum(log((x[1:i]+s)/(x[i+1]+s)))/i
    }
  }
  # ppwm - there is other
  if (method=="ppwmtest"){
    n <- length(x)
    x <- sort(x, decreasing = TRUE)
    h1 <- (cumsum(x)/(1:n))[-1]
    h2 <- (cumsum((0:(n-1))*x))/(1:n)
    h2 <- h2[-1]/(1:(n-1))
    est <- 1-h2/(h1-h2)    
  }
  # MOP
  if (method=="mop"){
    n <- length(x[x>0])
    pp <- param[1]
    #phi <- (1-x.rho/2)-sqrt((1-x.rho/2)^2-0.5)    
    #m1 <- ((x.beta*(1-phi))/(1-x.rho-phi))*(n/(1:(n-1)))^(x.rho)   
    if (pp == 0){
      x <- log(x[x>0]) 
      #xlog<-rev(x)
      tmp <- cumsum(x)/(1:n)  
      est <- tmp[-n]-x[-1]
    }
    else {
      x <- (x[x>0])^pp 
      #x<-x0#(rev(x0))
      tmp<-(cumsum(x)/(1:n))
      est <- (1-(tmp[-n]/x[-1])^(-1))/pp      
    }   
  }
  # Partially reduced bias (previosus rbmop)
  if (method=="PRB"){
    n <- length(x[x>0])
    pp <- param[1]
    phi <- (1-x.rho/2)-sqrt((1-x.rho/2)^2-0.5)    
    m1 <- ((x.beta*(1-phi))/(1-x.rho-phi))*(n/(1:(n-1)))^(x.rho)   
    if (pp == 0){
      x<-log(x[x>0]) 
      #xlog<-rev(x)
      tmp<-cumsum(xlog)/(1:n)  
      tmp <- (tmp[-n]-xlog[-1])
    }
    else {
      x0<-(x[x>0])^pp 
      x<-x0#(rev(x0))
      tmp<-(cumsum(x)/(1:n))
      tmp <- (1-(tmp[-n]/x[-1])^(-1))/pp      
    }   
    est <- tmp*(1-m1) 
  }
  #
  if (method=="PRBstar"){
    k0h <- threshold.hill(x, method="minmse", shape2=x.rho, scale2=x.beta)$threshold
    n<-length(x[x>0])
    evi <- evindex(x, log.x=FALSE, method="chill", threshold=k0h)
    phi <- (1-x.rho/2)-sqrt((1-x.rho/2)^2-0.5)
    pp <- phi/evi
    print(c(pp,pp))
    m1 <- ((x.beta*(1-phi))/(1-x.rho-phi))*(n/(1:(n-1)))^(x.rho)   
    if (pp==0){
      x<-log(x[x>0]) 
      xlog<-x#rev(x)
      tmp<-cumsum(xlog)/(1:n)  
      tmp <- (tmp[-n]-xlog[-1])
    }
    else {
      x0<-(x[x>0])^pp 
      x<-x0#(rev(x0))
      tmp<-(cumsum(x)/(1:n))
      tmp <- (1-(tmp[-n]/x[-1])^(-1))/pp      
    }   
    est <- tmp*(1-m1) 
  }
  #
  if (method=="CHp"){
    n <- length(x[x>0])
    pp <- param[1]
    if (pp == 0){
      x <- log(x[x>0]) 
      xlog<-x#rev(x)
      tmp <- cumsum(xlog)/(1:n)  
      tmp <- (tmp[-n]-xlog[-1])
    }
    else {
      x0 <- (x[x>0])^pp 
      x <- x0#(rev(x0))
      tmp <- (cumsum(x)/(1:n))
      tmp <- (1-(tmp[-n]/x[-1])^(-1))/pp      
    }
    m1 <- ((x.beta*(1-pp*tmp))/(1-x.rho-pp*tmp))*(n/(1:(n-1)))^(x.rho)   
    est <- tmp*(1-m1) 
  }
  if (method=="CHstar"){
    k0h <- threshold.hill(x, method="minmse", shape2=x.rho, scale2=x.beta)$threshold
    n<-length(x[x>0])
    evi <- evindex(x, log.x=FALSE, method="chill", threshold=k0h)
    phi <- (1-x.rho/2)-sqrt((1-x.rho/2)^2-0.5)
    pp <- phi/evi
    n <- length(x[x>0])
    if (pp == 0){
      x <- log(x[x>0]) 
      xlog<-x #rev(x)
      tmp <- cumsum(xlog)/(1:n)  
      tmp <- (tmp[-n]-xlog[-1])
    }
    else {
      x0 <- (x[x>0])^pp 
      x <- x0 #(rev(x0))
      tmp <- (cumsum(x)/(1:n))
      tmp <- (1-(tmp[-n]/x[-1])^(-1))/pp      
    }
    m1 <- ((x.beta*(1-pp*tmp))/(1-x.rho-pp*tmp))*(n/(1:(n-1)))^(x.rho)   
    est <- tmp*(1-m1) 
  }
  #
  if (method=="ppwm"){
    n <- length(x)
    a0 <- (cumsum(as.numeric(x))/(1:n))[-1]
    #a1 <- ((cumsum((0:(n-1))*x))[-1])/(as.numeric(2:n)*as.numeric(1:(n-1)))
    a1 <- ((cumsum(as.numeric((0:(n-1))*x)))[-1])/(2:n)
    a1 <- a1/(1:(n-1))
    est <- 1-a1/(a0-a1)
  }
  if (method=="plpwm"){
    if (log.x==FALSE){x <- log(x[x>0])}
    n <- length(x)
    h1 <- (cumsum(x)/(1:n))[-1]
    h2 <- ((cumsum((0:(n-1))*x))[-1])/(1:(n-1))
    est <- 2*h1-4*(h2/(2:n))
  }
  if (method=="ML"){
    if (log.x==FALSE){x<-log(x[x>0])}  
    n <- length(x)
    xlog <- sort(x, decreasing = T)
    ui <- (xlog[-n]-xlog[-1])
    m0 <- cumsum((1:(n-1))*ui)/(1:(n-1))
    m1 <- cumsum((1:(n-1))^(1-x.rho)*ui)/((1:(n-1))^(1-x.rho))
    est <- m0-x.beta*(n/(1:(n-1)))^(x.rho)*m1
  }
  if (method=="MLb"){
    if (log.x==FALSE){x<-log(x[x>0])}  
    n <- length(x)
    xlog <- sort(x, decreasing = T)
    ui <- (xlog[-n]-xlog[-1])
    est <- cumsum((exp(-x.beta*(n/(1:(n-1)))^x.rho))*(1:(n-1))*ui)/(1:(n-1))
  }    
  if (method=="NML"){
    if (log.x==FALSE){x<-log(x[x>0])}  
    n <- length(x)
    xlog <- sort(x, decreasing = T)
    ui <- (xlog[-n]-xlog[-1])
    est <- cumsum((2-exp(x.beta*(n/(1:(n-1)))^x.rho))*(1:(n-1))*ui)/(1:(n-1))
  }
  if (method=="moment"){
    if (log.x==FALSE){x<-log(x[x>0])}
    n <- length(x)
    x2 <- x^2
    t1 <- cumsum(x)/(1:n)
    t2 <- cumsum(x2)/(1:n)
    m1 <- t1[-n]-x[-1]
    m2 <- t2[-n]-(2*t1[-n]*x[-1])+x2[-1]
    est <- m1+1-0.5*(1-m1^2/m2)^(-1)
  }
  ### WTC
  if (method=="bbtv"){
    n <- length(x)
    a0 <- (cumsum(as.numeric(x))/(1:n))[-n]
    est <- (a0/as.numeric(x)[-1]-1)
  }
  # G WTC (numerator=Hill)
  if (method=="gwtc"){
    if (log.x==FALSE) x <- log(x[x>0])
    n <- length(x)
    tmp <- cumsum(x)/(1:n)
    a0 <- tmp[-n]-x[-1]
    a1 <- (n+1)/((1:(n-1)))
    a1 <- cumsum(log(log(a1)))/(1:(n-1))
    a2 <- (n+1)/((2:n))
    a2 <- log(log(a2))
    est <- a0/(a1-a2)
  }
  
  #print(est)
  if (is.nan(est[1])) print (c(x.rho,x.beta))
  ifelse(is.na(threshold), return(est), return(est[threshold]))
}


#
# Adaptive estimation of the threshold and the extreme value index
#
threshold.hill <- function(x, method="gh", shape2=NA, scale2=NA, ccrit=1.25, plot=FALSE,  legend = TRUE, ...){ 
  stopifnot(is.numeric(x))
  # Direct estimation of the optimal threshold computed in Hall (1982)
  if (method=="minmse"){
    x <- sort(x[x>0])
    n <- length(x)
    if (is.na(shape2)){shape2 <- adapt.shape2(x)}
    if (is.na(scale2)){scale2 <- scale2(x, shape2=shape2)}
    thr <- (1-shape2)^2*n^(-2*shape2)/(-2*shape2*scale2^2)
    thr <- floor(thr^(1/(1-2*shape2)))
    xx<-log(x)    
    ui<-xx[-1]-xx[-length(xx)]
    ui<-(1:length(ui))*rev(ui)
    hill <- cumsum(ui)/(1:length(ui))
  }
  # method suggested in Guillou and Hall (2001)
  if (method=="gh"){
    xx<-log(sort(x[x>0]))
    nn<-length(xx)
    ui<-xx[-1]-xx[-length(xx)]
    ui<-(1:length(ui))*rev(ui)
    hill <- cumsum(ui)/(1:length(ui))
    s0 <- sqrt(3/((1:(nn-1))^1))
    s2 <- cumsum(ui)
    s1 <- rep(NA,nn-1)
    for (k in 1:(nn-1)){
      s1[k] <- sum((k-2*(1:k)+1)*ui[1:k])
    }  
    tk <- abs(s0*s1/s2)
    tk2 <- tk^2
    qk <- rep(NA,nn-1)
    for (k in 1:(nn-1)){
      qk[k] <- sum(tk2[(k-floor(k/2)):(k+floor(k/2))])/(2*floor(k/2)+1)
      #mean(tk2[(k-floor(k/2)):(k+floor(k/2))])
    }
    # Remove NAs
    qk <- sqrt(na.omit(qk))
    if (plot){
      plot(qk, xlab="k", ylab="", type="l",...)
      abline(h=ccrit, col=4, lty=2)
      if (legend) legend("topleft", legend=c(expression(Q[n](k)),expression(c[crit])),col=c(1,4), lty=1:2, bty="n", cex=1.4)
    }
    thr <- NA
    for (i in length(qk):1){
      if (min(qk[i:length(qk)])>ccrit){thr <- i}
    } 
    if (is.na(thr)) hill <- NA
  }
  return(list(threshold=thr,hill=hill[thr]))  
}


#
# Adaptive estimation of the threshold based on sample path stability
#
threshold.ps <- function(x, method="hill", shape2=NA, scale2=NA){
  stopifnot(is.numeric(x))
  # STEP 1
  evi <- evindex(x, method=method)
  # STEP 2
  # Obtain the minimum value of dp, such that the rounded values, to dp decimal places, of the estimates are distinct.
  dp <- 0
  cont <- TRUE
  if (round(max(evi), digits=dp) != round(min(evi), digits=dp)){cont <- FALSE}
  while(cont){
    dp <- dp+1
    if (round(max(evi), digits=dp) != round(min(evi), digits=dp)){cont <- FALSE}
  }
  #
  # STEP 3
  evi.l <- lrun(evi, dp)
  #
  # STEP 4
  evi2 <- round(evi[evi.l[1]:evi.l[2]], dp+2)
  evi.m <- sample.mode(evi2)
  evi0 <- evi[evi.l[1]+evi.m$max.mode-1]
  return(list(threshold=evi.l[1]+evi.m$max.mode-1,evi=evi0))  
}



#
# Adaptive estimation of the threshold and the extreme value index using the bootstrap
#
threshold.boot <- function(x, method="bootstrap", n1=floor(length(x)^(.955)), r=1, b1=250, rseed=NA, evi="hill", k.aux=2*sqrt(length(x[x>0])), RB=F){  
  rmse <- NULL
  output.k <- output.k2 <- NULL
  output.e <- output.e2 <- NULL
  log.x <- F
  if (evi=="hill" | evi=="chill"){
    x <- log(x[x>0])
    log.x <- T
  }
  x.n<-length(x)
  # estimation of the second order parameters
  x.rho <- adapt.shape2(x, log.x=log.x, threshold=floor(x.n^(.995)))
  x.beta <- scale2(x, log.x=log.x, threshold=floor(x.n^(.995)), shape2=x.rho)
  x.evi <- evindex(x, log.x=log.x, method=evi, shape2=x.rho, scale2= x.beta)
  if (max(n1) > x.n){stop("n1 must be smaller or equal to n")}
  bsn <- ifelse(max(n1) > x.n, max(n1), x.n)
  c <- ifelse(RB, 2, 1)
  if (method=="hall"){
    x.evi <- evindex(x, log.x=log.x, method="hill")  
    x.evi0 <- x.evi[k.aux]
    for (i in n1){
      if (is.na(rseed) == F){set.seed(rseed)}
      k0.r <- NULL
      for (j in 1:r){
        tmse <- tbias <- 0
        for (k in 1:b1){
          bss <- sample(x, bsn, replace=T)         
          bss1 <- sort(bss[1:i])
          t1 <- evindex(bss1, log.x=log.x, method="hill")
          tmse <- tmse+(t1-x.evi0)^2
          tbias <- tbias+(t1-x.evi0)         
        }
        k1 <- which.min(tmse[-1]/b1)+1
        k0 <- floor(k1*(x.n/i)^(-2*c*x.rho/(1-2*c*x.rho)))
        k0.r <- c(k0.r,k0)
      }
      output.k <- cbind(output.k, k0.r) 
      output.e <- cbind(output.e, x.evi[k0.r]) 
    }
    return(list(n1=n1, 
                k0=output.k, k0mean=apply(output.k, 2, mean), 
                est=output.e, estmean=apply(output.e, 2, mean), 
                rho=x.rho, c=c)) 
  }  
###  
  if (method=="bootstrap" | method=="doublebootstrap"){  
  for (i in n1){
    if (is.na(rseed) == F){set.seed(rseed)}
    k0.r <- k0.r2 <- NULL
    for (j in 1:r){
      tmse <- tmse2 <- tbias <- 0
      for (k in 1:b1){
        bss <- sample(x, bsn, replace=T)         
        bss1 <- sort(bss[1:i])
        t1 <- evindex(bss1, log.x=log.x, method=evi, shape2=x.rho, scale2= x.beta)
        t2 <- t1[1:(length(t1))/2]
        t2 <- c(t2[1], t2)
        tmse <- tmse+(t2-t1)^2
        tbias <- tbias+(t2-t1) 
        if (method=="doublebootstrap"){
          n2 <- floor(i^2/x.n)+1
          bss2 <- sort(bss[1:n2])
          t1 <- evindex(bss2, log.x=log.x, method=evi, shape2=x.rho, scale2= x.beta)
          t2 <- t1[1:(length(t1))/2]
          t2 <- c(t2[1], t2)
          tmse2 <- tmse2+(t2-t1)^2
        }         
      }
      k1 <- which.min(tmse[-1]/b1)
      #x.rho<-log(k1)/(2*log(k1/i))
      k0 <- floor(k1*(1-2^(c*x.rho))^(2/(1-2*c*x.rho))*(x.n/i)^(-2*c*x.rho/(1-2*c*x.rho)))+1
      #print(c(x.rho,x.beta,c))
      #k0 <- k1
      # double bootstrap
      if (method == "doublebootstrap"){
        k2 <- which.min(tmse2[-1]/b1)
        k0d <- floor((1-2^(c*x.rho))^(2/(1-2*c*x.rho))*(k1^2)/k2)+1
        if (k0d>x.n) k0d<-x.n-1
      }
      # new
      if (method=="new"){
        tmp<-1-(x.n/i)^(2*x.rho)/((2^x.rho-1)^2)
        #print(x.n)
        amse<-tmse/b1-((tbias/b1)^2)*tmp
        #plot(2:length(amse), amse[-1], type="l", ylim=c(0,median(amse)))
        k0 <- which.min(amse[-1])+1
      }
      k0.r <- c(k0.r,k0)
      k0.r2 <- c(k0.r2,k0d)
    }
    output.k <- cbind(output.k, k0.r) 
    output.e <- cbind(output.e, x.evi[k0.r]) 
    output.k2 <- cbind(output.k2, k0.r2) 
    output.e2 <- cbind(output.e2, x.evi[k0.r2]) 
    
  }
  #return(list(kk=n1, k11=k0bs1.1, k12=k0bs1.2, bs11=bs1.1, bs12=bs1.2))
  #colnames(output.k) <- n1
  #colnames(output.e) <- n1 
  return(list(n1=n1, 
              k0=output.k, k0mean=apply(output.k, 2, mean), 
              est=output.e, estmean=apply(output.e, 2, mean), 
              k0d=output.k2, k0meand=apply(output.k2, 2, mean), 
              estd=output.e2, estmeand=apply(output.e2, 2, mean), 
              rho=x.rho, c=c)) 
}  

}


#
# Adaptive estimation of the threshold and the extreme value index using the bootstrap
#
threshold.spboot <- function(x, est="hill", thr=2:(length(x)-1), r=1, b1=250, rseed=NA){
  stopifnot(is.numeric(x))
  output.k <- NULL
  output.e <- NULL  
  x<-sort(x)
  n<-length(x)
  log.x <- F
  x.n<-length(x)
  if (max(thr) > n-1){stop("threshold must be smaller than n")}
  # estimation of parameters
  x.rho <- adapt.shape2(x, log.x=log.x, threshold=floor(x.n^(.995)))
  x.beta <- scale2(x, log.x=log.x, threshold=floor(x.n^(.995)), shape2=x.rho)
  x.evi <- evindex(x, log.x=log.x, method=est, shape2=x.rho, scale2= x.beta)
  output.mse<-NULL
  method<-"caers"
  if (method=="caers"){
    #requires ismev
    sc <- NULL
    sh <- NULL
    for (i in 1:length(thr)){
      f <- gpd.fit(x, threshold=x[n-thr[i]], show=F)$mle
      sc <- c(sc, f[1])
      sh <- c(sh, f[2])
    }
    for (i in 1:length(thr)){
      if (!is.na(rseed)){set.seed(rseed)}
      mse.r <- NULL
      for (j in 1:r){
        t0 <- rep(NA,b1)
        for (k in 1:b1){
          un <- runif(n)
          unt <- un[un>1-thr[i]/x.n]
          unt <- 1-(1-unt)*(x.n/thr[i])
          # Generate tail sample from the GPD model
          bss <- x[n-thr[i]]+sc[i]*(unt^(-x.evi[thr[i]])-1)/x.evi[thr[i]]       
          # Generate the remaining values from the sample
          bss0 <- sample(x[1:(x.n-thr[i])], x.n-length(bss), replace=T)
          bss1 <- c(bss0, bss)
          t0[k] <- evindex(bss1, log.x=log.x, method=est, threshold=thr[i], shape2=x.rho, scale2= x.beta)
        }
        tbias <- mean(t0)-x.evi[thr[i]]
        tvar <- var(t0)
        #mse.r <- c(mse.r,tbias^2+tvar)
        mse.r <- c(mse.r,mean((t0-x.evi[thr[i]])^2))
        #k1 <- which.min(tmse[-1]/b1)+1
        #k0 <- floor(k1*(x.n/i)^(-2*c*x.rho/(1-2*c*x.rho)))
        
      }
      output.mse <- cbind(output.mse, mse.r)
    }
    k0.r <-  apply(output.mse, 1, which.min)#which.min(output.mse)
    output.k <- cbind(output.k, thr[k0.r]) 
    output.e <- cbind(output.e, x.evi[thr[k0.r]])
    colnames(output.mse) <- NULL
    colnames(output.k) <- NULL
    colnames(output.e) <- NULL
    plot(thr,apply(output.mse, 2, mean), type="l")
    return(list(thr=thr,
                mse.r=output.mse, 
                mse.mean=apply(output.mse, 2, mean),
                k0=output.k, k0mean=apply(output.k, 2, mean),
                est=output.e, est.mean=apply(output.e, 2, mean) 
                #rho=x.rho, c=c)
                )
                ) 
  }  
  
}