bias <- function(type, xi, rho, zeta=1, r=1){
  if (type=="Hill") {b <- 1/(1-rho)}
  if (type=="MOM") {b <- (xi*(1-rho)+rho)/(xi*(1-rho)^2)}
  if (type=="MM") {b <- ((1+xi)*(xi+rho))/(xi*(1-rho)*(1+xi-rho))}
  if (type=="MOMRB") {
    #rho <- -xi/(1-xi)
    b <- (xi*(1-2*rho)+2*rho)*zeta/((1-2*rho)^2)-(rho*(1-2*rho-rho^2))/((1-2*rho)*(1-rho)^4)
    b <- b/(xi^2)
  }
  if (type=="MMRB") {
    rho <- -xi
    b <- (xi+2*rho)*zeta/((1-2*rho)*(1+xi-2*rho))
    #b <- b+((1-rho)*(1+2*xi)*(1+xi-rho)-(1+xi-2*rho)*(1+2*xi-rho))/((1+xi-rho)^2*(1+xi-2*rho)*(1-rho))
    #((1-rho)*(1+2*xi)*(1+xi-rho)-(1+xi-2*rho)*(1+2*xi+rho))
    b <- b-(rho*(1-3*rho+2*xi*(1+xi-rho)))/((1+xi-rho)^2*(1+xi-2*rho)*(1-rho))
    b <- (1+xi)*b/(xi^2)
  }  
  if (type=="GHRB") {
    #rho <- -xi/(1-xi)
    b <- (xi*(1-2*rho)+2*rho)*zeta/((1-2*rho)^2)-(rho)/((1-2*rho)*(1-rho)^2)
    b <- b/(xi^2)
  }  
  if (type=="CH") {b <- (zeta/(1-2*rho)-1/(1-rho)^2)/xi}
  if (type=="NCH") {b <- (zeta/(1-2*rho)-1.5/(1-rho)^2)/xi}
  if (type=="PPWM") {b <- (1-xi)*(2-xi)/((1-xi-rho)*(2-xi-rho))}
  if (type=="GPPWM"){b <- (xi+rho)*(1-xi)*(2-xi)/(xi*(1-xi-rho)*(2-xi-rho))}
  if (type=="PLPWM"){b <- 2/((1-rho)*(2-rho))}
  if (type=="MOP"){b <- (1-r)/(1-rho-r)}
  if (type=="RBMOP"){
    phi <- 1-rho/2-sqrt((1-rho/2)^2-1/2)
    b <- rho*(r-phi)/((1-r-rho)*(1-rho-phi))
  }
  if (type=="RB"){
    phi <- 1-rho/2-sqrt((1-rho/2)^2-1/2)
    b0 <- (1-phi)/(xi*(1-phi-rho)^2)
    b <- b0*((r*(1-phi-rho)^2+phi*rho-(1-phi)*(1-phi-2*rho))/(1-phi-2*rho))
  }
  if (type=="RB2"){
    phi <- 1-rho/2-sqrt((1-rho/2)^2-1/2)
    b0 <- (1-phi)/(xi*(1-phi-rho)^3)
    b <- b0*((r*(1-phi-rho)^3+phi*rho^2-(1-phi)*(1-phi-rho)*(1-phi-2*rho))/(1-phi-2*rho))
  }
  if (type=="KernelP") b <- r/(r-rho) #(1+r)/(1+r-rho)
  if (type=="RBKernelP") b <- (zeta*r/(r-2*rho)-r^2/(r-rho)^2)/xi
  if (type=="RBKernelL") b <- (zeta/(1-2*rho)^r-1/(1-rho)^(2*r))/xi
  if (type=="Lehmer" | type=="KernelL") {b <- 1/(1-rho)^r}
  if (type=="RBLehmer") {b <- ((zeta*rho+1-2*rho)/(1-2*rho)^r-(1-rho*(1-rho)^r)/(1-rho)^(2*r))/(xi*rho)}
  #if (type=="RBLehmer") {b <- ((zeta*rho+1-2*rho)/(1-2*rho)^r-((1-rho)^r+rho^2)/(1-rho)^(2*r))/(xi*rho)}
  if (type=="GM") {b <- (1-(1-rho)^r)/(r*rho*(1-rho)^r)}
    return(b)
}
#bias(type="CH", xi=1, rho=-1, r=1, w=0)
#bias(type="KernelP", xi=1, rho=-1, r=1, w=0)
sigma2 <- function(type, xi, rho, r=1){
  if(type=="Hill" | type=="CH"| type=="NCH") {s2 <- xi^2}
  if(type=="MOM" | type=="MOMRB" | type=="GHRB")  {s2 <- xi^2+1}
  if(type=="MM" | type=="MMRB")  {s2 <- (xi+1)^2}
  if(type=="PPWM") {s2 <- xi^2*(1-xi)*(2-xi)^2/((1-2*xi)*(3-2*xi))}
  if(type=="GPPWM"){s2 <- (1-xi+2*xi^2)*(1-xi)*(2-xi)^2/((1-2*xi)*(3-2*xi))}
  if(type=="PLPWM"){s2 <- xi^2*4/3}
  if(type=="MOP")  {s2 <- xi^2*(1-r)^2/(1-2*r)}
  if(type=="RBMOP"){s2 <- xi^2*(1-r)^2/(1-2*r)}
  if(type=="RB" | type=="RB2"){
    phi <- 1-rho/2-sqrt((1-rho/2)^2-1/2)
    s2 <- xi^2*(1-phi)^2/(1-2*phi)}
  if (type=="KernelP" | type=="RBKernelP") s2 <-  xi^2*r^2/(2*r-1) # xi^2*(1+r)^2/(1+2*r)
  if(type=="Lehmer" | type=="RBLehmer" | type=="KernelL" | type=="RBKernelL"){s2 <- xi^2*gamma(2*r-1)/(gamma(r)^2)}
  if(type=="GM"){s2 <- (xi^2/r^2)*(gamma(2*r+1)/gamma(r+1)^2-1)}
    return(s2)
}

#x<-10:80/20; plot(x, sigma2("Lehmer",xi=1, rho=0, r=x))
################ OLD
#areff<- function(est1, est2, xi, rho, r=1){
#  tmp <- ((sigma2(est2,xi,rho,r)/sigma2(est1,xi,rho,r))^(-rho)*abs(bias(est2,xi,rho,r)/bias(est1,xi,rho,r)))^(1/(1-2*rho))
#  tmp <- ifelse(tmp>10^6, 10^6, tmp)
#  return(tmp)
#}
#areff2<- function(est1, est2, xi, rho, r=1){
#  tmp <- ((sigma2(est2,xi,rho,r)/sigma2(est1,xi,rho,r))^(-2*rho)*abs(bias(est2,xi,rho,r)/bias(est1,xi,rho,r)))^(1/(1-4*rho))
#  tmp <- ifelse(tmp>10^6, 10^6, tmp)
#  return(tmp)
#}
################ OLD
#areff<- function(r, rho, est1, est2, xi, zeta, RB=TRUE){
#  cc <- 1
#  if (RB) cc<-2
#  tmp <- ((sigma2(est2,xi,rho,r)/sigma2(est1,xi,rho,r))^(-cc*rho)*abs(bias(est2,xi,rho,zeta,r)/bias(est1,xi,rho,zeta,r)))^(1/(1-2*cc*rho))
#  tmp <- ifelse(tmp>10^3, 10^3, tmp)
#  return(tmp)
#}

areff <- function(x, y, est1, est2, p1="r", p2="rho", xi=1, rho=-1, zeta=1, r=1, RB=TRUE){
  if (p1 == "r") r <- x
  if (p1 == "xi") xi <- x
  if (p1 == "rho") rho <- x
  if (p1 == "zeta") zeta <- x
  if (p2 == "r") r <- y
  if (p2 == "xi") xi <- y
  if (p2 == "rho") rho <- y
  if (p2 == "zeta") zeta <- y  
  if (p2 == "zeta") zeta <- y 
  if (p2 == "xi+rho") {
    xi <- y
    rho <- -xi/(1-xi)} 
  if (p2 == "MMxi+rho") {
    xi <- y
    rho <- -xi}   
  cc <- 1
  if (RB) cc<-2
  tmp <- ((sigma2(est2,xi,rho,r)/sigma2(est1,xi,rho,r))^(-cc*rho)*abs(bias(est2,xi,rho,zeta,r)/bias(est1,xi,rho,zeta,r)))^(1/(1-2*cc*rho))
  tmp <- ifelse(tmp>10^3, 10^3, tmp)
  return(tmp)
}




# AREFF (GM|H)
bname<-"fig_areff_gm-h"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x <- (51:200)/100
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="GM", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(alpha), xlim=c(0.5,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=1, cex.lab=1.0, 
      main=expression(AREFF[paste(K(alpha)," | ",H,sep="")]))
abline(h=seq(-1.5,-0.5, by=0.5), v=0, lty=3, col=gray(0.5))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x <- (101:500)/200
y <- -(400:1)/200
z <- outer(x, y, "areff", est1="GM", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x, y, z, add=TRUE, levels=c(0.5, 0.7, 0.8, 0.9, 0.96, 0.98, 0.99, 1.01, 1.02, 1.03, 1.04, 1.06, 1.08, 1.1), col="black", labcex=0.9, lwd=1)  #
contour(x, y, z, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

x <- (201:500)/200
y <- -(400:1)/200
z <- outer(x, y, "areff", est1="GM", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x, y, z, add=TRUE, levels=c(0.88, 0.92, .94), col="black", labcex=0.9, lwd=1, method="flattest")  #
#curve((1/2-x^2)/(1-x), from = 0.5, to = 1, add = TRUE, col = "green3")

#abline(v=seq(0.5,0.45, by=0.05), lty=3, col=gray(0.5))
axis(1, at=(6:24)/10,  labels=rep("",19))
box()

dev.off()











# Lehmer
#filename<-paste("lehmer",".pdf",sep="")
#pdf(file = filename, width=7, height=6)#, horizontal=F)
x<-(101:600)/200
y<--(400:1)/200
z <- outer(x, y, "areff", est1="Lehmer", est2="Hill", p1="r", p2="rho", RB=F)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(gamma), xlim=c(0.5,2), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[Lehmer/Hill]))
image(x, y, z, col=gray(0.99), axes=TRUE, zlim=c(0.999,1), add=TRUE)
image(x, y, z, col=gray(0.93), axes=TRUE, zlim=c(1.00,1.01), add=TRUE)
image(x, y, z, col=gray(0.86), axes=TRUE, zlim=c(1.01,1.02), add=TRUE)
image(x, y, z, col=gray(0.82), axes=TRUE, zlim=c(1.02,1.03), add=TRUE)
image(x, y, z, col=gray(0.75), axes=TRUE, zlim=c(1.03,1.04), add=TRUE)
image(x, y, z, col=gray(0.70), axes=TRUE, zlim=c(1.04,1.05), add=TRUE)
contour(x, y, z, add=TRUE, levels=c(0.95,0.995, 0.998, 1, 1.01, 1.02,1.03, 1.04), col="black", labcex=0.75, lwd=1, method="flattest")  #
abline(h=seq(-2.0,-0.5, by=0.5), v=0, lty=3, col=gray(0.3))
abline(v=seq(0.1,0.4, by=0.1), lty=3, col=gray(0.3))
#abline(-0.5,0,lty=3, col=gray(0.3))
box()
#dev.off()


# RB Lehmer v2
filename<-paste("rblehmer",".pdf",sep="")
pdf(file = filename, width=7, height=6)#, horizontal=F)
x <- c((101:195)/200, (1951:2049)/2000, (103:150)/100)
#x <- c((101:195)/200, (1951:2049)/2000, (103:150)/100)
#y<--c((101:195)/200, (1951:2049)/2000, (205:400)/200)
par(mar=c(5,4,4,2))
y<--(660:1)/300
z <- outer(x, y, "areff", est1="RBLehmer", est2="CH", p1="r", p2="rho", zeta=1, RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(alpha), xlim=c(0.5,1.5), ylim=c(-2.0,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[RB-Lehmer/CH]))
image(x, y, z, col=gray(0.99), axes=TRUE, zlim=c(0.999,1), add=TRUE)
image(x, y, z, col=gray(0.93), axes=TRUE, zlim=c(1.00,1.003), add=TRUE)
image(x, y, z, col=gray(0.91), axes=TRUE, zlim=c(1.003,1.007), add=TRUE)
image(x, y, z, col=gray(0.89), axes=TRUE, zlim=c(1.005,1.01), add=TRUE)
image(x, y, z, col=gray(0.86), axes=TRUE, zlim=c(1.01,1.02), add=TRUE)
image(x, y, z, col=gray(0.84), axes=TRUE, zlim=c(1.02,1.03), add=TRUE)
image(x, y, z, col=gray(0.82), axes=TRUE, zlim=c(1.03,1.04), add=TRUE)
image(x, y, z, col=gray(0.80), axes=TRUE, zlim=c(1.04,1.05), add=TRUE)
image(x, y, z, col=gray(0.78), axes=TRUE, zlim=c(1.05,1.06), add=TRUE)
image(x, y, z, col=gray(0.76), axes=TRUE, zlim=c(1.06,1.07), add=TRUE)
image(x, y, z, col=gray(0.74), axes=TRUE, zlim=c(1.07,1.10), add=TRUE)
image(x, y, z, col=gray(0.72), axes=TRUE, zlim=c(1.10,1.2), add=TRUE)
image(x, y, z, col=gray(0.70), axes=TRUE, zlim=c(1.2,1.35), add=TRUE)
image(x, y, z, col=gray(0.68), axes=TRUE, zlim=c(1.35,1.4), add=TRUE)
image(x, y, z, col=gray(0.66), axes=TRUE, zlim=c(1.4,1.5), add=TRUE)
image(x, y, z, col=gray(0.64), axes=TRUE, zlim=c(1.5,1.7), add=TRUE)
image(x, y, z, col=gray(0.62), axes=TRUE, zlim=c(1.7,1.8), add=TRUE)
image(x, y, z, col=gray(0.60), axes=TRUE, zlim=c(1.75,2), add=TRUE)
image(x, y, z, col=gray(0.58), axes=TRUE, zlim=c(2,3), add=TRUE)
contour(x, y, z, add=TRUE, levels=c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95,0.99, 1, 1.005, 1.01, 1.02,1.04, 1.1,1.25, 1.5, 2.5), col="black", labcex=0.75, lwd=1, method="flattest")  #
abline(h=seq(-2.5,0, by=0.1), v=0, lty=3, col=gray(0.3))
abline(v=seq(0.5,1.5, by=0.1), lty=3, col=gray(0.3))
#abline(-0.5,0,lty=3, col=gray(0.3))
box()
dev.off()


# RB Lehmer
filename<-paste("rblehmer",".pdf",sep="")
pdf(file = filename, width=7, height=6)#, horizontal=F)
x <- c((101:195)/200, (1951:2049)/2000, (103:150)/100)
#y<--c((101:195)/200, (1951:2049)/2000, (205:400)/200)
par(mar=c(5,4,4,2))
y<--(600:1)/300
z <- outer(x, y, "areff", est1="RBLehmer", est2="CH", p1="r", p2="rho", zeta=1, RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(alpha), xlim=c(0.5,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[RB-Lehmer/CH]))
image(x, y, z, col=gray(0.99), axes=TRUE, zlim=c(0.999,1), add=TRUE)
image(x, y, z, col=gray(0.93), axes=TRUE, zlim=c(1.00,1.01), add=TRUE)
image(x, y, z, col=gray(0.86), axes=TRUE, zlim=c(1.01,1.02), add=TRUE)
image(x, y, z, col=gray(0.82), axes=TRUE, zlim=c(1.02,1.04), add=TRUE)
image(x, y, z, col=gray(0.75), axes=TRUE, zlim=c(1.04,1.07), add=TRUE)
image(x, y, z, col=gray(0.725), axes=TRUE, zlim=c(1.07,1.10), add=TRUE)
image(x, y, z, col=gray(0.70), axes=TRUE, zlim=c(1.10,1.25), add=TRUE)
image(x, y, z, col=gray(0.675), axes=TRUE, zlim=c(1.25,1.5), add=TRUE)
image(x, y, z, col=gray(0.65), axes=TRUE, zlim=c(1.5,1.75), add=TRUE)
image(x, y, z, col=gray(0.62), axes=TRUE, zlim=c(1.75,2), add=TRUE)
image(x, y, z, col=gray(0.575), axes=TRUE, zlim=c(2,3), add=TRUE)
image(x, y, z, col=gray(0.55), axes=TRUE, zlim=c(3,100), add=TRUE)
contour(x, y, z, add=TRUE, levels=c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95,0.99, 1, 1.01, 1.02,1.04, 1.5, 2.5), col="black", labcex=0.75, lwd=1, method="flattest")  #
abline(h=seq(-2.0,0, by=0.1), v=0, lty=3, col=gray(0.3))
abline(v=seq(0.5,1.5, by=0.1), lty=3, col=gray(0.3))
#abline(-0.5,0,lty=3, col=gray(0.3))
box()
dev.off()


abline(h=-1.11, lwd=1.1)
x<-101:800/200; plot(x, areff(x=x, y=-1, est1="RBLehmer", est2="CH", p1="r", p2="rho", zeta=5/3, RB=T))

x<-seq(0.55,1.5,0.0001)
y<-areff(x,-1,est1="RBLehmer", est2="CH", p1="r", p2="rho", zeta=5/3, RB=T)
plot(x,y)


#############################################################
### Computation of \arg \max_p {\rm AREFF}_{\rm L}(p) ###
#############################################################
lmse<-function(p,rho=-1){-((gamma(p)^2/gamma(2*p-1))^(-rho)*(1-rho)^(p-1))^(1/(1-2*rho))}
optim(parm <- 1, fn=lmse, rho=-1, method="BFGS", hessian=F, control=list(maxit=2000))$par


lehmerp<-function(rho){
  lmse<-function(p,rho=-1){-((gamma(p)^2/gamma(2*p-1))^(-rho)*(1-rho)^(p-1))^(1/(1-2*rho))}
  pp<-rep(NA,length(rho))
  for (i in 1:length(rho))
  pp[i] <- optim(parm <- 1, fn=lmse, rho=rho[i], method="BFGS", hessian=F, control=list(maxit=20000))$par
  return(pp)
}

rr <- -rev(1:80/40)
out<-lehmerp(rr)
cbind(rr,out)
out2<-.1487*rr^3+.6537*rr^2+1.0621*rr+1.949
plot(out-out2)




#
# Mathematical Methods in Applied Sciences
#areff(x=0, y=-1, est1="KernelP", est2="Hill", p1="r", p2="rho", xi=1, rho=-1, zeta=1, r=1, RB=F)

  
#Kernel
# AREFF (KernelP|H)
filename <- "fig_areff_pk-h"
filename <- "figure2_above_areff_pk-h"
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript(file = paste(filename,".eps",sep=""), width=9, height=6)
#pdf(file = paste(filename,".pdf",sep=""), width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x <- (51:150)/100
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="KernelP", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(omega), xlim=c(0.5,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=1, cex.lab=1.0, 
      main=expression(AREFF[paste(P(omega)," | ",H,sep="")]))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2<-(101:300)/200
y2 <- -(400:1)/200
z2 <- outer(x2, y2, "areff", est1="KernelP", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x2, y2, z2, add=TRUE, levels=c(0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.01, 1.015, 1.02, 1.025), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

x2<-(200:300)/200
y2 <- -(400:1)/200
z2 <- outer(x2, y2, "areff", est1="KernelP", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x2, y2, z2, add=TRUE, levels=c(0.91, 0.92, 0.93, 0.94,0.96,0.97), col="black", labcex=0.9, lwd=1, method="flattest")  #

#curve((1/2-x^2)/(1-x), from = 0.5, to = 1, add = TRUE, col = "green3")
abline(h=seq(-1.5,-0.5, by=0.5), v=0, lty=3, col=gray(0.5))
#abline(v=seq(0.5,0.45, by=0.05), lty=3, col=gray(0.5))
axis(1, at=(5:15)/10,  labels=rep("",11))
box()
dev.off()


# Solve[(a^2/(2a-1))^(-r)*(a*(1-r)/(a-r))==1 && r=-1/2,a]



# AREFF (KernelL|H)
filename <- "fig_areff_lk-h"
filename <- "figure2_below_areff_lk-h"
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript(file = paste(filename,".eps",sep=""), width=9, height=6)
#pdf(file = paste(filename,".pdf",sep=""), width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x <- (51:200)/100
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="KernelL", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(omega), xlim=c(0.5,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=1, cex.lab=1.0, 
      main=expression(AREFF[paste(L(omega)," | ",H,sep="")]))
abline(h=seq(-1.5,-0.5, by=0.5), v=0, lty=3, col=gray(0.5))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x <- (101:500)/200
y <- -(400:1)/200
z <- outer(x, y, "areff", est1="KernelL", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x, y, z, add=TRUE, levels=c(0.5, 0.7, 0.8, 0.9, 0.96, 0.98, 0.99, 1.01, 1.02, 1.025, 1.03, 1.035, 1.04, 1.042), col="black", labcex=0.9, lwd=1)  #
contour(x, y, z, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

x <- (201:500)/200
y <- -(400:1)/200
z <- outer(x, y, "areff", est1="KernelL", est2="Hill", p1="r", p2="rho", zeta=1, RB=F)
contour(x, y, z, add=TRUE, levels=c(0.88, 0.92, .94), col="black", labcex=0.9, lwd=1, method="flattest")  #

#curve((1/2-x^2)/(1-x), from = 0.5, to = 1, add = TRUE, col = "green3")
#abline(v=seq(0.5,0.45, by=0.05), lty=3, col=gray(0.5))
axis(1, at=(6:24)/10,  labels=rep("",19))
box()
dev.off()




# AREFF (RBKernelP|CH)
filename <- "fig_areff_rbkp-ch"
filename <- "figure5_above_areff_rbkp-ch"
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript(file = paste(filename,".eps",sep=""), width=9, height=6)
#pdf(file = paste(filename,".pdf",sep=""), width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x <- (51:200)/100
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="RBKernelP", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)
par(mar=c(5,4,4,2))
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlim=c(0.5,2), ylim=c(-2,0), las=1,
      xlab=expression(omega),  ylab=expression(rho), cex.axis=1, cex.lab=1.0, 
      main=expression(AREFF[paste(tilde(P)(omega)," | ",CH,sep="")]))
abline(h=seq(-1.5,-0.5, by=0.5), v=0, lty=3, col=gray(0.5))
#abline(v=seq(0.25,2.25, by=0.25), lty=3, col=gray(0.5))
# lmt<-c(0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1, 1.05, 1.5, 2, 3, 4,5, 10, 100, 1e3)
# lmt<-c(0:32/32, 1.01, 1.025, 1.05, 1.075, 1:25/20*2+1, 4:9, 1:8*10, 1e6)
# nn<-length(lmt)
# n0<-length(lmt[lmt<1])
# mycol <- colorRampPalette(c('blue','cyan','green','yellow','red'), interpolate="spline")(nn)
# mygray <- 1-(1:(nn-1)/(2.8*nn+2))-c(rep(0,n0), rep(0.06,nn-1-n0))
# for (i in (1:(nn-1))){
#   image(x, y, z, col=gray(mygray[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
# }

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <- c((251:449)/500, (900:1100)/1000, 551:1000/500)
y2 <- -c((2000:1000)/1000, (199:1)/200)
z2 <- outer(x2, y2, "areff", est1="RBKernelP", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)
contour(x2, y2, z2, add=TRUE, levels=c(0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 1.01, 1.02, 1.05, 1.1, 1.2, 1.3,1.5,1.75, 2, 3, 4, 5), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")
#drawlabels
#
axis(1, at=(5:20)/10,  labels=rep("",16))
#abline(-0.5,0,lty=3, col=gray(0.3))
box()
dev.off()



# AREFF (RBKernelL|CH)
filename <- "fig_areff_rbkl-ch"
filename <- "figure5_below_areff_rbkl-ch"
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript(file = paste(filename,".eps",sep=""), width=9, height=6)
#pdf(file = paste(filename,".pdf",sep=""), width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x <- (51:150)/100
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="RBKernelL", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)
par(mar=c(5,4,4,2))
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlim=c(0.5,1.5), ylim=c(-2,0), las=1,
      xlab=expression(omega),  ylab=expression(rho), cex.axis=1, cex.lab=1.0, 
      main=expression(AREFF[paste(tilde(L)(omega)," | ",CH,sep="")]))
abline(h=seq(-1.5,-0.5, by=0.5), v=0, lty=3, col=gray(0.5))
#abline(v=1, col="blue")
# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <- c((251:449)/500, (900:1100)/1000, 551:750/500)
y2 <- -c((2000:1000)/1000, (199:1)/200)
z2 <- outer(x2, y2, "areff", est1="RBKernelL", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)
contour(x2, y2, z2, add=TRUE, levels=c(0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, .99, 1.01, 1.02, 1.05, 1.1, 1.2, 1.3,1.5,1.75, 2, 3, 4, 5), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")
#
x3 <- (401:600)/400
y3 <- -(500:1)/250
z3 <- outer(x3, y3, "areff", est1="RBKernelL", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)
contour(x3, y3, z3, add=TRUE, levels=c(0.92, 0.96, 0.995, 1.001, 1.002, 1.003), col="black", labcex=0.9, lwd=1, method="flattest")  #
#drawlabels
#
axis(1, at=(5:20)/10,  labels=rep("",16))
box()
dev.off()





areff(x=1, y=-2, est1="RBKernelP", est2="CH", p1="r", p2="rho", zeta=1, RB=TRUE)

# w optimo - vários estimadores RB
lmse <- function(w, rho=-1, zeta=1, type="RBKernelP"){
  v1 <- sigma2(type=type, xi=1, rho=rho, r=w)
  b1<-bias(type=type, xi=1, rho=rho, zeta=zeta, r=w)
  return((-4*rho)*log(v1)+2*log(abs(b1)))
}
#optim(parm <- 1, fn=lmse, rho=-1.75, type="RBkernelL", method="BFGS", hessian=F, control=list(maxit=5000))
w0RB <- function(rho=-1, zeta=1, type="RBKernelL"){
  w <- rep(NA, length(rho))
  for (i in 1:length(rho)){w[i] <- optim(parm <- 1, fn=lmse, rho=rho[i], type=type, method="L-BFGS-B", lower=0.5001, upper=100000, hessian=F, control=list(maxit=100000))$par}
  return(w)
}

w0RB(rho=-3, zeta=1, type="RBKernelL")
rr <- c(-3,-2.5,-2,-1.5,-1,-0.75,-0.5,-0.25,-0.1)
w0RBP <- w0RB(rr, zeta=1, type="RBKernelP")
w0RBL <- w0RB(rr, zeta=1, type="RBKernelL")

round(w0RBP,3)
round(w0RBL,3)

w0RB <- function(rho=-1, zeta=1, type="RBKernelP"){
  w <- rep(NA, length(rho))
  for (i in 1:length(rho)){w[i] <- optim(parm <- 1, fn=lmse, rho=rho[i], type=type, method="BFGS", hessian=F, control=list(maxit=100000))$par}
  return(w)
}
rr <- c(-3,-2.5,-2,-1.5,-1,-0.75,-0.5,-0.25,-0.1)
w0RBP <- w0RB(rr, zeta=1, type="RBKernelP")









cl <-c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", #"#F7FBFF",
       #"#FFF5F0", 
       "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
cl<-colorRampPalette(c("#08519C","white","#A50F15"), space = "Lab")(30)[-14:16]
cl<-colorRampPalette(c("blue","white","red"), space = "Lab")(40)
cl <- rev(heat.colors(40))
#Blue and Red
cl<-colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF","#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")), space = "Lab")(40)
filled.contour(x, y, z, 
               #color = col.l,
               #nlevels=40,
               plot.axes = { axis(1); axis(2); points(10, 10) },
               #levels=c(seq(from=0.3, to=1, length=9),seq(from=1, to=1.0422, length=9)[-1]),
               levels=c(seq(from=0.2, to=1, length=21),seq(from=1, to=1.0422, length=20)[-1]),
               col=cl
)
min(z)

seq(from=1, to=1.0422, length=9)[-1]





cl<-colorRampPalette(c("#08519C","white","#A50F15"), space = "Lab")(30)
barplot(rep(1,30), col=cl)

cl<-colorRampPalette(c('red','white','blue'), space = "Lab")(30)
barplot(rep(1,30), col=cl)

cl <- rev(heat.colors(30))
barplot(rep(1,30), col=cl)

# Colors from "RColorBrewer
cl <-c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", #"#F7FBFF",
       #"#FFF5F0", 
       "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
barplot(rep(1,16), col=cl)


#install.packages("RColorBrewer")
#barplot(rep(1,25), col=brewer.pal(n = 25, name = "RdBu"))
#barplot(rep(1,9), col=brewer.pal(n = 9, name = "Reds"))
#barplot(rep(1,9), col=brewer.pal(n = 9, name = "Blues"))


barplot(rep(1,30), col=rainbow(30))
barplot(rep(1,30), col=heat.colors(30))
barplot(rep(1,30), col=terrain.colors(30))
barplot(rep(1,30), col=topo.colors(30))
barplot(rep(1,30), col=cm.colors(30))

color=
barplot(rep(1,30), col=colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30))
#ou
clredblue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF","#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")
barplot(rep(1,30), col=colorRampPalette(rev(clredblue))(30))

#Oranges
barplot(rep(1,40), col=colorRampPalette(brewer.pal(11, "Oranges"))(40))












#Ivanilda - SPE2015 & ICCMSE2016
# AREFF (NCH|CH)
bname<-"areff_nch-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x<-c((-100:400)/200)
y <- -(300:1)/200
#m<-matrix(1,nrow = length(x), ncol = length(y))
#for (j in 1:length(y)){m[,j]<-areff(est1="NCH",est2="CH",xi=1,rho=y[j],r=x)}
z <- outer(x, y, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
par(mar=c(5,4,4,2))
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(-.5,2), ylim=c(-1.5,0), las=1, ylab=expression(rho), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(NCH," | ",CH,sep="")]))
#lmt <- quantile(m[m>=1],(1:25)/25)
lmt<-c(0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1, 1.05, 1.5, 2, 3, 4,5, 10, 100, 1e3)
lmt<-c(0:32/32, 1.01, 1.025, 1.05, 1.075, 1:25/20*2+1, 4:9, 1:8*10, 1e6)
nn<-length(lmt)
n0<-length(lmt[lmt<1])
mycol <- colorRampPalette(c('blue','cyan','green','yellow','red'), interpolate="spline")(nn)
mygray <- 1-(1:(nn-1)/(2.8*nn+2))-c(rep(0,n0), rep(0.06,nn-1-n0))
for (i in (1:(nn-1))){
  image(x, y, z, col=gray(mygray[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
  #image(x, y, m, col=rainbow(nn+2, start=i/(nn+5), end=(i+1)/(nn+5)), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
  #image(x, y, m, col=mycol[i], axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2<-c((-125:500)/250)
y2 <- -(450:2)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.5, 0.6, 0.7, 0.8, 0.9,1.1, 1.2, 1.5,2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=2, method="flattest")
#
x2<-c((-125:150)/250)
y2 <- -(450:3)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.75,.8+2*1:4/100, 1:9/100+.9), col="black", labcex=0.8, lwd=1, method="flattest")  #

x2<-c((150:165)/250)
y2 <- -(450:100)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=0.75, col="black", labcex=0.8, lwd=1, method="flattest", drawlabels = F) 
#
#falta 0.1, 3
x2<-c((500:1800)/1000)
y2 <- -(3000:2)/2000#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.1, 0.3, 3, 5), col="black", labcex=0.8, lwd=1, method="flattest")  #

x <- c((640:1000)/1000)
y <- -(1500:2)/500
z <-matrix(-1, nrow = length(x), ncol = length(y))
x0 <- round((1-2*y)/(1-y)^2,3)
for (i in 1:length(y)) z[which(x0[i]==x),i] <- 0
contour(x, y, z, add=TRUE, levels=0, col="red", labcex=0.8, lwd=2, method="flattest")  #

x <- c((800:1500)/1000)
y <- -(1500:2)/500
z <-matrix(-1, nrow = length(x), ncol = length(y))
x0 <- round(1.5*(1-2*y)/(1-y)^2,3)
for (i in 1:length(y)) z[which(x0[i]==x),i] <- 0
contour(x, y, z, add=TRUE, levels=0, col="green3", labcex=0.8, lwd=2, method="flattest")  #


# parte 3 - detalhes finais
abline(h=seq(-1.4,-0.2, by=0.2), v=0, lty=3, col=gray(0.3))
abline(v=seq(-0.25,1.75, by=0.25), lty=3, col=gray(0.3))
box()
legend("bottomright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("red","blue","green3"), ncol=3, cex=.85)
dev.off()


#v2 ficheiro menos pesado
# AREFF (NCH|CH)
bname<-"areff_nch-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# parte 1 - zonas a cinzento
x<-c((-75:300)/150)
y <- -(300:1)/200
#m<-matrix(1,nrow = length(x), ncol = length(y))
#for (j in 1:length(y)){m[,j]<-areff(est1="NCH",est2="CH",xi=1,rho=y[j],r=x)}
z <- outer(x, y, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
par(mar=c(5,4,4,2))
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(-.5,2), ylim=c(-1.5,0), las=1, ylab=expression(rho), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(NCH," | ",CH,sep="")]))
#lmt <- quantile(m[m>=1],(1:25)/25)
lmt<-c(0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1, 1.05, 1.5, 2, 3, 4,5, 10, 100, 1e3)
lmt<-c(0:32/32, 1.01, 1.025, 1.05, 1.075, 1:25/20*2+1, 4:9, 1:8*10, 1e6)
nn<-length(lmt)
n0<-length(lmt[lmt<1])
mycol <- colorRampPalette(c('blue','cyan','green','yellow','red'), interpolate="spline")(nn)
mygray <- 1-(1:(nn-1)/(2.8*nn+2))-c(rep(0,n0), rep(0.06,nn-1-n0))
for (i in (1:(nn-1))){
  image(x, y, z, col=gray(mygray[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
  #image(x, y, m, col=rainbow(nn+2, start=i/(nn+5), end=(i+1)/(nn+5)), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
  #image(x, y, m, col=mycol[i], axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2<-c((-125:500)/250)
y2 <- -(450:2)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5,2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=2, method="flattest")
#
x2<-c((-125:150)/250)
y2 <- -(450:3)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.75,.8+2*1:4/100, 1:9/100+.9), col="black", labcex=0.8, lwd=1, method="flattest")  #

x2<-c((150:165)/250)
y2 <- -(450:100)/300#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=0.75, col="black", labcex=0.8, lwd=1, method="flattest", drawlabels = F) 
#
#falta 0.1, 3
x2<-c((500:1800)/1000)
y2 <- -(3000:2)/2000#-(300:1)/200
z2 <- outer(x2, y2, "areff", est1="NCH", est2="CH", p1="zeta", p2="rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.1, 0.3, 3, 5), col="black", labcex=0.8, lwd=1, method="flattest")  #

x <- c((640:1000)/1000)
y <- -(1500:2)/500
z <-matrix(-1, nrow = length(x), ncol = length(y))
x0 <- round((1-2*y)/(1-y)^2,3)
for (i in 1:length(y)) z[which(x0[i]==x),i] <- 0
contour(x, y, z, add=TRUE, levels=0, col="red", labcex=0.8, lwd=2, method="flattest")  #

x <- c((800:1500)/1000)
y <- -(1500:2)/500
z <-matrix(-1, nrow = length(x), ncol = length(y))
x0 <- round(1.5*(1-2*y)/(1-y)^2,3)
for (i in 1:length(y)) z[which(x0[i]==x),i] <- 0
contour(x, y, z, add=TRUE, levels=0, col="green3", labcex=0.8, lwd=2, method="flattest")  #


# parte 3 - detalhes finais
abline(h=seq(-1.4,-0.2, by=0.2), v=0, lty=3, col=gray(0.3))
abline(v=seq(-0.25,1.75, by=0.25), lty=3, col=gray(0.3))
box()
legend("bottomright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("red","blue","green3"), ncol=3, cex=.85)
dev.off()



#ICRA5 - Caeiro & Gomes
#AREFF (PLPWM|H)
filename<-paste("caeirogomes_areff1",".pdf",sep="")
pdf(file = filename, width=7, height=6)#, horizontal=F)
x<-seq(-5,0,0.01)
y<-((3/4)^(-x)*abs(1-x/2))^(1/(1-2*x))
plot(x,y,type="l", las=1, lwd=3, col=1, xaxs="i", yaxs="i", ylim=c(.9,1.1),font.lab = 3, cex.axis=0.9, cex.lab=1.0, xlab=expression(rho), ylab=expression(AREFF), main=expression(AREFF[PLPWM/H]))
#legend("top", legend=c(expression(sigma[rho]^CG),expression(sigma[rho]^FAGH)),col=c(1,"gray"),lty=c(1,1),lwd=c(2,2),bty="n", horiz = TRUE, cex=1.2)
abline(h=0.9+0.05*(1:8),lty=3, col=gray(0.35))
abline(v=-5+1*(1:8),lty=3, col=gray(0.35))
dev.off()

# AREFF (PLPWM|PPWM)
filename<-paste("caeirogomes_areff2",".pdf",sep="")
pdf(file = filename, width=7, height=6)#, horizontal=F)
x<-(1:150)/300
y<--(400:1)/200
z <- outer(x, y, "areff", est1="PLPWM", est2="PPWM", p1="xi", p2="rho", RB=F)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(gamma), xlim=c(0,0.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[PLPWM/PPWM]))
image(x, y, z, col=gray(0.99), axes=TRUE, zlim=c(0.999,1), add=TRUE)
image(x, y, z, col=gray(0.93), axes=TRUE, zlim=c(1.00,1.01), add=TRUE)
image(x, y, z, col=gray(0.86), axes=TRUE, zlim=c(1.01,1.02), add=TRUE)
image(x, y, z, col=gray(0.82), axes=TRUE, zlim=c(1.02,1.05), add=TRUE)
image(x, y, z, col=gray(0.75), axes=TRUE, zlim=c(1.05,1.10), add=TRUE)
image(x, y, z, col=gray(0.70), axes=TRUE, zlim=c(1.10,1.25), add=TRUE)
image(x, y, z, col=gray(0.64), axes=TRUE, zlim=c(1.25,2), add=TRUE)
image(x, y, z, col=gray(0.59), axes=TRUE, zlim=c(2,20), add=TRUE)
contour(x, y, z, add=TRUE, levels=c(0.95,0.995, 0.998, 1, 1.01, 1.02,1.05, 1.10,1.25,2), col="black", labcex=0.75, lwd=1, method="flattest")  #
abline(h=seq(-2.0,-0.5, by=0.5), v=0, lty=3, col=gray(0.3))
abline(v=seq(0.1,0.4, by=0.1), lty=3, col=gray(0.3))
#abline(-0.5,0,lty=3, col=gray(0.3))
box()
dev.off()









# Artigo da Extremes
# AREFF (PRBMOP|H)
# parte 1 - zonas a cinzento
x<-c((1:200)/400)
y <- -(300:1)/200
m<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){m[,j]<-areff(est1="RBMOP",est2="Hill",xi=1,rho=y[j],r=x)}
par(mar=c(5,4,4,2))
image(x, y, m, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(p,xi,sep="")), xlim=c(0,0.5), ylim=c(-1.5,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[paste(p," | ",0,sep="")]))
lmt <- quantile(m[m>=1],(1:25)/25)
lmt[1:6]<-c(0.5, 1, 1.25, 1.5, 1.75, 2)
lmt[9]<-2.5
lmt[11]<-3
lmt[15]<-4
lmt[17]<-5
lmt[20]<-10
lmt[23]<-20
lmt[24]<-120
nn<-length(lmt)
for (i in (1:(nn-1))){
  image(x, y, m, col=gray(1-(i-0.4)/46), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}

# parte 2 - curvas de n??vel (necessita de uma matriz com mais linhas e colunas)
x2<-c((1:50)/500, (201:500)/2000, (251:400)/1000, (401:500)/1000)
y2 <- -c((6000:1501)/3000, (500:1)/1000)
m2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){m2[,j]<-areff(est1="RBMOP",est2="Hill",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, m2, add=TRUE, levels=c(0.01,0.5,1,1.5,2, 2.5,3:5,10,20, 50, 120), col="black", labcex=0.75, lwd=1, method="flattest")  #

# parte 3 - detalhes finais
abline(h=seq(-1.4,-0.2, by=0.2), v=0, lty=3, col=gray(0.3))
abline(v=seq(0.05,0.45, by=0.05), lty=3, col=gray(0.3))
box()








bname<-"areff_mop"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# AREFF (PRBMOP|H)
# parte 1 - zonas a cinzento
x<-(-124:124)/250
y <- -(300:1)/150
mop<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){mop[,j]<-areff(est1="MOP",est2="Hill",xi=1,rho=y[j],r=x)}
par(mar=c(5,4,4,2))
image(x, y, mop, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab="a", xlim=c(-0.5,0.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[paste(a," | ",0,sep="")]^{paste("*",sep="")}))
lmt<-c(0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.99, 0.999, 0.9999, 1, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01, 1.0125, 1.015, 1.0175, 1.02, 1.0225, 1.025, 1.1)
nn<-length(lmt)
for (i in (1:(nn-1))){
  image(x, y, mop, col=gray(1-(i/(2.8*nn+10))), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}
# 0 and 1; zero indicates "black" and one indicates "white".


# parte 2 - curvas de n??vel (necessita de uma matriz com mais linhas e colunas)
x2<- (-299:299)/600
y2 <- -(800:1)/400
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff(est1="MOP",est2="Hill",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.6, 0.8, 0.9, 0.95, 0.97, 0.99, 1, 1.005, 1.01, 1.015, 1.02, 1.023, 1.025), col="black", labcex=0.75, lwd=1, method="flattest")  #

# parte 3 - detalhes finais
abline(h=seq(-1.75,-0.25, by=0.25), v=0, lty=3, col=gray(0.3))
abline(v=seq(-0.45,0.45, by=0.05), lty=3, col=gray(0.3))
box()

dev.off()























bname<-"areff_rb-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# AREFF (RB|CH)
# parte 1 - zonas a cinzento
x<-(-100:400)/200
y <- -(300:1)/150
#x<-(-249:999)/500
#y <- -(800:1)/400
mop<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){mop[,j]<-areff2(est1="RB",est2="CH",xi=1,rho=y[j],r=x)}
par(mar=c(5,4,4,2))
image(x, y, mop, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(zeta), xlim=c(-0.5,2), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[paste(bar(RB)," | ",CH,sep="")]))
#lmt<-c(0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.99, 0.999, 0.9999, 1, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01, 1.0125, 1.015, 1.0175, 1.02, 1.0225, 1.025, 2)
#nn<-length(lmt)
lmt1<-c(0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.99,0.995, 0.9975, 0.999)
#lmt2<-c(1, 1+(1:19)/1000, 1+(8:39)/400, 1+(2:20)/20, 3:12,100, 200, 1000)
lmt2<-c(1, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01, 1.0125, 1.015, 1.0175, 1.02, 1.0225, 1.025, 1+5:10/100, 1.5, 2, 100, 300, 500, 10^6)
lmt <- c(lmt1, lmt2)
n1 <-length(lmt1)
n2 <-length(lmt2)
g1 <-(1-(1:n1/(4*n1+20)))
g2 <-(1-((n1+1):(n1+n2))/(4*n1+20))-.05
for (i in (1:(n1+n2-1))){
  image(x, y, mop, col=gray(c(g1,g2)[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}
# 0 and 1; zero indicates "black" and one indicates "white".


# parte 2 - curvas de n??vel (necessita de uma matriz com mais linhas e colunas)
x2<- (-399:1599)/800
y2 <- -(1599:1)/800
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0, 0.8, 0.9, 0.95, 0.975, 0.99, 1, 1.005, 1.007, 1.01, 1.015, 1.02, 1.025, 1.05, 1.1), col="black", labcex=0.75, lwd=1, method="flattest")  #

# parte 3 - detalhes finais
abline(h=seq(-1.75,-0.25, by=0.25), v=0, lty=3, col=gray(0.3))
abline(v=seq(-1,2, by=0.25), lty=3, col=gray(0.3))
box()
dev.off()








bname<-"areff_rb2-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# AREFF (RB2|CH)
# parte 1 - zonas a cinzento
x <- c(-74:299)/150
  #c((-199:100)/200, (251:499)/500, (100:199)/100)
y <- -(399:1)/200
mop<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){mop[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y[j],r=x)}
par(mar=c(5,4,4,2))
image(x, y, mop, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(zeta), xlim=c(-0.5,2), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[paste(bar(bar(RB))," | ",CH,sep="")]))
#lmt<-c(0, (25:49)/50, 0.99,0.995, 0.9975 , 0.999, 1, 1.001, 1.002, (51:80)/50, 2, 3, 4, 10,5000)
lmt1<-c(0, (25:49)/50, 0.99,0.995, 0.9975, 0.999)
lmt2<-c(1, 1+(1:19)/1000, 1+(8:39)/400, 1+(2:20)/20, 3:12,100, 200, 1000)
lmt <- c(lmt1, lmt2)
n1 <-length(lmt1)
n2 <-length(lmt2)
g1 <-(1-(1:n1/(4*n1+20)))
g2 <-(1-((n1+1):(n1+n2))/(4*n1+20))-.05

nn<-length(lmt)
for (i in (1:(nn-1))){
  #image(x, y, mop, col=gray(1-(i/(2.5*nn+10))), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
  image(x, y, mop, col=gray(c(g1,g2)[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}
# 0 and 1; zero indicates "black" and one indicates "white".


# parte 2 - curvas de n??vel (necessita de uma matriz com mais linhas e colunas)
#x2<- (-299:599)/300
x2<- c((-149:149)/300, (500:750)/1000, (226:599)/300)
y2 <- -(1999:1)/1000
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.75, 1, 1.2), col="black", labcex=0.75, lwd=1, method="flattest")  #
#0.9, 1.02

# parte 2.1
x2<- c((200:399)/200)
y2 <- -(399:1)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.95, 0.99,0.995), col="black", labcex=0.75, lwd=1, method="flattest")  #
x2<- (150:200)/200
y2 <- -(399:200)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.95), col="black", labcex=0.75, lwd=1, drawlabels = F)  #
x2<- (125:150)/200
y2 <- -(200:150)/100
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.95), col="black", labcex=0.75, lwd=1, drawlabels = F)  #
x2<- (180:220)/200
y2 <- -(110:90)/100
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0.95), col="black", labcex=0.75, lwd=1, drawlabels = F)  #


#x2<- (-149:399)/200
x2<- (300:399)/200
y2 <- -(47:1)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.01, 1.02), col="black", labcex=0.7, lwd=1, method="flattest")  #
x2<- (-99:199)/200
y2 <- -(20:1)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.01, 1.02, 1.03, 1.04, 1.05), col="black", labcex=0.7, lwd=1, method="flattest")  #
x2<- (-99:110)/200
y2 <- -(200:20)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.03, 1.04, 1.05), col="black", labcex=0.7, lwd=1, drawlabels = F)  #
x2<- (-99:100)/200
y2 <- -(400:201)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.03), col="black", labcex=0.7, lwd=1, method="flattest")  #
contour(x2, y2, mop2, add=TRUE, levels=c(1.04), col="black", labcex=0.7, lwd=1, drawlabels = F)  #
x2<- (50:99)/200
y2 <- -(400:200)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.05), col="black", labcex=0.7, lwd=1, drawlabels = F)  #
##################
x2<- (100:150)/200
y2 <- -(50:20)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.01, 1.02, 1.03, 1.04, 1.05), col="black", labcex=0.7, lwd=1, drawlabels = F)  #

x2<- (-100:50)/200
y2 <- -(399:300)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.01, 1.02), col="black", labcex=0.75, lwd=1, method="flattest")  #
x2<- (-100:0)/200
y2 <- -(299:200)/200
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="CH",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(1.02), col="black", labcex=0.75, lwd=1, drawlabels = F)  #


# parte 3 - detalhes finais
abline(h=seq(-1.75,-0.25, by=0.25), v=0, lty=3, col=gray(0.3))
abline(v=seq(-0.25,1.75, by=0.25), lty=3, col=gray(0.3))
box()
dev.off()






bname<-"areff_rb2-rb"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=9, height=6)
par(mar=c(5,4,4,2))

# AREFF (RB2|RB)
# parte 1 - zonas a cinzento
x<-(-124:499)/250
y <- -(300:1)/150
mop<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){mop[,j]<-areff2(est1="RB2",est2="RB",xi=1,rho=y[j],r=x)}
par(mar=c(5,4,4,2))
image(x, y, mop, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(zeta), xlim=c(-0.5,2), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[paste(bar(bar(RB))," | ",bar(RB),sep="")]))
lmt<-c(0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.99, 0.999, 0.9999, 1, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01, 1.0125, 1.015, 1.0175, 1.02, 1.0225, 1.025, 2)
nn<-length(lmt)
lmt1<-c(0, (25:49)/50, 0.99,0.995, 0.9975, 0.999)
lmt2<-c(1, 1+(1:19)/1000, 1+(8:39)/400, 1+(2:20)/20, 3:12,100, 150,200)
lmt <- c(lmt1, lmt2)
n1 <-length(lmt1)
n2 <-length(lmt2)
g1 <-(1-(1:n1/(4*n1+20)))
g2 <-(1-((n1+1):(n1+n2))/(4*n1+20))-.05
for (i in (1:(n1+n2-1))){
  #image(x, y, mop, col=gray(1-(i/(2.8*nn+10))), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE) 
  image(x, y, mop, col=gray(c(g1,g2)[i]), axes=TRUE, zlim=c(lmt[i],lmt[i+1]), add=TRUE)  
}
# 0 and 1; zero indicates "black" and one indicates "white".


# parte 2 - curvas de n??vel (necessita de uma matriz com mais linhas e colunas)
x2<- (-399:1199)/600
y2 <- -(800:1)/400
mop2<-matrix(1,nrow = length(x2), ncol = length(y2))
for (j in 1:length(y2)){mop2[,j]<-areff2(est1="RB2",est2="RB",xi=1,rho=y2[j],r=x2)}
contour(x2, y2, mop2, add=TRUE, levels=c(0, 0.6, 0.8, 0.9, 0.95, 0.97, 0.99, 1, 1.005, 1.01, 1.015, 1.02, 1.025, 1.05, 1.1, 2), col="black", labcex=0.75, lwd=1, method="flattest")  #

# parte 3 - detalhes finais
abline(h=seq(-1.75,-0.25, by=0.25), v=0, lty=3, col=gray(0.3))
abline(v=seq(-1,2, by=0.25), lty=3, col=gray(0.3))
box()
dev.off()


round(areff2(est1="RB2",est2="RB",xi=1,rho=-2,r=(0:30-10)/10),2)
round(areff2(est1="RB2",est2="CH",xi=1,rho=-2,r=(0:30-10)/10),2)
round(areff2(est1="RB2",est2="CH",xi=1,rho=-(0:25)/10,r=.5),2)
round(areff2(est1="RB2",est2="CH",xi=1,rho=-(0:25)/10,r=-.5),2)













bname<-"npa_fdens"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=7, height=6)
par(mar=c(5,4,4,2))

x<-seq(2,10, by=0.05)
y1<-1-8/(x^3)
y1<-24/(x^4)
x2<-seq(0,2, by=0.05)
y2<-rep(0, length(x2))
plot(x,y1,type="l",lty=1, col=4, ylab="", lwd=2, xlim=c(0,10), main="f(x)", las=1)
abline(h=seq(0,1.5, by=0.1), v=0, lty=3, col=gray(0.15))
abline(v=seq(0,10, by=0.5), lty=3, col=gray(0.15))
lines(x2,y2,lty=1,col=4, lwd=2)
dev.off()


x <- c(-149:599)/300
#c((-199:100)/200, (251:499)/500, (100:199)/100)
y <- -(399:1)/200
mop<-matrix(1,nrow = length(x), ncol = length(y))
for (j in 1:length(y)){mop[,j]<-areff2(est1="RB2",est2="RB",xi=1,rho=y[j],r=x)}
persp(x, y, mop,   theta = 25, phi = 30,   shade = 0.75, 
      ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "AREFF", zlim=c(0,2), col=hcl(240,100,80), border="black")

theta = 45, phi = 30)


png("~/persp.png", 500,500)
persp(u,v,z, ticktype="detailed", col=hcl(240,100,80), border="white", lwd=0.5)
dev.off()




bk<-c(0:30/30, 1:25/25*2+1, 1:8*10, 1e6)
pal.1=colorRampPalette(c("black", "red", "yellow"), space="rgb")(length(bk)-1)
pal.2=colorRampPalette(c("black", "blue", "cyan"), space="rgb")
mycol <- colorRampPalette(c('black','blue','cyan','green','yellow','red'), interpolate="spline")(length(bk)-1)
#par(mar=c(1,1,1,1))
x<-c((-100:400)/200)
y <- -(300:1)/200
#falta o z 
image(x, y, z, 
      col = rev(heat.colors(length(bk)-1)), 
      #col = mycol, 
      #col=pal.1,
      breaks=bk, 
      axes = T)
#grid()
par(mar=c(1,1,1,3))
image.scale(m, col=heat.colors(length(bk)-1), breaks=bk, horiz=FALSE, yaxt="n")



x <- seq(pi/4, 5 * pi, length.out = 100)
y <- seq(pi/4, 5 * pi, length.out = 100)
r <- as.vector(sqrt(outer(x^2, y^2, "+")))
grid <- expand.grid(x=x, y=y)
grid$z <- cos(r^2) * exp(-r/(pi^3))
levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = FALSE, region = TRUE)


x <- c((-100:400)/200)
y <- -(300:1)/200
m <- as.vector(outer(x, y, "areff", est1="NCH", est2="CH", xi=1))
grid <- expand.grid(x=x, y=y)
grid$z <- as.vector(outer(x, y, "areff", est1="NCH", est2="CH", xi=1))
levelplot(z~x*y, grid, cuts = 20, 
          at= c(0, .5, .6, .7, 0.75, .8, 0.85, .9, 0.95, 1, 1.025, 1.05, 1.075, 1.1, 1.2, 1.5, 2:5),
          xlab="", ylab="", 
          main="Weird Function", sub="with log scales",
          col.regions = rainbow(20), 
          colorkey = T, region = TRUE)


levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = FALSE, region = TRUE)


rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
heat.colors(n, alpha = 1)
terrain.colors(n, alpha = 1)
topo.colors(n, alpha = 1)
cm.colors(n, alpha = 1)


filled.contour(m, color = terrain.colors, asp = 1) # simple





#
#Linstat 2018
#
y<-seq(from=-10, to=0, by=0.1)
x<- -y/(1-y)
plot(x,y)

areff(x=.5,y=-1,est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)



# MOM/Hill
filename<-paste("areff_mom-hill",".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4,4,2))
x <- (1:225)/150
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
#zz<-sort(z, decreasing=T)
#zz<-round(log(zz[zz>1]),1)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(xi), xlim=c(0,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[list(M,GH)/H]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(-2.0,-0.25, by=0.25), v=0, lty=2, col=gray(0.4))
abline(v=seq(0,2, by=0.1), lty=2, col=gray(0.4))

lv<-seq(from=1, to=4, length.out=30)
lv[1]<-1.001
lv[30]<-3000000
#for (i in 1:29)  image(x, y, z, col=gray(0.80+i/50), axes=TRUE, zlim=c(lv[i],lv[i+1]), add=TRUE)

#image(x, y, z, col=gray(0.70), axes=TRUE, zlim=c(1.04,1.05), add=TRUE)
#contour(x, y, z, add=TRUE, levels=c(0.2,0.5,0.9,0.95,0.99, 1, 1.01, 1.02,1.03, 1.04, 1.1, 1.4,2,5), col="black", labcex=0.75, lwd=1, method="flattest")  #
x <- (1:30)/300
y <- -(800:200)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.1), col="black", labcex=0.8, lwd=1)  #
y <- -(200:1)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.1), col="black", labcex=0.75, lwd=1,drawlabels = F)  #
#
x <- (45:450)/300
y <- -(800:1)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.1,0.25, 0.5,.75,2), col="black", labcex=0.8, lwd=1)  #
contour(x, y, z, add=TRUE, levels=c(1), col="blue", labcex=0.8, lwd=1)  #
x <- (1:60)/400
y <- -(800:1)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.25, 0.5,.75,2), col="black", labcex=0.75, lwd=1,drawlabels = F)  #
contour(x, y, z, add=TRUE, levels=c(1), col="blue", labcex=0.75, lwd=1,drawlabels = F)  #
x <- (251:600)/400
y <- -(800:1)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.95,.97,.99,1.025,1.035, 1.05), col="black", labcex=0.8, lwd=1)  #
x <- (1:150)/300
y <- -(35:1)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.025,1.035,1.05), col="black", labcex=0.75, lwd=1, drawlabels = F)  #
x <- (150:251)/400
y <- -(150:5)/400
z <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.025,1.035,1.05), col="black", labcex=0.75, lwd=1, drawlabels = F)  #

#x <- (1:700)/1000
#y <- -(2000:1)/1000
#z <-matrix(-1, nrow = length(x), ncol = length(y))
#x0 <- round(-y/(1-y),3)
#for (i in 1:length(y)) z[which(x0[i]==x),i] <- 0
#contour(x, y, z, add=TRUE, levels=0, col="red3", labcex=0.8, lwd=2, method="flattest",drawlabels = F)  #
#z <-matrix(-1, nrow = length(x), ncol = length(y))
curve((-x)/(1-x), from = 0, to = .24, add = TRUE, col = "red3")
curve((-x)/(1-x), from = 0.28, to = 2/3, add = TRUE, col = "red3")
text(x=.26, y = -1/3-.01, labels=expression(infinity), col="red3", cex=1.1)
box()
#axis(1, at=(5:15)/10,  labels=rep("",11))
dev.off()


# MM/Hill
filename<-paste("areff_mm-hill",".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4,4,2))
x <- (1:300)/200
y <- -(200:1)/100
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
#zz<-sort(z, decreasing=T)
#zz<-round(log(zz[zz>1]),1)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.95,0.999), xlab=expression(xi), xlim=c(0,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main=expression(AREFF[list(MM)/H]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(-2.0,-0.25, by=0.25), v=0, lty=2, col=gray(0.4))
abline(v=seq(0,2, by=0.1), lty=2, col=gray(0.4))
#
x <- (1:600)/400
y <- -(600:230)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.1,0.2,.3,.4,.5,.6,.75), col="black", labcex=0.8, lwd=1)  #
x <- (1:200)/400
y <- -(230:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(0.1,0.2,.3,.4,.5,.6,.75), col="black", labcex=0.8, lwd=1,drawlabels = F)  #

x <- (60:560)/400
y <- -(600:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(.9,1.2,1.5,2), col="black", labcex=0.8, lwd=1)  #
contour(x, y, z, add=TRUE, levels=c(1), col="blue", labcex=0.8, lwd=1)  #

x <- (1:60)/400
y <- -(600:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(.9,1.2,1.5,2), col="black", labcex=0.8, lwd=1,drawlabels = F)  #
contour(x, y, z, add=TRUE, levels=c(1), col="blue", labcex=0.8, lwd=1,drawlabels = F)  #
x <- (560:600)/400
y <- -(600:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(.9,1.2,1.5,2), col="black", labcex=0.75, lwd=1,drawlabels = F)  #
contour(x, y, z, add=TRUE, levels=c(1), col="blue", labcex=0.8, lwd=1,drawlabels = F)  #
#
x <- (161:600)/400
y <- -(150:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.01,1.02,1.05,1.1), col="black", labcex=0.8, lwd=1, drawlabels = T)  #
x <- (80:160)/400
y <- -(75:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.01,1.02,1.05,1.1), col="black", labcex=0.8, lwd=1, drawlabels = F)  #
x <- (1:80)/400
y <- -(25:1)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.01,1.02,1.05,1.1), col="black", labcex=0.8, lwd=1, drawlabels = F)  #
x <- (401:600)/400
y <- -(260:150)/300
z <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
contour(x, y, z, add=TRUE, levels=c(1.01,1.02,1.05,1.1), col="black", labcex=0.8, lwd=1, drawlabels = F)  #

curve(-x, from = 0, to = .3, add = TRUE, col = "red3")
curve(-x, from = .35, to = 1.5, add = TRUE, col = "red3")
text(x=.33, y = -.32, labels=expression(infinity), col="red3", cex=1.1)

box()
#axis(1, at=(5:15)/10,  labels=rep("",11))
dev.off()






# AREFF (MOMRB|CH)
bname<-"areff_mom-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(-75:225)/150
y <- (1:149)/150
z <- outer(x, y, "areff", est1="MOMRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(-.5,1.5), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(MOM," | ",CH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,1.5, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(-400:1200)/800
y2 <- (1:795)/800
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(.1,.3,.45, 0.6, 0.7, 0.8, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

curve(sqrt(1-x), from = 0, to = 1, add = TRUE, col = "green3")
#curve((3*sqrt(3)+sqrt(27*x^2-28*x-2)-27*x+14)^(1/3)*(3*2^(2/3)), from = 0, to = 1, add = TRUE, col = "red3")
y <- seq(0,1, by=0.005)
x <- (1+y)*(1-2*y^2)
lines(x,y, col="red3")
box()
legend("topright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()


#final
# AREFF (MOMRB|CH)
bname<-"areff_mom-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(0:300)/150
y <- (1:149)/150
z <- outer(x, y, "areff", est1="MOMRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(0,2), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(MOM," | ",CH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,2, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(0:800)/400
y2 <- (1:595)/600
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(.1,.3,.45, 0.6, 0.7, 0.8, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

curve(sqrt(1-x), from = 0, to = 1, add = TRUE, col = "green3")
#curve((3*sqrt(3)+sqrt(27*x^2-28*x-2)-27*x+14)^(1/3)*(3*2^(2/3)), from = 0, to = 1, add = TRUE, col = "red3")
y <- seq(0,1, by=0.005)
x <- (1+y)*(1-2*y^2)
lines(x,y, col="red3")
box()
legend("topright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()
#areff(x=1, y=.366, est1="MOMRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)




# AREFF (GHRB|CH)
bname<-"areff_gh-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(-75:225)/150
y <- (1:149)/150
z <- outer(x, y, "areff", est1="GHRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(-.5,1.5), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(GH," | ",CH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,1.5, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(-300:900)/600
y2 <- (1:795)/800
z2 <- outer(x2, y2, "areff", est1="GHRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(.1,.3,.45, 0.6, 0.7, 0.8, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

curve(sqrt(1-x), from = 0, to = 1, add = TRUE, col = "green3")
curve(x-1, from = 0, to = 2, add = TRUE, col = "red3")
#y <- seq(0,1, by=0.005)
#x <- (1+y)*(1-2*y^2)
#lines(x,y, col="red3")
box()
legend("topright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()

#final
# AREFF (GHRB|CH)
bname<-"areff_gh-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(0:300)/150
y <- (1:149)/150
z <- outer(x, y, "areff", est1="GHRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(0,2), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(GH," | ",CH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,2, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(0:1200)/600
y2 <- (1:596)/600
z2 <- outer(x2, y2, "areff", est1="GHRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(.1,.3,.45, 0.6, 0.7, 0.8, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

curve(sqrt(1-x), from = 0, to = 1, add = TRUE, col = "green3")
curve(x-1, from = 0, to = 2, add = TRUE, col = "red3")
#y <- seq(0,1, by=0.005)
#x <- (1+y)*(1-2*y^2)
#lines(x,y, col="red3")
box()
legend("top", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()
#areff(x=5/3, y=.5, est1="GHRB", est2="CH", p1="zeta", p2="xi+rho", RB=T)



# AREFF (MOMRB|GHRB)
bname<-"areff_mom-gh"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(-50:150)/100
y <- (1:199)/200
z <- outer(x, y, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(-.5,1.5), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(MOM," | ",GH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,1.5, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(-400:1200)/800
y2 <- (1:715)/800
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.25, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

x2 <-(-400:1200)/800
y2 <- (715:798)/800
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
#contour(x2, y2, z2, add=TRUE, levels=c(.25, 0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.25, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest", drawlabels = F)

curve(x-1, from = 0, to = 2, add = TRUE, col = "green3")
y <- seq(0,1, by=0.005)
x <- (1+y)*(1-2*y^2)
lines(x,y, col="red3")
box()
legend("topright", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()

#final
# AREFF (MOMRB|GHRB)
bname<-"areff_mom-gh"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(0:200)/100
y <- (1:199)/200
z <- outer(x, y, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(0,2), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,"/(",1-rho,")",sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(MOM," | ",GH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,2, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(0:1000)/500
y2 <- (1:715)/800
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.25, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")

x2 <-(0:1000)/500
y2 <- (715:798)/800
z2 <- outer(x2, y2, "areff", est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)
#contour(x2, y2, z2, add=TRUE, levels=c(.25, 0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.25, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest", drawlabels = F)
contour(x2, y2, z2, add=TRUE, levels=c(.25, 0.5, 0.7, 0.8, 0.9, 0.95, 1.05, 1.25, 2), col="black", labcex=0.9, lwd=1, drawlabels = F)  #

curve(x-1, from = 0, to = 2, add = TRUE, col = "green3")
y <- seq(0,1, by=0.005)
x <- (1+y)*(1-2*y^2)
lines(x,y, col="red3")
box()
legend("top", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()
#areff(x=1, y=.62, est1="MOMRB", est2="GHRB", p1="zeta", p2="xi+rho", RB=T)




# AREFF (MMRB|CH)
bname<-"areff_mm-ch"
filename<-paste(bname,".pdf",sep="")
pdf(file = filename, width=8, height=6)
par(mar=c(4.5,4.2,4,2))

# parte 1 - zonas a cinzento
x <-(0:300)/150
y <- (1:149)/150
z <- outer(x, y, "areff", est1="MMRB", est2="CH", p1="zeta", p2="MMxi+rho", RB=T)
image(x, y, z, col=gray(1), axes=TRUE, zlim=c(0.01,1), xlab=expression(paste(zeta," = ",beta,"'/",beta,sep="")), xlim=c(0,2), ylim=c(0,1), las=1, ylab=expression(paste(xi==-rho,sep="")), cex.axis=1, cex.lab=1.0, main=expression(AREFF[paste(MM," | ",CH,sep="")]))
image(x, y, z, col=gray(0.85), axes=TRUE, zlim=c(1.001,10000), add=TRUE)
abline(h=seq(0,1, by=0.1), v=0, lty=2, col=gray(0.4))
abline(v=seq(-0.5,2, by=0.25), lty=2, col=gray(0.4))

# parte 2 - curvas de nivel (necessita de uma matriz com mais linhas e colunas)
x2 <-(0:800)/400
y2 <- (1:599)/600
z2 <- outer(x2, y2, "areff", est1="MMRB", est2="CH", p1="zeta", p2="MMxi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=c(.2,.3,.4, .5, 0.6, 0.8, 2), col="black", labcex=0.9, lwd=1, method="flattest")  #
x2 <-(700:1200)/600
y2 <- (1:399)/400
z2 <- outer(x2, y2, "areff", est1="MMRB", est2="CH", p1="zeta", p2="MMxi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest")
x2 <-(600:700)/600
y2 <- (1:399)/400
z2 <- outer(x2, y2, "areff", est1="MMRB", est2="CH", p1="zeta", p2="MMxi+rho", RB=T)
contour(x2, y2, z2, add=TRUE, levels=1, col="blue", labcex=0.8, lwd=1, method="flattest",drawlabels = F)
#curve(sqrt(1-x), from = 0, to = 1, add = TRUE, col = "green3")
y <- seq(0,1, by=0.00625)
x <- (1+2*y)/(1+y)^2
lines(x,y, col="green3")
#curve((3*sqrt(3)+sqrt(27*x^2-28*x-2)-27*x+14)^(1/3)*(3*2^(2/3)), from = 0, to = 1, add = TRUE, col = "red3")
#y <- seq(0,2, by=0.005)
#x <- (1+4*y)/(1+2*y)
#lines(x,y, col="red3")
curve((1-x)/(2*x-4), from = 1, to = 2, add = TRUE, col = "red3")
box()
legend("topleft", legend=c(0,1,expression(infinity)), lty = c(1,1,1), 
       lwd=c(2,2,2), col=c("green3","blue","red3"), ncol=3, cex=.85)
dev.off()
#areff(x=5/3, y=1, est1="MMRB", est2="CH", p1="zeta", p2="MMxi+rho", RB=T)

(1-5/3)/(2*5/3-4)


# M/GH/MM/Hill
#deepskyblue2
filename<-paste("areff_m_gh_mm_h",".pdf",sep="")
pdf(file = filename, width=8, height=6)
#c(bottom, left, top, right)
par(mar=c(4.5,4,1,2))
x <- (1:375)/250
y <- -(400:1)/200
z1 <- outer(x, y, "areff", est1="MOM", est2="Hill", p1="xi", p2="rho", RB=F)
z2 <- outer(x, y, "areff", est1="MM", est2="Hill", p1="xi", p2="rho", RB=F)
#z1[z1<1]<-0
z1<-ifelse(z1<1,0,z1)
z2<-ifelse(z2<1,0,z2)
z1<-ifelse(z1>z2,z1,0)
z2<-ifelse(z2>z1,z2,0)
image(x, y, z1, col="lightsteelblue1", axes=TRUE, zlim=c(.99,10^3), xlab=expression(xi), xlim=c(0,1.5), ylim=c(-2,0), las=1, ylab=expression(rho), cex.axis=0.9, cex.lab=1.0, main="")
image(x, y, z2, col="steelblue3", axes=TRUE, zlim=c(.99,10^3), add=TRUE)
abline(h=seq(-1.75,-0.25, by=0.25), v=0, lty=2, col=gray(0.4))
abline(v=seq(0.1,1.4, by=0.1), lty=2, col=gray(0.4))
text(0.25,-0.9,"H", cex=1.4)
text(1.05,-1.65,"H", cex=1.4)
text(1.05,-0.67,"MM", cex=1.4)
text(1.31,-0.18,"M, GH", cex=1.4)
text(0.61,-1.18,"M, GH", cex=1.4)
box()
dev.off()



### Ivanilda

ARE1 <- function(a,rho){
  xi <- 1
  S1 <- (((xi)^2)*((1-a)^2))/(1-2*a)
  S2 <- (xi)^2
  B1 <- (1-a)/(1-a-rho)
  B2 <- 1/(1-rho)
  Qs <- ((S2)/(S1))^(-rho)
  Qb <- abs((B2)/(B1))
  Q <- (Qs * Qb)^(1/(1-2*rho))
  return(Q)
}

x <- seq(from=0, to=0.48, by=0.02)
y <- seq(from=-2, to=0, by=0.1)
z <- outer(x, y, "ARE1")

#names(x) <- x
#names(y) <- y
#round(z,2)

contour(x, y, z, add=F, levels=c(0.95, 1, 1.01), col="black", labcex=0.75, lwd=1, method="flattest", main=expression(AREFF[.../H])) 
abline(h=seq(-2.0,0, by=0.1), v=0, lty=3, col=gray(0.5))
abline(v=seq(0,0.45, by=0.05), lty=3, col=gray(0.5))


ARE1(a=0.49,rho=c(-0.1,-2))
