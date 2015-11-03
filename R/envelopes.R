#' Simulated Evelope
#' @description A normal plot with simulated envelope of the residual is produced.
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param alfa value giving the size of the envelope. Default is 0.05 which is equivalent to a 95\% band.
#' @param res type of residuals to be extracted. Default is deviance.
#'
#' @export


envelope <- function(model,k=19,alfa=0.05,res="deviation")
{
  alfa1 <- ceiling(k*alfa)
  alfa2 <- ceiling(k*(1-alfa))
  n <- model$N
  td  <- residuals(model,residual=res)
  sigma <- model$sigma.fv
  mu <- model$mu.fv
  re <- matrix(0,n,k)
  X <- model$mu.x[,-1]
  Z <- model$sigma.x[,-1]

  for(i in 1:k)
  { #

    y1 <- mapply(rBS,n=1,mu=mu,sigma=sigma)
    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
    model1 <- gamlss(formula=y1~X,sigma.formula= ~ Z ,family=BS(mu.link="identity",sigma.link="sqrt"),method=CG(),control = conh0)

    rdf <- residuals(model1,residual=res)

    re[,i]=sort(rdf)
  }
  e1 = numeric(n)
  e2 = numeric(n)
  for(l in 1:n)
  {
    eo = sort(re[l,])
    e1[l] = eo[alfa1]
    e2[l] = eo[alfa2]
  }
  a<-  qqnorm(e1,plot.it=FALSE)$x
  a1<-  qqnorm(e1,plot.it=FALSE)$y
  b<-  qqnorm(e2,plot.it=FALSE)$x
  b1<-  qqnorm(e2,plot.it=FALSE)$y
  r<-  qqnorm(td,plot.it=FALSE)$x
  r1<-  qqnorm(td,plot.it=FALSE)$y
  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))

  xb = apply(re,1,mean)
  faixa = range(td,e1,e2)
  par(mar=c(4., 4., 0.1, 0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=TRUE)
  polygon(xx,yy,col="gray80",border=NA)
  par(new=TRUE)
  qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col="gray")
  par(new=TRUE)
  qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col="gray")
  par(new=TRUE)
  qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
  par(new=TRUE)
  qqnorm(td,xlab="perc",main="",ylab="Residuals",ylim=faixa,pch=8,cex=1,lwd=1)
}


envelope.bs=function(model,k=19,alfa=0.05,res="deviation")
{
  alfa1 = ceiling(k*alfa)
  alfa2 = ceiling(k*(1-alfa))
  n=model$N
  td  = residuals(model, residual= res)
  sigma = model$sigma.fv[1]
  mu = model$mu.fv
  alpha=sqrt(2/sigma)
  bet.a= (sigma*mu)/(sigma+1)
  re = matrix(0,n,k)
  X = model$mu.x[,-1]


  for(i in 1:k)
  {

    y1 <- mapply(rBS,n=1,mu=mu,sigma=sigma)

    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
    model1 <- gamlss(formula=y1~X ,family=BS(mu.link="identity"),method=CG(),control = conh0)
    rdf=residuals(model1,residual=res)
    re[,i]=sort(rdf)
  }
  e1 = numeric(n)
  e2 = numeric(n)
  for(l in 1:n)
  {
    eo = sort(re[l,])
    e1[l] = eo[alfa1]
    e2[l] = eo[alfa2]
  }

  a<-  qqnorm(e1,plot.it=FALSE)$x
  a1<-  qqnorm(e1,plot.it=FALSE)$y
  b<-  qqnorm(e2,plot.it=FALSE)$x
  b1<-  qqnorm(e2,plot.it=FALSE)$y
  r<-  qqnorm(td,plot.it=FALSE)$x
  r1<-  qqnorm(td,plot.it=FALSE)$y
  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))

  xb = apply(re,1,mean)
  faixa = range(td,e1,e2)
  par(mar=c(4., 4., 0.1, 0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=TRUE)
  polygon(xx,yy,col="gray80",border=NA)
  par(new=T)
  qqnorm(e1,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col="gray")
  par(new=TRUE)
  qqnorm(e2,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=1,lwd=1,col="gray")
  par(new=TRUE)
  qqnorm(xb,axes=F,xlab="",ylab="",main="",type="l",ylim=faixa,lty=2,lwd=1,col="black")
  par(new=TRUE)
  qqnorm(td,xlab="perc",main="",ylab="Residuals",ylim=faixa,pch=8,cex=1,lwd=1)
}

