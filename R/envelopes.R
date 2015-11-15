#' Simulated Evelope
#' @description A normal plot with simulated envelope of the residual is produced.
#'
#'@usage envelope(model,k=19,alpha=0.05,res="deviance", precision = c("fixed","varying"),
#' dist = RBS(mu.link = "identity",sigma.link = "identity") )
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param alpha value giving the size of the envelope. Default is 0.05 which is equivalent to a 95\% band.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#' @param dist The function RBS() defines the RBS distribution.
#'
#'@return A simulated envelope of the class RBS.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'
#'
#' @examples
#'
#'## Fixed Precision
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'set.seed(2015)
#'envelope(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'set.seed(2015)
#'envelope(model)
#'
#'##
#'library(alr3)
#'data(landrent)
#'attach(landrent)
#'resp <- I(Y/X1)
#'y1 <-  split(resp, X4)$"1"
#'x21 <-  split(X2, X4)$"1"
#'
#'##Fixed Precision
#'fit0 <- gamlss(y1 ~ x21, family=RBS(mu.link="identity"),method=CG()  )
#'summary(fit0)
#'set.seed(2015)
#'envelope(fit0,alpha=0.01, precision="fixed",res="quantile",dist=RBS(mu.link="identity"))
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'summary(fit11)
#'set.seed(2015)
#'envelope(fit1,alpha=0.01, precision="varying",res="quantile",dist=RBS(mu.link="identity",sigma.link="sqrt"))
#' @export


envelope <- function(model,k=19,alpha=0.05,res="deviance", precision = c("fixed","varying"), dist = RBS(mu.link = "identity",sigma.link = "identity") )
{

  if(precision!= "fixed")
  {
  alfa1 <- ceiling(k*alpha)
  alfa2 <- ceiling(k*(1-alpha))
  n <- model$N
  td  <- residuals(model,residual=res)
  sigma <- model$sigma.fv
  mu <- model$mu.fv
  re <- matrix(0,n,k)
  X <- model$mu.x
  Z <- model$sigma.x
  p<- ncol(X)
  q<- ncol(Z)

  for(i in 1:k)
  {

    y1 <- mapply(rRBS,n=1,mu=mu,sigma=sigma)
    nresp <- y1

    if(p == 1) form <- nresp ~ 0 + X
    else form <- nresp ~ X[,-1]

    if(q == 1) form1 <- nresp ~ 0 + Z
    else  form1 <- nresp ~ Z[,-1]

    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
    model1 <- gamlss(formula=form,sigma.formula = form1 ,family=dist,method=CG(),control = conh0)
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
  qqnorm(td,xlab="Theorical Quantile",main="",ylab="Empirical Quantile",ylim=faixa,pch=20,cex=1,lwd=1)
  } else{
    alfa1 <- ceiling(k*alpha)
    alfa2 <- ceiling(k*(1-alpha))
    n<-model$N
    td  <- residuals(model, residual= res)
    sigma <- model$sigma.fv[1]
    mu <- model$mu.fv
    alpha<-sqrt(2/sigma)
    bet.a<- (sigma*mu)/(sigma+1)
    re <- matrix(0,n,k)
    X <- model$mu.x
    p <- ncol(X)
    mulink <- model$mu.link


    for(i in 1:k)
    {

      y1 <- mapply(rRBS,n=1,mu=mu,sigma=sigma)
      nresp <- y1
      if(p == 1) form <- nresp ~ 0 + X
      else  form <- nresp ~ X[,-1];

      conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula=form ,family=dist,method=CG(),control = conh0)
      rdf <- residuals(model1,residual=res)
      re[,i] <- sort(rdf)
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
    qqnorm(td,xlab="Theorical Quantile",main="",ylab="Empirical Quantile",ylim=faixa,pch=20,cex=1,lwd=1)
    }
  }





