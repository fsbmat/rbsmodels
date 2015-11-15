#' Diagnostic Analysis - Local Influnce
#' @description Diagnostics for the RBS model
#'
#'@usage diag.bs(model)
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#'
#'@return Local influence measures.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'
#'@examples
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'fit = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(fit)
#'set.seed(2015)
#'envelope(fit)
#'Cib <- diag.bs(fit)$Ci.beta
#'plot(Cib,ylim=c(0,1),pch=19)
#'abline(h=2*mean(Cib),type=2)
#'Cia <- diag.bs(fit)$Ci.alpha
#'plot(Cia,ylim=c(0,1),pch=19)
#'abline(h=2*mean(Cia),type=2)
#'@export

diag.bs=function(model,mu.link = "identity",sigma.link = "identity")
{

  x <- model$mu.x
  z  <- model$sigma.x
  y <- model$y
  p<-ncol(x)
  q<-ncol(z)

  linkstr <- mu.link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- sigma.link
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta


  B=function(Delta,I,M)
  {
    #I ? o inverso da matrix hessiana
    #M determina sobre qual sobre quais cojuntos de par?metros deseja-se estudar o influ?ncia
    #Delta ? a curvatura que determina qual o tipo de perturba??o deseja-se estudar
    B=(t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  loglik <- function(vP)
  {
    betab = vP[1:p]
    alpha = vP[-(1:p)]
    eta   = as.vector(x%*%betab)
    tau   = as.vector(z%*%alpha)
    mu    = linkinv(eta)
    sigma = sigma_linkinv(tau)


    f <- 0.5*sigma -0.5*log((sigma+1)) - 0.5*log(mu) - 1.5*log(y) + log((sigma*y) + y + (sigma*mu)) - y*(sigma+1)/(4*mu) - (sigma*sigma*mu)/(4*y*(sigma+1)) - 0.5*log(16*pi)
    return(sum(f))
  }

  muest <- model$mu.coefficients
  sigmaest <- model$sigma.coefficients
  x0<- c(muest,sigmaest)
  h0 <- hessian(loglik,x0)

  Ldelta= h0[(p+1):(p+q),(p+1):(p+q)]
  Lbeta=h0[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, q))
  b12=cbind(matrix(0, q, p), solve(Ldelta))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(solve(Lbeta), matrix(0, p, q))
  b212= cbind(matrix(0, q, p), matrix(0, q, q))
  B2=rbind(b211,b212)  # parameter delta

  b311 =cbind(matrix(0, p, p), matrix(0, p, q))
  b312= cbind(matrix(0, q, p), matrix(0, q, q))
  B3=rbind(b311,b312)  # parameter theta

  mu <- model$mu.fv
  sigma <- model$sigma.fv
  eta <- linkfun(mu)
  ai <- mu.eta(eta)
  dmu <- (-1/(2*mu)) + sigma/((y*sigma) + y + (sigma*mu)) +  ((sigma+1)*y)/(4*(mu^2)) - (sigma^2)/(4*y*(sigma+1)) #ok!
  Deltamu <- crossprod(x,diag(ai*dmu))


  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  dsigma <- (y+ mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (sigma*(sigma+2)*mu)/(4*(sigma+1)*(sigma+1)*y) + sigma/(2*(sigma+1))
  Deltasigma <- crossprod(z,diag(bi*dsigma))

  Delta <- rbind(Deltamu,Deltasigma)

  ##################theta#########################
  BT<-B(Delta,solve(h0),B3)
  autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
  vetorpcthetaPC<- eigen(BT,symmetric=TRUE)$vec[,1]
  dmaxG.theta<-abs(vetorpcthetaPC)
  vCithetaPC<-2*abs(diag(BT))
  Cb0<-vCithetaPC
  Cb.theta<-Cb0/sum(Cb0)
  ######################betas########################
  BM<-B(Delta,solve(h0),B1)
  autovmaxbetaPC<-eigen(BM,symmetric=TRUE)$val[1]
  vetorpcbetaPC<-eigen(BM,symmetric=TRUE)$vec[,1]
  dmaxG.beta<-abs(vetorpcbetaPC)
  vCibetaPC<-2*abs(diag(BM))
  Cb1<-vCibetaPC
  Cb.beta<-Cb1/sum(Cb1)
  ####################alphas#########################
  BD<-B(Delta,solve(h0),B2)
  autovmaxdeltaPC<-eigen(BD,symmetric=TRUE)$val[1]
  vetordeltaPC<-eigen(BD,symmetric=TRUE)$vec[,1]
  dmaxG.alpha=abs(vetordeltaPC)
  vCideltaPC=2*abs(diag(BD))
  Cb2=vCideltaPC
  Cb.alpha=Cb2/sum(Cb2)

  result <- list(dmax.beta = dmaxG.beta,
                 dmax.alpha = dmaxG.alpha,
                 dmax.theta = dmaxG.theta,
                 Ci.beta = Cb.beta,
                 Ci.alpha = Cb.alpha,
                 Ci.theta = Cb.theta)
return(result)
}
