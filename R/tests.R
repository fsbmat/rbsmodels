#'Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage lr.test(modelh0,modelh1)
#'grad.test(modelh0,modelh1)
#'wald.test(modelh1)
#'score.test(modelh0,modelh1)
#'
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples
#'
#'##
#'library(alr3)
#'data(landrent)
#'attach(landrent)
#'resp <- I(Y/X1)
#'y1 <-  split(resp, X4)$"1"
#'x21 <-  split(X2, X4)$"1"
#'##Fixed Precision
#'fit0 <- gamlss(y1 ~ x21, family=RBS(mu.link="identity"),method=CG()  )
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'#Precision Test
#'lr.test(fit0,fit1)
#'score.test(fit0,fit1)
#'grad.test(fit0,fit1)
#'wald.test(fit1)
#'@export

lr.test <- function(modelh0,modelh1)
{
  log.like = function(y,mu,delta)
  {
    lmd.i = (1/2)*delta - (1/2)*log((delta+1)) - (1/2)*log(mu) - (3/2)*log(y) + log(((delta*y) + y + (delta*mu))) - (y*(delta+1))/(4*mu) - (delta*delta*mu)/(4*y*(delta+1)) - (1/2)*log(16*pi)

    log.theta = sum(lmd.i)
    return(log.theta)
  }

  METHOD = "Likelihood Ratio test"
  DNAME  = deparse(substitute(modelh0))
  DNAME  = paste(DNAME, "vs", deparse(substitute(modelh1)))
  y=modelh0$y
  mu0 = modelh0$mu.fv
  sigma0 = modelh0$sigma.fv
  lH0 = log.like(y,mu0,sigma0)
  mu1 = modelh1$mu.fv
  sigma1 = modelh1$sigma.fv
  lH1 = log.like(y,mu1,sigma1)
  LR=2*(lH1 - lH0)
  q = modelh1$sigma.df
  gl = q-1
  PVAL = pchisq(LR,gl,lower.tail= F)
  names(gl) = "df"
  names(LR) = "LR"
  RVAL <- list(statistic = LR, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage lr.test(modelh0,modelh1)
#'grad.test(modelh0,modelh1)
#'wald.test(modelh0,modelh1)
#'score.test(modelh0,modelh1)
#'
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto
#'
#'@export

grad.test=function(modelh0,modelh1)
{

  UalphaH0.gam = function(modelh0,modelh1)
  {
    phi_linkstr = modelh1$sigma.link  #link function precision
    phi_linkobj = make.link(phi_linkstr)
    phi_mu.eta  = phi_linkobj$mu.eta # h'(delta)
    tau = phi_linkobj$linkfun
    tauH0 = tau(modelh0$sigma.fv)
    muH0 = modelh0$mu.fv # mu under H0
    deltaH0 = modelh0$sigma.fv  #delta under H0
    z = modelh1$sigma.x # matrix Z
    vt = modelh0$y # response y
    dH0 = 0.5 - (0.5/(deltaH0+1)) + ((vt+muH0)/((deltaH0*vt) + vt + (deltaH0*muH0))) - (vt/(4*muH0)) - ((deltaH0*(deltaH0+2)*muH0)/( 4*(((deltaH0+1)^2)*vt)))
    rval = dH0*phi_mu.eta(tauH0)*z
    colSums(rval)
  }

  METHOD = "Gradient test"
  DNAME  = deparse(substitute(modelh0))
  DNAME  = paste(DNAME, "vs", deparse(substitute(modelh1)))
  Ua = UalphaH0.gam(modelh0,modelh1)[-1]  #sem o intercepto
  alphah1 = modelh1$sigma.coefficients[-1]  #sem o intercepto
  G = t(Ua)%*%alphah1
  q = modelh1$sigma.df
  gl = q-1
  PVAL = pchisq(G,q-1,lower.tail= F)
  names(gl) = "df"
  names(G) = "G"
  RVAL <- list(statistic = G, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage lr.test(modelh0,modelh1)
#'grad.test(modelh0,modelh1)
#'wald.test(modelh0,modelh1)
#'score.test(modelh0,modelh1)
#'
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto
#'
#'@export
#'
wald.test=function(modelh1) #ok! est? implementado corretamente
{

  vcov.gam=function(model)
  {

    mu = model$mu.fv
    sigma = delta = model$sigma.fv

    x=model$mu.x
    z=model$sigma.x

    linkstr = model$mu.link
    linkobj = make.link(linkstr)
    mu.eta = linkobj$mu.eta
    eta = linkobj$linkfun
    eta = eta(model$mu.fv)

    phi_linkstr = model$sigma.link
    phi_linkobj = make.link(phi_linkstr)
    phi_mu.eta = phi_linkobj$mu.eta
    tau = phi_linkobj$linkfun
    tau = tau(model$sigma.fv)

    integral=function(aest)
    {
      fu=function(u)
      {
        w1 = (1 / ((1 +u^2)*(u^2)))
        w2 = (exp((-1 /(2*aest^2) )*((u - 1/u)^2)))
        (w1*w2)
      }
      return(integrate(fu,0,Inf)$value)
    }


    const = function(alpha,beta)
    {
      const = 1/(alpha*beta*beta*sqrt(2*pi))
      return(const)
    }

    esp = function(mu=1,sigma=1)
    {
      alpha = sqrt(2/sigma)
      beta = (mu*sigma)/(sigma+1)
      e = const(alpha,beta)*integral(alpha)
      return(e)
    }

    intemod =  mapply(esp,mu,sigma)


    v = (delta/(2*(mu^2))  + ((delta/(delta+1))^2)*intemod)*(mu.eta(eta)^2)
    s = (1/(2*mu*(delta+1)) + ((delta*mu)/((delta+1)^3))*intemod)*	mu.eta(eta)*phi_mu.eta(tau)
    u= ((((delta^2) +(3*delta)+1)/(2*((delta+1)^2)*(delta^2))) + (((mu^2)/((delta+1)^4))*intemod))*(phi_mu.eta(tau)^2)

    kbb = crossprod(v*x,x)
    kaa = crossprod(u*z,z)
    kba = crossprod(s*x,z)

    fisher = cbind(rbind(kbb, t(kba)), rbind(kba, kaa))

    vcov = solve(fisher)

    return(vcov)
  }


  METHOD = "Wald test"
  DNAME  = deparse(substitute(modelh1))
  p=modelh1$mu.df
  p1=p+1
  vcov = vcov.gam(modelh1)
  varalpha = vcov[-(1:p1),-(1:p1)]
  alphah1 = modelh1$sigma.coefficients[-1]
  W = t(alphah1)%*%solve(varalpha)%*%alphah1
  q = modelh1$sigma.df
  gl = q-1
  PVAL = pchisq(W,q-1,lower.tail= F)
  names(gl) = "df"
  names(W) = "W"
  RVAL <- list(statistic = W, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'Precision test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage lr.test(modelh0,modelh1)
#'grad.test(modelh0,modelh1)
#'wald.test(modelh0,modelh1)
#'score.test(modelh0,modelh1)
#'
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto
#'
#'@export

score.test=function(modelH0,modelH1)
{

  UalphaH0.gam = function(modelh0,modelh1)
  {
    phi_linkstr = modelh1$sigma.link  #link function precision
    phi_linkobj = make.link(phi_linkstr)
    phi_mu.eta  = phi_linkobj$mu.eta # h'(delta)
    tau = phi_linkobj$linkfun
    tauH0 = tau(modelh0$sigma.fv)
    muH0 = modelh0$mu.fv # mu under H0
    deltaH0 = modelh0$sigma.fv  #delta under H0
    z = modelh1$sigma.x # matrix Z
    vt = modelh0$y # response y
    dH0 = 0.5 - (0.5/(deltaH0+1)) + ((vt+muH0)/((deltaH0*vt) + vt + (deltaH0*muH0))) - (vt/(4*muH0)) - ((deltaH0*(deltaH0+2)*muH0)/( 4*(((deltaH0+1)^2)*vt)))
    rval = dH0*phi_mu.eta(tauH0)*z
    colSums(rval)
  }




  vcovH0=function(modelh0,modelh1)
  {

    mu = modelh0$mu.fv #estimativa mu sob H0
    delta = sigma = modelh0$sigma.fv  #estimativa delta sob H0

    x=modelh1$mu.x        #matriz X
    z=modelh1$sigma.x        #matriz Z

    linkstr = modelh0$mu.link
    linkobj = make.link(linkstr)
    mu.eta = linkobj$mu.eta
    eta = linkobj$linkfun
    etaH0 = eta(modelh0$mu.fv)

    phi_linkstr = modelh1$sigma.link
    phi_linkobj = make.link(phi_linkstr)
    phi_mu.eta = phi_linkobj$mu.eta
    tau = phi_linkobj$linkfun
    tauH0 = tau(modelh0$sigma.fv)

    #integral I(delta_i)
    integral=function(aest)
    {
      fu=function(u)
      {
        w1 = (1 / ((1 +u^2)*(u^2)))
        w2 = (exp((-1 /(2*aest^2) )*((u - 1/u)^2)))
        (w1*w2)
      }
      return(integrate(fu,0,Inf)$value)
    }

    const = function(alpha,beta)
    {
      const = 1/(alpha*beta*beta*sqrt(2*pi))
      return(const)
    }

    esp = function(mu=1,sigma=1)
    {
      alpha = sqrt(2/sigma)
      beta = (mu*sigma)/(sigma+1)
      e = const(alpha,beta)*integral(alpha)
      return(e)
    }
    intemod =  mapply(esp,mu,sigma)

    v = (delta/(2*(mu^2))  + ((delta/(delta+1))^2)*intemod)*(mu.eta(etaH0)^2)
    s = (1/(2*mu*(delta+1)) + ((delta*mu)/((delta+1)^3))*intemod)*mu.eta(etaH0)*phi_mu.eta(tauH0)
    u= ((((delta^2) +(3*delta)+1)/(2*((delta+1)^2)*(delta^2))) + (((mu^2)/((delta+1)^4))*intemod))*(phi_mu.eta(tauH0)^2)

    kbb = crossprod(v*x,x)
    kaa = crossprod(u*z,z)
    kba = crossprod(s*x,z)

    fisher = cbind(rbind(kbb, t(kba)), rbind(kba, kaa)) #matriz de informa?ao de fisher

    vcov = solve(fisher) # matriz de covari?ncias

    return(vcov)
  }


  METHOD = "Rao score test"
  DNAME  = deparse(substitute(modelh0))
  DNAME  = paste(DNAME, "vs", deparse(substitute(modelh1)))
  p = modelH0$mu.df
  p1= p+1
  varalpha = vcovH0(modelH0,modelH1)[-(1:p1),-(1:p1)]
  Ua.= UalphaH0.gam(modelH0,modelH1)[-1]
  SC = t(Ua.)%*%varalpha%*%Ua.
  q = modelH1$sigma.df
  gl = q-1
  PVAL = pchisq(SC,q-1,lower.tail= F)
  names(gl) = "df"
  names(SC) = "SC"
  RVAL <- list(statistic = SC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#'Goodness of Fit Tests
#'
#'@description Performs the Anderson-Darlin and Cramer-von Mises tests
#'
#'@usage ad.testg(x,cdf)
#'cvm.testg(x,cdf)
#'
#'
#' @param x a numeric vector of data values, the number of which must be greater than 7. Missing values are allowed.
#' @param cdf cumlative distribution function.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto
#'
#'@references
#'Stephens, M.A. (1986): Tests based on EDF statistics. In: D'Agostino, R.B. and Stephens, M.A., eds.: Goodness-of-Fit Techniques. Marcel Dekker, New York.
#'
#'Thode Jr., H.C. (2002): Testing for Normality. Marcel Dekker, New York.
#'
#'@examples
#'x<- rRBS(1000)
#'fit <- gamlss(x1~1,family=BS(),method=CG())
#'mu<- fit$mu.coefficients ; mu
#' sigma <- fit$sigma.coefficients ; sigma
#' cdf <- function(x) pBS(x,mu,sigma)
#' ad.testg(x,cdf)
#' cvm.testg(x,cdf)
#'@export
ad.testg = function(x,cdf)
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  u <- sapply(x,cdf)
  y <- qnorm(u)
  test <- ad.test(y)
  stat <- test$statistic
  pval<-test$p.value
  RVAL <- list(statistic = stat, p.value = pval, method = "Anderson-Darling test for F distribution",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

#'Goodness of Fit Tests
#'
#'@description Performs the Anderson-Darlin and Cramer-von Mises tests
#'
#'@usage ad.testg(x,cdf)
#'cvm.testg(x,cdf)
#'
#'
#' @param x a numeric vector of data values, the number of which must be greater than 7. Missing values are allowed.
#' @param cdf cumlative distribution function.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto
#'
#'@references
#'Stephens, M.A. (1986): Tests based on EDF statistics. In: D'Agostino, R.B. and Stephens, M.A., eds.: Goodness-of-Fit Techniques. Marcel Dekker, New York.
#'
#'Thode Jr., H.C. (2002): Testing for Normality. Marcel Dekker, New York.
#'
#'@examples
#'x<- rRBS(1000)
#'fit <- gamlss(x~1,family=RBS(),method=CG())
#'mu<- fit$mu.coefficients ; mu
#'sigma <- fit$sigma.coefficients ; sigma
#' cdf <- function(x) pRBS(x,mu,sigma)
#' ad.testg(x,cdf)
#' cvm.testg(x,cdf)
#'@export

cvm.testg = function(x,cdf)
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  u <- sapply(x,cdf)
  y <- qnorm(u)
  test <- cvm.test(y)
  stat <- test$statistic
  pval<-test$p.value
  RVAL <- list(statistic = stat, p.value = pval, method = "Cramer-von Mises test for F distribution",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
