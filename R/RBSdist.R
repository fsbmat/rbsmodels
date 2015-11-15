#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A, Barros, M. (2014)  Birnbaum-Saunders statistical modelling: a new approach. \emph{Statistical Modelling}, v. 14, p. 21-48, 2014.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#' ##install.packages(faraway)
#'library(faraway)
#'data(cpd)
#'attach(cpd)
#'model0 = gamlss(actual ~ projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model0)
#'model = gamlss(actual ~ 0+projected, family=RBS(mu.link="identity"),method=CG())
#'summary(model)
#'@export
#'

RBS <- function (mu.link = "identity" , sigma.link="identity")
{
  mstats = checklink("mu.link", "RBS", substitute(mu.link),c("sqrt","log","identity"))
  dstats = checklink("sigma.link", "RBS", substitute(sigma.link),c("sqrt", "log", "identity"))

  structure(list(family = c("RBS"),
                 parameters = list(mu=TRUE,sigma=TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 #the first derivative of the likelihood with respect to the location parameter mu
                 dldm = function(y,mu,sigma) #first derivate of log-density respect to mu
                 {
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + sigma/((sigma*y) + y + (sigma*mu))

                   dldm = ystart-mustart
                   dldm
                 },
                 #the expected second derivative of the likelihood with respect to the location parameter mu
                 d2ldm2 = function(y,mu,sigma) {        #expected of second derivate of log-density respect to mu
                   Ims = mapply(esp,mu,sigma)
                   d2ldm2 =  - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims
                   d2ldm2
                 },

                 #the first derivative of the likelihood with respect to the scale parameter sigma
                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma
                   sigmastart  = -(sigma)/(2*(sigma+1))
                   y2start   = (y+mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = y2start-sigmastart
                   dldd
                 },
                 #the expected second derivative of the likelihood with respect to the scale parameter sigma ok
                 d2ldd2 = function(y,mu,sigma) {      #expected of second derivate log-density respect to sigma
                   Ims = mapply(esp,mu,sigma)
                   lss =  ((sigma^2) + (3*sigma) + 1)/(2*sigma*sigma*((sigma+1)^2))
                   d2ldd2 = -lss - ((mu^2)/((sigma+1)^4))*Ims
                   d2ldd2
                 },
                 #the expected cross derivative of the likelihood with respect to both the location mu and scale parameter sigma
                 d2ldmdd = function(y,mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma
                   Ims = mapply(esp,mu,sigma)
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd = - lms - ((mu*sigma)/((sigma+1)^3))*Ims
                   d2ldmdd
                 },


                 G.dev.incr = function(y,mu,sigma,...) -2*dRBS(y,mu,sigma,log=TRUE),

                 rqres = expression(resrbs(y=y,mu=mu,sigma=sigma)),

                 mu.initial = expression({mu = mean(y)}),
                 sigma.initial = expression({sigma = rep(sigmatil(y),length(y))}),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0)),
            class = c("gamlss.family","family"))
}



#' sigmatil
#'
sigmatil=function(y)
{
  s = mean(y)
  r = 1/mean(1/y)
  alphatil = (2*( (s/r)^(1/2)  - 1))^(1/2)
  dest = 2/(alphatil^2 )
  return(dest)
}

#' esp
esp = function(mu=1,sigma=1)
{

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

  alpha = sqrt(2/sigma)
  beta = (mu*sigma)/(sigma+1)
  e = const(alpha,beta)*integral(alpha)
  return(e)
}

#'Residuals

resrbs=function(y,mu,sigma)
{
  Ims = mapply(esp,mu,sigma)
  z = -1/(2*mu) - (sigma^2)/(4*(sigma+1)*y) + ((sigma+1)*y)/(4*mu*mu) + sigma/((sigma*y) + y + (sigma*mu))
  v = sigma/(2*mu*mu) + ((sigma*sigma)/((sigma+1)*(sigma+1)))*Ims
  res = z/sqrt(v)
  return(res)
}


#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export


dRBS<-function(x, mu=1, sigma=1, log=FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))  stop(paste("x must be positive", "\n", ""))
  log.lik =  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1)))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}

#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export


pRBS <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{       if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  cdf1 <- pnorm((1/a)*(sqrt(q/b) - sqrt(b/q)))
  cdf <- cdf1

  ## the problem with this approximation is that it is not working with
  ## small sigmas and produce NA's. So here it is a solution
  if (any(is.na(cdf)))
  {
    index <- seq(along=q)[is.na(cdf)]
    for (i in index)
    {
      cdf[i] <- integrate(function(x)
        dcbs(x, alpha = a[i], beta = b[i], log=FALSE), 0.001, q[i] )$value
    }
  }

  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}

#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#'
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export

qRBS = function (p, mu = 0.5, sigma = 1, lower.tail = TRUE,
                log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = phi
    stop(paste("sigma must be positive", "\n", ""))


  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p <= 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)

  suppressWarnings(q <- qcbs(p ,alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE))
  q
}


#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export


rRBS = function(n, mu=1, sigma=1)
{
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  r = rcbs(n, alpha = a, beta = b)
  r
}


#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export
plotRBS = function(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, title="title", ...)
{
  y = seq(from = 0.001, to = to, length.out = n)
  pdf = dRBS(y, mu = mu, sigma = sigma)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf)), type = "l",lwd=3)

}

#'Reparameterized Birnbaum-Saunders (RBS) distribution for fitting a GAMLSS
#'
#'@description The fuction \code{RBS()} defines the BS distribution, a two paramenter
#'distribution, for a gamlss.family object to be used in GAMLSS fitting using using the
#'function \code{gamlss()}, with mean equal to the parameter \code{mu} and \code{sigma}
#'equal the precision parameter. The functions \code{dRBS}, \code{pRBS}, \code{qRBS} and
#'\code{rBS} define the density, distribution function, quantile function and random
#'genetation for the \code{RBS} parameterization of the RBS distribution.
#'
#'@usage RBS(mu.link = "identity", sigma.link = "identity")
#'dRBS(x, mu = 1, sigma = 1, log = FALSE)
#'pRBS(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'qRBS(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
#'rRBS(n, mu = 1, sigma = 1)
#'plotRBS(mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, ...)
#'meanRBS(obj)
#'
#' @param mu.link object for which the extraction of model residuals is meaningful.
#' @param sigma.link type of residual to be used.
#' @param x,q vector of quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param obj a fitted RBS object
#' @param ... other graphical parameters for plotting
#' @details  The parametrization of the normal distribution given in the function RBS() is
#'
#' \deqn{f_{Y}(y;\mu,\delta)=\frac{\exp\left(\delta/2\right)\sqrt{\delta+1}}{4\sqrt{\pi\mu}\,y^{3/2}}
#'\left[y+\frac{\delta \mu}{\delta+1}\right] \exp\left(-\frac{\delta}{4}
#'                                                     \left[\frac{y\{\delta+1\}}{\delta\mu}+\frac{\delta\mu}{y\{\delta+1\}}\right]\right) y>0.}
#'
#'@return returns a \code{gamlss.family} object which can be used to fit a normal distribution in the \code{gamlss()} function.
#'
#'@note For the function RBS(), mu is the mean and sigma is the precision parameter of the Birnbaum-Saunders distribution.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotRBS()
#'dat <- rRBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=RBS(),method=CG())
#'meanRBS(fit)
#'
#'@export

meanRBS = function (obj)
{
  if (obj$family[1] != "RBS")
    stop("the object do not have a RBS distribution")
  meanofY = fitted(obj, "mu")
  meanofY
}




