# ---------------------------------------------------------------------------------------
# Reparameterized Birnbaum-Saunders distribution (RBS)
# created by Leiva, V., Santos-Neto, M., Cysneiros, F.J.A. and Barros, M.
# manoel.ferreira@ufcg.edu.br
# ---------------------------------------------------------------------------------------

Birnbaum.Saunders = BS = function (mu.link = "identity" , sigma.link="identity")
{
  mstats = checklink("mu.link", "BS", substitute(mu.link),c("sqrt","log","identity"))
  dstats = checklink("sigma.link", "BS", substitute(sigma.link),c("sqrt", "log", "identity"))

  structure(list(family = c("BS","Birnbaum.Saunders"),
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


                 G.dev.incr = function(y,mu,sigma,...) -2*dBS(y,mu,sigma,log=TRUE),

                 rqres = expression(resbs(y=y,mu=mu,sigma=sigma)),

                 mu.initial = expression({mu = mean(y)}),
                 sigma.initial = expression({sigma = rep(sigmatil(y),length(y))}),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0)),
            class = c("gamlss.family","family"))
}

sigmatil=function(y)
{
  s = mean(y)
  r = 1/mean(1/y)
  alphatil = (2*( (s/r)^(1/2)  - 1))^(1/2)
  dest = 2/(alphatil^2 )
  return(dest)
}


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



########## Density  function of Birnbaum-Saunders ##########
dBS<-function(x, mu=1, sigma=1, log=FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))  stop(paste("x must be positive", "\n", ""))
  log.lik =  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1)))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}

########## Acumulate  function distribution of Birnbaum-Saunders ##########
pBS <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
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
        dbs(x, alpha = a[i], beta = b[i], log=FALSE), 0.001, q[i] )$value
    }
  }

  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}


qBS = function (p, mu = 0.5, sigma = 1, lower.tail = TRUE,
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

  suppressWarnings(q <- qbs(p ,alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE))
  q
}

rBS = function(n, mu=1, sigma=1)
{
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  r = rbs(n, alpha = a, beta = b)
  r
}


resbs=function(y,mu,sigma)
{
  Ims = mapply(esp,mu,sigma)
  z = -1/(2*mu) - (sigma^2)/(4*(sigma+1)*y) + ((sigma+1)*y)/(4*mu*mu) + sigma/((sigma*y) + y + (sigma*mu))
  v = sigma/(2*mu*mu) + ((sigma*sigma)/((sigma+1)*(sigma+1)))*Ims
  res = z/sqrt(v)
  return(res)
}


########## plot the function density of Birnbaum-Saunders ##########
plotBS = function (mu = .5, sigma = 1, from = 0, to = 0.999, n = 101, title="title",
                   ...)
{
  y = seq(from = 0.001, to = to, length.out = n)
  pdf = dBS(y, mu = mu, sigma = sigma)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf,pr0)), type = "l",lwd=3)

}


#----------------------------------------------------------------------------------------
#calculates the expected value of the response for a Birnbaum-Saunders fitted model
meanBS = function (obj)
{
  if (obj$family[1] != "BS")
    stop("the object do not have a BS distribution")
  meanofY = fitted(obj, "mu")
  meanofY
}




