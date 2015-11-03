# ---------------------------------------------------------------------------------------
# The Zero-Inflated Reparameterized Birnbaum-Saunders distribution (ZIRBS)
# created by Leiva, V., Santos-Neto, M., Cysneiros, F.J.A. and Barros, M.
# manoel.ferreira@ufcg.edu.br
# ---------------------------------------------------------------------------------------
ZIRBS <- function (mu.link = "log", sigma.link = "log", nu.link = "logit")
{
  mstats = checklink("mu.link", "ZIRBS", substitute(mu.link),
                     c("sqrt", "log", "identity"))
  dstats = checklink("sigma.link", "ZIRBS", substitute(sigma.link),
                     c("sqrt", "log", "identity"))
  vstats = checklink("nu.link", "ZIRBS", substitute(nu.link),
                     c("logit", "probit", "cloglog", "log", "own"))

  structure(list(family = c("ZIRBS", "Zero-Inflated Reparameterized Birnbaum-Saunders"),
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
                 nopar = 3,
                 type = "Mixed",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link = as.character(substitute(nu.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 nu.linkfun = vstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 nu.linkinv = vstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 nu.dr = vstats$mu.eta,
                 dldm = function(y,mu,sigma) #first derivate of log-density respect to mu
                 {
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + (sigma/(sigma+1))*(1/( y + ((mu*sigma)/(sigma+1)) ))

                   dldm <- ifelse((y == 0) , 0, ystart - mustart)
                   dldm
                 },

                 d2ldm2 = function(y,mu,sigma) {        #expected of second derivate of log-density respect to mu
                   Ims = mapply(esp,mu,sigma)
                   d2ldm2 = ifelse((y == 0) , 0, - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims )
                   d2ldm2
                 },


                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma
                   sigmastart  = -(sigma+2)/(2*(sigma+1))
                   ystart.sigma   = (mu/((sigma+1)^2))*(1/( y + ((mu*sigma)/(sigma+1)) )) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = ifelse((y==0),0,(ystart.sigma - sigmastart))
                   dldd
                 },

                 d2ldd2 = function(y,mu,sigma) {      #expected of second derivate log-density respect to sigma
                   Ims = mapply(esp,mu,sigma)
                   lss =  ((sigma^2) + 3*sigma + 1)/(2*((sigma+1)^2)*(sigma^2))
                   d2ldd2 <- ifelse((y==0),0, -lss - ((mu^2)/((sigma+1)^4))*Ims)
                   d2ldd2
                 },

                 dldv = function(y,nu) {       #first derivate log-density respect to nu
                   dldv <- ifelse(y == 0, 1/nu, -1/(1 - nu))
                   dldv
                 },


                 d2ldv2 = function(nu) {         #expected of second derivate log-density respect to nu
                   d2ldv2 <- -1/(nu * (1 - nu))
                   d2ldv2
                 },


                 d2ldmdd = function(y,mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma
                   Ims = mapply(esp,mu,sigma)
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd <- ifelse((y == 0), 0, - lms - ((mu*sigma)/((sigma+1)^3))*Ims )
                   d2ldmdd
                 },

                 d2ldmdv = function(y) {  #partial derivate of log-density respect to mu and nu
                   d2ldmdv <- rep(0,length(y))
                   d2ldmdv
                 },

                 d2ldddv = function(y) {   #partial derivate of log-density respect to sigma and nu
                   d2ldddv <- rep(0,length(y))
                   d2ldddv
                 },


                 G.dev.incr = function(y, mu, sigma, nu, ...){  # Global deviance
                   -2 * dZIRBS(y, mu, sigma, nu, log = TRUE)
                 },

                 rqres = expression({     # (Normalize quantile) residuals
                   uval <- ifelse(y == 0, nu * runif(length(y), 0, 1),
                                  (1 - nu) * pZIRBS(y, mu, sigma, nu))
                   rqres <- qnorm(uval)
                 }),
                 mu.initial = expression(mu <- mean(y[y>0])),
                 sigma.initial = expression(sigma <- rep(1, length(y))),
                 nu.initial = expression(nu <- rep(0.3, length(y))),
                 mu.valid = function(mu) TRUE,
                 sigma.valid = function(sigma) all(sigma > 0),
                 nu.valid = function(nu) all(nu > 0) && all(nu < 1),
                 y.valid = function(y) all(y >= 0)),
            class = c("gamlss.family", "family"))

}


#-------------------------------------------------------------------------------
# Expected value de (T + mu*sigma/(sigma+1))^-2
esp <- function(mu=1,sigma=1)
{
  a <- sqrt(2/sigma)
  b <- (mu*sigma)/(sigma+1)
  E <- function(x)(1/((x+beta.)^2))*dbs(x,alpha=a,beta=b)
  return(integrate(E,0,Inf)$value)
}


#----------------------------------------------------------------------------------------
########## Density  function of ZIRBS ##########
dZIRBS<-function(x, mu=1, sigma=1, nu=.1, log=FALSE)
{
  if (any(mu < 0) || any(is.na(mu)))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0) || any(is.na(sigma)))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x < 0))  stop(paste("x must be positive", "\n", ""))
  log.lik <- ifelse(x==0, log(nu), log(1-nu) +  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1))))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}

#----------------------------------------------------------------------------------------
########## Acumulate  function distribution of ZIRBS ##########
pZIRBS <- function(q, mu=1, sigma=1, nu=0.1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  a <- sqrt(2/sigma)
  b <- (mu*sigma)/(sigma+1)
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
  cdf <- ifelse((q==0), nu, nu+(1-nu)*cdf)
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}

#----------------------------------------------------------------------------------------
########## Quantile function of ZIRBS ##########

qZIRBS <- function (p, mu = 0.5, sigma = 1, nu = 0.1, lower.tail = TRUE,
                    log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = delta
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = p
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  a <- sqrt(2/sigma)
  b <- (mu*sigma)/(sigma+1)

  suppressWarnings(q <- ifelse((nu >= p),0, qbs((p - nu)/(1-nu),alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE)))
  q
}


#----------------------------------------------------------------------------------------
######## Random generation function of ZIRBS########

rZIRBS <- function (n, mu = 0.5, sigma = 1, nu = 0.1)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = delta
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = p
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qZIRBS(p, mu = mu, sigma = sigma, nu = nu)
  r
}


#----------------------------------------------------------------------------------------
# plot the function density of ZIRBS ##########
plotZIRBS = function (mu = .5, sigma = 1, nu = 0.1, from = 0, to = 0.999, n = 101, title="title",
                      ...)
{
  y <- seq(from = 0.001, to = to, length.out = n)
  pdf <- dZIRBS(y, mu = mu, sigma = sigma, nu = nu)
  pr0 <- c(dZIRBS(0, mu = mu, sigma = sigma, nu = nu))
  po <- c(0)
  plot(pdf ~ y, main=title, ylim = c(0, max(pdf,pr0)), type = "l",lwd=3)
  points(po, pr0, type = "h",lwd=3)
  points(po, pr0, type = "p", col = "red",lwd=3)
}

#----------------------------------------------------------------------------------------
#calculates the expected value of the response for a ZIRBS fitted model
meanZIRBS <- function (obj)
{
  if (obj$family[1] != "ZIRBS")
    stop("the object do not have a ZIRBS distribution")
  meanofY <- (1 - fitted(obj, "nu")) * fitted(obj, "mu")
  meanofY
}


