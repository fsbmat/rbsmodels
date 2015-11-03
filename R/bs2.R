#------------------------------------------------------------------------------
# BS PDF, CDF, quantiles, random numbers
#------------------------------------------------------------------------------
# PDF
dbs <- function(x, alpha = 1, beta = 1, log = FALSE){
  if(!is.numeric(x)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- x
  c   <-(1 / sqrt(2 * pi))
  u   <- (alpha^(-2)) * ((x / beta) + (beta / x) - 2)
  e   <- exp((-1 / 2 ) * u)
  du  <- (x^(-3 / 2) * (x + beta)) /
    (2 * alpha * sqrt (beta))
  pdf <- c * e * du
  if (log==TRUE){pdf <-log(pdf)}
  return(pdf)
}

# CDF

pbs <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if(!is.numeric(q)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- q
  s   <- (x / beta)
  a   <- ((1 / alpha) * (s^(1 / 2) - s^(-1 / 2)))
  cdf <- pnorm(a, 0, 1)
  if (lower.tail == FALSE){cdf <-(1 - cdf)}
  if (log.p == TRUE){cdf <-log(cdf)}
  return(cdf)
}

# Quantile function

qbs <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  if (log.p == TRUE){p  <- log(p)}
  if (lower.tail == FALSE){p <- (1 - p)}
  q   <- beta * (((alpha * qnorm(p, 0, 1) / 2) +
                    sqrt(((alpha * qnorm(p, 0, 1) / 2)^2) +
                           1)))^2
  return(q)
}

# Random number generator

rbs <- function(n, alpha = 1, beta = 1){
  if (!is.numeric(n)||!is.numeric(alpha)||!is.numeric(beta))
  {stop("non-numeric argument to mathematical function")}
  if (n == 0){stop("value of n must be greater or equal then 0")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  z   <- rnorm(n, 0, 1)
  t   <- beta * (1 + (((alpha^2) * (z^2)) / 2) +
                   (alpha * z *sqrt((((alpha^2) * (z^2))/ 4)+ 1)))
  return(t)
}

