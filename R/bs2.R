#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
dcbs <- function(x, alpha = 1, beta = 1, log = FALSE){
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

#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
pcbs <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
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

#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'
qcbs <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  if (log.p == TRUE){p  <- log(p)}
  if (lower.tail == FALSE){p <- (1 - p)}
  q   <- beta * (((alpha * qnorm(p, 0, 1) / 2) +
                    sqrt(((alpha * qnorm(p, 0, 1) / 2)^2) +
                           1)))^2
  return(q)
}

#'The Classical Birnbaum-Saunders (BS) distribution
#'
#'@description Density, distribution function, quantile function and random generation
#'for the normal distribution with mean equal to \code{alpha} and standard deviation equal to \code{beta}.
#'
#'@usage dcbs(x, alpha = 1, beta = 1, log = FALSE)
#'pcbs(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'qcbs(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#'rcbs(n, alpha = 1, beta = 1)

#'
#' @param x,q vector of quantiles
#' @param alpha vector of scale parameter values
#' @param beta vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#'
#'
#' @details  Birnbaum and Saunders (1969) proposed the two-parameter Birnbaum-Saunders
#' distribution with density
#' \deqn{
#'  f_{T}(t) = \frac{1}{\sqrt{2\pi}} \exp\left[-\frac{1}{2\alpha^{2}}
#'                                             \left(\frac{t}{\beta}+\frac{\beta}{t}-2\right) \right]
#'  \frac{t^{-\frac{3}{2}} (t+\beta)}{2\alpha\sqrt{\beta}}; \ t>0,
#'  \alpha > 0, \beta > 0,}
#'as a failure time distribution for fatigue failure caused under cyclic loading. The parameters
#'alpha and beta are the shape and the scale parameters, respectively. In their derivation,
#'it was assumed that the failure is due to the development and growth of a
#' dominant crack.
#'
#'@return \code{dcbs} gives the density, \code{pcbs} gives the distribution function,
#'\code{qcbs} gives the quantile function, and \code{rcbs} generates random deviates.
#'
#'@references
#'Birnbaum, Z. W. and Saunders, S. C. (1969). A new family of life distributions. J. Appl. Probab. 6(2): 637-652.
#'
#'Chang D. S. and Tang, L. C. (1994). Random number generator for the Birnbaum-Saunders distribution. Computational and Industrial Engineering, 27(1-4):345-348.
#'
#'Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2006). Random number generators for the generalized Birnbaum-Saunders distribution. Submitted to Publication.
#'
#'Rieck, J. R. (2003). A comparison of two random number generators for the Birnbaum-Saunders distribution. Communications in Statistics - Theory and Methods, 32(5):929-934.
#'
#'@author
#'Víctor Leiva \email{victor.leiva@uv.cl}, Hugo Hernández \email{hugo.hernande@msn.com}, and Marco Riquelme \email{mriquelm@ucm.cl}.
#'
#'@examples
#'
#'## density for the Birnbaum-Saunders distribution
#'## with parameters alpha=0.5 y beta=1.0 in x=3.
#'dcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## cdf for the Birnbaum-Saunders distribution
#' ## with parameters alpha=0.5 y beta=1.0 in x=3.
#' pcbs(3,alpha=0.5,beta=1.0,log=FALSE)
#'
#' ## quantil function for p=0.5 in the Birnbaum-Saunders distribution
## with parameters alpha=0.5 y beta=1.0.
#' qcbs(0.5,alpha=0.5,beta=1.0,log=FALSE)
#'
#'## Examples for simulations
#'rcbs(n=6,alpha=0.5,beta=1.0)
#'sample<-rcbs(n=100,alpha=0.5,beta=1.0)
#'## Higtogram for sample
#'hist(sample)
#'@export
#'

rcbs <- function(n, alpha = 1, beta = 1)
  {
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

