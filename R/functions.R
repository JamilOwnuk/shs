##########################################################################################
#SKEW HYPERBOLIC SECANT TYPE I and II DISTRIBUTION (d,p,q and r functions)
##########################################################################################
#' @export
dshsI = function(y, mu = 0, sigma = 1, nu = 1)
{
  if(length(y) == 0) stop("y must be provided.")
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(nu <= 0) stop("nu must be a positive number.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  sapply(X=y, FUN=function(y, mu, sigma, nu){(2*gamma((nu+1)/2)/(sigma*sqrt(nu*pi)*gamma(nu/2)))*exp((y - (mu ))/sigma)*(1+exp(2*((y - (mu))/sigma))/nu)^-((nu+1)/2)}, mu=mu,sigma=sigma,nu=nu)
}
##########################################################################################
#' @export
dshsII = function(y, mu = 0, sigma = 1, alpha = 0)
{
  if(length(y) == 0) stop("y must be provided.")
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(abs(alpha) >= 1 ) stop("alpha must be between -1 and 1.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  sapply(X=y, FUN=function(y, mu, sigma, alpha){cos((pi/2)*alpha) * pi ^-1 * sigma ^ -1 * exp(alpha * ((y-mu)/sigma)) * cosh((y - mu)/sigma) ^ -1}, mu=mu,sigma=sigma,alpha=alpha)
}

##########################################################################################
#' @export
pshsI = function(q, mu = 0, sigma = 1, nu = 1,lower.tail=TRUE)
{
  if(length(q) == 0) stop("q must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(nu <= 0) stop("nu must be a positive number.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")

  cdf = pbeta(((exp(2*(q-mu)/sigma))/((nu)+exp(2*(q-mu)/sigma))), 0.5, nu/2)
  ifelse(test=lower.tail == TRUE,yes=return(cdf),no=return(1-cdf))
}

##########################################################################################
#' @export
pshsII = function(q, mu = 0, sigma = 1, alpha = 0,lower.tail=TRUE)
{
  if(length(q) == 0) stop("q must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(abs(alpha) >= 1 ) stop("alpha must be between -1 and 1.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")

  cdf = pbeta(exp(2*q)/(1+exp(2*q)), (1+alpha)/2, (1-alpha)/2)
  ifelse(test=lower.tail == TRUE,yes=return(cdf),no=return(1-cdf))
}

##########################################################################################
#' @export
qshsI = function(prob,mu=0,sigma=1,nu=1,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("prob must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(nu <= 0) stop("nu must be a positive number.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in [0,1].")

  q = 0.5*log((qbeta(prob, 0.5, nu/2, lower.tail = lower.tail) *nu)/(1-qbeta(prob, 0.5, nu/2, lower.tail = lower.tail)))
  return(sigma*q+mu)
}

##########################################################################################
#' @export
qshsII = function(prob,mu=0,sigma=1,alpha=0,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("prob must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(abs(alpha) >= 1 ) stop("alpha must be between -1 and 1.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in [0,1].")

  q = 0.5*log((qbeta(prob, (1+alpha)/2, (1-alpha)/2, lower.tail = lower.tail))/(1-qbeta(prob, (1+alpha)/2, (1-alpha)/2, lower.tail = lower.tail)))
  return(sigma*q+mu)
}

##########################################################################################
#' @export
rshsI = function(n, mu = 0, sigma = 1, nu = 1)
{
  if(length(n) == 0) stop("The sample size n must be provided.")
  if(n <= 0 || n%%1 !=0) stop("The sample size n must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(nu <= 0) stop("nu must be a positive number.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  library(truncdist)
  rnghs <- n
  rtt   <- rtrunc(n, spec = "t", 0, Inf,nu)
  r     <- sigma*log(rtt)+mu
  return(r)
}

##########################################################################################
#' @export
rshsII = function(n, mu = 0, sigma = 1, alpha= 0)
{
  if(length(n) == 0) stop("The sample size n must be provided.")
  if(n <= 0 || n%%1 !=0) stop("The sample size n must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(abs(alpha) >= 1 ) stop("alpha must be between -1 and 1.")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  r <- rep(0,n)
  for (i in 1:n) {
    R  <- rbeta(1, (1 + alpha)/2 , (1 - alpha)/2)
    if(R == 1){

    }
    else {
      r[i] <- mu + sigma * 0.5 * log(R/(1 - R))
    }
  }
  return(r)
}

##########################################################################################
#' ML fit of Linear regression model with skew heperbolic secant type I distributed errors by using \code{optim} function.
#'
#' ML fit of Linear regression model with SHS type I distributed errors by using \code{optim} function.
#' @param y response vector.
#' @param x design matrix.
#' @param initial vector of initial values.
#' @param method The method to be used. See \code{optim} for Detailes.
#' @param lower,upper Bounds on the variables for the "L-BFGS-B" method, or bounds in which to search for method "Brent".
#' @param control A list of control parameters. See \code{optim} for Detailes.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @return \item{coffitients}{estimation of coffitients.}
#' @return \item{sigma}{estimation of sigma.}
#' @return \item{nu}{estimation of nu.}
#' @return \item{loglik}{value of log-likelihood function at estimated parameters.}
#' @return \item{hessian}{estimated hessian matrix.}
#' @return \item{counts}{See \link{optim} for Detailes.}
#' @return \item{convergents}{See \link{optim} for Detailes.}
#' @author Jamil Ownuk <\email{JamilOwnuk@yahoo.com}>, Hossein Baghishani <\email{hbaghishani@yahoo.com}> and Ahmad Nezakati <\email{anezakati591@yahoo.com}>
#' @examples
#' beta = c(1, 2, 3); sigma = 2; nu = 6; n = 10000;
#' e    = rshsI(n,mu = 0,sigma = sigma,nu = nu)  #A random sample
#' x1   = rnorm(n); x2 = rnorm(n); #covariates
#' y    = beta[1] + beta[2] * x1 + beta[3] * x2 + e
#' regfit = mleshsIreg(y, x=cbind(rep(1, length(y)), x1, x2),
#'                     initial = c(beta, sigma, nu), method = "Nelder-Mead")
#' cbind(c(beta, sigma, nu), c(regfit$coefficient, regfit$sigma, regfit$nu))
#' @export
mleshsIreg <- function(y, x=NULL, initial, method = "BFGS",
                       lower = -Inf, upper = Inf,
                       control = list(), hessian = FALSE){
  if(length(y) == 0) stop("y must be provided.")
  if(is.null(x)){
    x=matrix(rep(1,n), ncol = 1, nrow = n)
  }
  n       <- length(y)
  p       <- ncol(x)
  LL1 <- function(b) {
    R=(2*gamma((b[p+2]+1)/2)/(b[p+1]*sqrt(b[p+2]*pi)*gamma(b[p+2]/2)))*
      exp((y - x %*% b[1:p])/b[p+1])*(1+exp(2*((y - x %*% b[1:p])/b[p+1]))/b[p+2])^
      -((b[p+2]+1)/2)
    #
    -sum(log(R))
  }

  LL2 <- function(b) {
    w1      <- -1 + ((b[p+2]+1)/b[p+2]) * (exp(2*(y - x %*% b[1:p])/(b[p+1]))) /
      (1 + b[p+2]^(-1) * exp(2*(y - x %*% b[1:p])/(b[p+1])))
    w2      <- ((b[p+2]+1)/b[p+2]) * (exp(2*(y - x %*% b[1:p])/(b[p+1]))) /
      (1 + b[p+2]^(-1) * exp(2*(y - x %*% b[1:p])/(b[p+1])))
    elbeta  <- 1/b[p+1] * (t(x) %*% w1)
    elsigma <- 1/b[p+1] * (t((y - x %*% b[1:p])/(b[p+1])) %*% w2) - 1/b[p+1] *
      rep(1, n) %*% ((y - x %*% b[1:p])/(b[p+1])) - n/b[p+1]
    elnu    <- 1/(2*b[p+2]) * (rep(1, n) %*% w2) - 0.5 *
      rep(1, n) %*% log(1 + b[p+2]^(-1) * exp(2*(y - x %*% b[1:p])/(b[p+1]))) -
      n * psigamma(b[p+2]/2)/2 + n * psigamma((b[p+2]+1)/2)/2 - n/(2*b[p+2])
    -c(elbeta, elsigma, elnu)
  }
  fit <- optim(par=initial, fn=LL1, gr = LL2, method = method,
               lower = lower, upper = upper,
               control = control, hessian = hessian)
  coefficient = fit$par[1:p]; sigma = fit$par[p+1]; nu = fit$par[p+2];
  loglik = -fit$value; hessian = fit$hessian; counts = fit$counts;
  convergence = fit$convergence;
  return(list(coefficient=coefficient, sigma=sigma, nu=nu, loglik=loglik,
              hessian=hessian, counts=counts, convergence=convergence))
}
###########################################################################################
#' ML fit of Linear regression model with skew heperbolic secant type II distributed errors by using \code{optim} function.
#'
#' ML fit of Linear regression model with SHS type II distributed errors by using \code{optim} function.
#' @param y response vector.
#' @param x design matrix.
#' @param initial vector of initial values.
#' @param method The method to be used. See \optim for Detailes.
#' @param lower,upper Bounds on the variables for the "L-BFGS-B" method, or bounds in which to search for method "Brent".
#' @param control A list of control parameters. See \optim for Detailes.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @return \item{coffitients}{estimation of coffitients.}
#' @return \item{sigma}{estimation of sigma.}
#' @return \item{alpha}{estimation of alpha}
#' @return \item{loglik}{value of log-likelihood function at estimated parameters.}
#' @return \item{hessian}{estimated hessian matrix.}
#' @return \item{counts}{See \link{optim} for Detailes.}
#' @return \item{convergents}{See \link{optim} for Detailes.}
#' @author Jamil Ownuk <\email{JamilOwnuk@yahoo.com}>, Hossein Baghishani <\email{hbaghishani@yahoo.com}> and Ahmad Nezakati <\email{anezakati591@yahoo.com}>
#' @examples
#' beta = c(1, 2, 3); sigma = 2; alpha = 0.3; n = 10000;
#' e    = rshsII(n,mu = 0,sigma = sigma,alpha = alpha)  #A random sample
#' x1   = rnorm(n); x2 = rnorm(n); #covariates
#' y    = beta[1] + beta[2] * x1 + beta[3] * x2 + e
#' regfit = mleshsIIreg(y, x=cbind(rep(1, length(y)), x1, x2),
#'                     initial = c(beta, sigma, alpha), method = "Nelder-Mead")
#' cbind(c(beta, sigma, alpha), c(regfit$coefficient, regfit$sigma, regfit$alpha))
#' @export
mleshsIIreg <- function(y, x=NULL, initial, gr = NULL, method = "BFGS",
                       lower = -Inf, upper = Inf,
                       control = list(), hessian = FALSE){
  if(length(y) == 0) stop("y must be provided.")
  if(is.null(x)){
    x=matrix(rep(1,n), ncol = 1, nrow = n)
  }
  n       <- length(y)
  p       <- ncol(x)
  LL1 <- function(b) {
    R=cos((pi/2)*b[p+2]) * pi ^-1 * b[p+1] ^ -1 * exp(b[p+2] * ((y - x %*% b[1:p])/b[p+1])) * cosh((y - x %*% b[1:p])/b[p+1]) ^ -1
    #
    -sum(log(R))
  }

  fit <- optim(par=initial, fn=LL1, gr = NULL, method = method,
               lower = lower, upper = upper,
               control = control, hessian = hessian)
  coefficient = fit$par[1:p]; sigma = fit$par[p+1]; alpha = fit$par[p+2];
  loglik = -fit$value; hessian = fit$hessian; counts = fit$counts;
  convergence = fit$convergence;
  return(list(coefficient=coefficient, sigma=sigma, alpha=alpha, loglik=loglik,
              hessian=hessian, counts=counts, convergence=convergence))
}
