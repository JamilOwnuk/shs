\name{SHSI}
\alias{SHSI}
\alias{dshsI}
\alias{pshsI}
\alias{qshsI}
\alias{rshsI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
skew heperbolic secant type I distribution
}
\description{
gives the density, cumulative, quantil function and random sample of shsI distribution.
}
\usage{
dshsI(y, mu = 0, sigma = 1, nu = 1)
pshsI(q, mu = 0, sigma = 1, nu = 1, lower.tail = TRUE)
qshsI(prob, mu = 0, sigma = 1, nu = 1, lower.tail = TRUE)
rshsI(n, mu = 0, sigma = 1, nu = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y, q}{a numeric vector.}
\item{prob}{vector of probabilities.}

\item{mu}{location parameter.}

\item{sigma}{scale parameter.}

\item{nu}{skewness parameter.}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].
}
}
\value{
dshsI, pshsI, qshsI,  gives the density, the distribution function and the quantile function, respectively and rshsI generates a random sample.
}
\author{
Jamil Ownuk <\email{JamilOwnuk@yahoo.com}>, Hossein Baghishani <\email{hbaghishani@yahoo.com}> and Ahmad Nezakati <\email{anezakati591@yahoo.com}>
}

\seealso{
\link{dshsII}, \link{pshsII}, \link{qshsII}, \link{rshsII}
}
\examples{
## Density fuction
library(shs)
sseq = seq(-10,10,0.01)
plot(sseq,dshsI(y=sseq,mu=1,sigma=2,nu=1),type = "l",lwd=2,col="black",
xlab="x",ylab="f(x)",ylim=c(0,0.25), main="SHS type I Density function")
lines(sseq, dshsI(y=sseq,mu=1,sigma=2,nu=3),type = "l",lwd=2,col="red")
lines(sseq, dshsI(y=sseq,mu=1,sigma=2,nu=0.75),type = "l",lwd=2,col="blue")

## Distribution fuction
sseq = seq(-10,10,0.01)
plot(sseq,pshsI(q=sseq, mu=1, sigma=2, nu=1),type="l",lwd=2,col="black",
xlab="x",ylab="F(x)", main="SHS type I Distribution function")
lines(sseq, pshsI(q=sseq,mu=1,sigma=2,nu=3),type = "l",lwd=2,col="red")
lines(sseq, pshsI(q=sseq,mu=1,sigma=2,nu=0.75),type = "l",lwd=2,col="blue")

## Quantile fuction
prob = seq(0,1,length.out = 1000)
plot(prob,qshsI(prob, mu=1, sigma=2, nu=1),type="l",lwd=2,col="black",
xlab="x",ylab=expression(F^{-1}~(x)), main="SHS type I Inverse Distribution function")
lines(prob, qshsI(prob,mu=1,sigma=2,nu=3),type = "l",lwd=2,col="red")
lines(prob, qshsI(prob,mu=1,sigma=2,nu=0.75),type = "l",lwd=2,col="blue")

## random sample generation
sample = rshsI(n=10000,mu=1,sigma=2,nu=2)
hist(sample, breaks=70, freq = FALSE, main="")
title(main="Histogram and True density of SHS type I")
sseq = seq(-20,15,0.01)
lines(sseq,dshsI(y=sseq,mu=1,sigma=2,nu=2),col="red",lwd=2)
}
