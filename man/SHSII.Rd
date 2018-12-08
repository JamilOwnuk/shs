\name{SHSII}
\alias{SHSII}
\alias{dshsII}
\alias{pshsII}
\alias{qshsII}
\alias{rshsII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
skew heperbolic secant type II distribution
}
\description{
gives the density, cumulative, quantil function and random sample of shsII distribution.
}
\usage{
dshsII(y, mu = 0, sigma = 1, alpha = 0)
pshsII(q, mu = 0, sigma = 1, alpha = 0, lower.tail = TRUE)
qshsII(prob, mu = 0, sigma = 1, alpha = 0, lower.tail = TRUE)
rshsII(n, mu = 0, sigma = 1, alpha = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y, q}{a numeric vector.}
\item{prob}{vector of probabilities.}

\item{mu}{location parameter.}

\item{sigma}{scale parameter.}

\item{alpha}{skewness parameter.}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].
}
}
\value{
dshsII, pshsII, qshsII,  gives the density, the distribution function and the quantile function, respectively and rshsII generates a random sample.
}
\author{
Jamil Ownuk <\email{JamilOwnuk@yahoo.com}>, Hossein Baghishani <\email{hbaghishani@yahoo.com}> and Ahmad Nezakati <\email{anezakati591@yahoo.com}>
}

\seealso{
\link{dshsI}, \link{pshsI}, \link{qshsI}, \link{rshsI}
}
\examples{
## Density fuction
library(shs)
sseq = seq(-10,10,0.01)
plot(sseq,dshsII(y=sseq,mu=1,sigma=2,alpha=0),type = "l",lwd=2,col="black",
xlab="x",ylab="f(x)",ylim=c(0,0.2), main="SHS type II Density function")
lines(sseq, dshsII(y=sseq,mu=1,sigma=2,alpha=0.5),type = "l",lwd=2,col="red")
lines(sseq, dshsII(y=sseq,mu=1,sigma=2,alpha=-0.5),type = "l",lwd=2,col="blue")

## Distribution fuction
sseq = seq(-10,10,0.01)
plot(sseq,pshsII(q=sseq,mu=1,sigma=2,alpha=0),type = "l",lwd=2,col="black",
xlab="x",ylab="F(x)", main="SHS type II Distribution function")
lines(sseq, pshsII(q=sseq,mu=1,sigma=2,alpha=0.5),type = "l",lwd=2,col="red")
lines(sseq, pshsII(q=sseq,mu=1,sigma=2,alpha=-0.5),type = "l",lwd=2,col="blue")

## Quantile fuction
prob = seq(0,1,length.out = 1000)
plot(prob,qshsII(prob, mu=1, sigma=2, alpha=0),type="l",lwd=2,col="black", xlab="x"
,ylab=expression(F^{-1}~(x)), main="SHS type II Inverse Distribution function")
lines(prob, qshsII(prob,mu=1,sigma=2,alpha=0.5),type = "l",lwd=2,col="red")
lines(prob, qshsII(prob,mu=1,sigma=2,alpha=-0.5),type = "l",lwd=2,col="blue")

## random sample generation
sample = rshsII(n=10000,mu=1,sigma=2,alpha=0.25)
hist(sample, breaks=70, freq = FALSE, main="")
title(main="Histogram and True density of SHS type II")
sseq = seq(-15,25,0.01)
lines(sseq,dshsII(y=sseq,mu=1,sigma=2,alpha=0.25),col="red",lwd=2)
}
