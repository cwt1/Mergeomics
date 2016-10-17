\name{tool.normalize}
\alias{tool.normalize}
\title{
Estimate statistical scores based on Gauss distribution
}
\description{
To estimate the both pre-liminary and final p-values, 
\code{tool.normalize} normalizes the given data, \code{x}, based on 
Gaussian distribution defined by \code{prm} if it is provided. If
\code{prm} is not provided \code{tool.normalize} utilizes the mean and 
std dev of \code{x}.
}
\usage{
tool.normalize(x, prm = NULL, inverse = FALSE)
}
\arguments{
\item{x}{data that is aimed to be normalized and produced by a simulation
process}
\item{prm}{normalization will take place according to the specified 
Gaussian distribution parameters, i.e. mean and std dev. If it is not 
specified, Gaussian statistics of \code{x} will be obtained and utilized}
\item{inverse}{specifies whether the normalization takes place in reverse
order}
}
\value{
\item{prm }{transformed (normalized) parameters for either enrichment 
score or p-values}
}
\examples{
set.seed(1)
## let us assume we have a set of simulated enrichment scores and 
## one observed score
x <- rnorm(10) ## obtained from 1st permutation test
obs <- rnorm(1)
## Estimate preliminary P-value:
param <- tool.normalize(x)
z <- tool.normalize(obs, param)
p <- pnorm(z, lower.tail=FALSE)

## Estimate final P-value.
y <- rnorm(10) ## obtained from 2nd permutation test
param <- tool.normalize(c(x, y))
z <- tool.normalize(obs, param)
p <- pnorm(z, lower.tail=FALSE)
p <- max(p, .Machine$double.xmin)
}
\author{
Ville-Petteri Makinen 
}
