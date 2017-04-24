\name{tool.normalize.quality}
\alias{tool.normalize.quality}
\title{
Check normalization quality
}
\description{
\code{tool.normalize.quality} checks transformation quality by using 
Kolmogorov-Smirnov Test. It seeks the best log transform parameter within 
the previously specified upper and lower limits, and applies the log 
transform with the best log parameter.
}
\usage{
tool.normalize.quality(g, z)
}
\arguments{
\item{g}{normalization quality control will take place according to the 
normal distribution parameters defined by \code{g}, e.g. it can be normal
distribution with 0-mean and std dev 1. }
\item{z}{transformed data, i.e. either p-value or enrichment score}
}
\value{
\item{res}{statitics of Kolmogorov-Smirnov Test result obtained for 
\code{z} values}
}
\examples{
set.seed(1)
## let us assume we have a set of normalized scores: 
z <- abs(rnorm(10)) ## it should be positive and at least 10 length-vector
z <- z/median(z)
## Find the best log transform.
gamma <- optim(par=1.0, fn=tool.normalize.quality, gr=NULL, z,
lower=-9, upper=9, control=list(reltol=1e-3))
## After finding the best log transform, apply transform:
z <- log(exp(gamma$par)*z + 1.0)
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.normalize}}
}
