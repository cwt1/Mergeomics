\name{tool.unify}
\alias{tool.unify}
\title{
Convert a distribution to uniform ranks 
}
\description{
\code{tool.unify} converts a distribution to uniform ranks with respect 
to a background distribution (or self if no background available).}
\usage{
tool.unify(xtrait, xnull = NULL)
}
\arguments{
\item{xtrait}{the distribution that will be standardized, i.e. uniformly
distributed}
\item{xnull}{background distribution to be used to distribute 
\code{xtrait} uniformly. If \code{xnull} is not specified, \code{xtrait}
will be used as background distr}
}
\value{
\item{y }{uniformly distributed form of \code{xtrait}}
}
\examples{
x <- rnorm(10)
y <- tool.unify(x) ## uniformly distributed form of x when null dist is x
z <- tool.unify(x, y) ## uniformly distributed form of x when null dist is y
}
\author{
Ville-Petteri Makinen 
}
