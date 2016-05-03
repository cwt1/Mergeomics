\name{tool.fdr.empirical}
\alias{tool.fdr.empirical}
\title{
Estimate Empirical False Discovery Rates
}
\description{
\code{tool.fdr.empirical} estimates empirical FDR for modules
}
\usage{
tool.fdr.empirical(p, f0)
}
\arguments{
\item{p}{p-values of modules}
\item{f0}{pre-defined threshold for FDR}
}
\value{
\item{res }{data list including the estimated false discovery rates of 
modules}
}
\examples{
## let us assume we have a set of pvalues
## and would like to find FDR values:
set.seed(1)
p <- abs(rnorm(10))*1e-2
f = 0.05 ## pre-defined threshold for FDR
FDRs <- tool.fdr.empirical(p, f) 
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.fdr}}, \code{\link{tool.fdr.bh}}
}
