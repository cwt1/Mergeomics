\name{tool.fdr.bh}
\alias{tool.fdr.bh}
\title{
Benjamini and Hochberg False Discovery Rate
}
\description{
\code{tool.fdr.bh} estimates FDRs of modules by using Benjamini and 
Hochberg method.  
}
\usage{
tool.fdr.bh(p)
}
\arguments{
\item{p}{p-values of modules}
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
FDRs <- tool.fdr.bh(p) ## the default method is already Benjamini Hochberg
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.fdr}}, \code{\link{tool.fdr.bh}}
}
