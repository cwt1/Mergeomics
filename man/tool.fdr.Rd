\name{tool.fdr}
\alias{tool.fdr}
\title{
Estimate False Discovery Rates (FDR)
}
\description{
\code{tool.fdr} estimates FDRs for modules as another module statistic.
}
\usage{
tool.fdr(p, f = NULL)
}
\arguments{
\item{p}{p-values of modules}
\item{f}{pre-defined threshold for FDR}
}
\details{
FDRs of modules can be obtained by using either empirical method or 
Benjamini and Hochberg method.
}
\value{
\item{res }{data list including the estimated false discovery rates of 
modules }
}
\examples{
## let us assume we have a set of pvalues
## and would like to find FDR values:
set.seed(1)
p <- abs(rnorm(10))*1e-2
FDRs <- tool.fdr(p) ## default method is Benjamini Hochberg
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.fdr.empirical}}, \code{\link{tool.fdr.bh}}
}
