\name{ssea.analyze.statistic}
\alias{ssea.analyze.statistic}
\title{
MSEA statistics for enrichment score
}
\description{
\code{\link{ssea.analyze.statistic}} estimates the enrichment score based
on observed and expected ones.
}
\usage{
ssea.analyze.statistic(o, e)
}
\arguments{
\item{o}{
observed enrichment score
}
\item{e}{
expected enrichment score
}
}
\value{
\item{score }{estimated enrichment score based on observed and expected 
scores}
}
\examples{
## O and E the observed and expected counts of positive findings 
## (enrichment scores) at a given cutoff:
set.seed(1)
o <- rnorm(1)
e <- rnorm(1)
## find the final enrichment score from the observed and estimated scores:
z <- ssea.analyze.statistic(o, e)
}
\references{
Shu L, Zhao Y, Kurt Z, Byars SG, Tukiainen T, Kettunen J, Orozco LD, 
Pellegrini M, Lusis AJ, Ripatti S, Zhang B, Inouye M, Makinen V-P, Yang X.
Mergeomics: multidimensional data integration to identify pathogenic 
perturbations to biological systems. BMC genomics. 2016;17(1):874.
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{ssea.analyze}}
}
