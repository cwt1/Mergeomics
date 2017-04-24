\name{tool.aggregate}
\alias{tool.aggregate}
\title{
Aggregate the entries 
}
\description{
\code{tool.aggregate} aggregates the entries with respect to the given 
feature. It first finds raw indices (either genes or markers), then sorts
them, and finds the blocks (segments) of identical entries.
}
\usage{
tool.aggregate(entries, limit = 1)
}
\arguments{
\item{entries}{
an array that will be sorted and aggregated within blocks
}
\item{limit}{
minimum block size to be included
}
}

\value{
\item{res }{a data list with the following components: \preformatted{
labels: shared values within blocks
lengths: numbers of entries in blocks
blocks: integer arrays of entry positions within blocks
ranks: entry positions included in blocks
}
}
}
\examples{
aa <- data.frame(MODULE=c("Mod1", "Mod1", "Mod2", "Mod2", "Mod3"), 
GENE=c("GeneA", "GeneC", "GeneB", "GeneC", "GeneA"))
## aggregate according to the module names:
bb <- tool.aggregate(aa$MODULE)
bb
## aggregate according to the gene names:
cc <- tool.aggregate(aa$GENE)
cc
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