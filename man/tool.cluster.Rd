\name{tool.cluster}
\alias{tool.cluster}
\title{
Hierarchical clustering of nodes 
}
\description{
\code{tool.cluster} performs agglomerative hierarchical clustering for 
nodes (genes) 
}
\usage{
tool.cluster(edges, cutoff = NULL)
}
\arguments{
\item{edges}{
edge (weight) list among two group, whose overlapping information 
(overlapping ratio based on shared entries of two groups, number of 
members in both group) had been assesed previously 
}
\item{cutoff}{cutting level of dendrogram for hierarchical clustering}
}
\details{
\code{tool.cluster} takes overlapping information between two groups, 
produces distance matrix based on 1-strength(overlap) ratio between 
two groups, and apply agglomerative hierarchical clustering based on the 
distance matrix.
}
\value{
\item{res }{data list including clustering results: \preformatted{
CLUSTER: cluster label
NODE: item (node) name
}
}
}
\examples{
## read the coexpr module file as an example:
moddata <- tool.read(system.file("extdata",
"modules.mousecoexpr.liver.human.txt", package="Mergeomics"))

## let us cluster the first 10 modules in the module file:
mod.names <- unique(moddata$MODULE)[1:min(length(unique(moddata$MODULE)),
10)]
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
## Find clusters.
rmax = 0.33
edges <- tool.overlap(items=moddata$GENE, groups=moddata$MODULE)
clustdat <- tool.cluster(edges, cutoff=rmax)
nclust <- length(unique(clustdat$CLUSTER))
nnodes <- length(unique(clustdat$NODE))
}

\references{
Shu L, Zhao Y, Kurt Z, Byars S, Tukiainen T, Kettunen J, Ripatti S, Zhang B,
Inouye M, Makinen VP, Yang X. Mergeomics: integration of diverse genomics
resources to identify pathogenic perturbations to biological systems.
bioRxiv doi: http://dx.doi.org/10.1101/036012
}
\author{
Ville-Petteri Makinen 
}
