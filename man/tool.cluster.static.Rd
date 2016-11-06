\name{tool.cluster.static}
\alias{tool.cluster.static}
\title{
Static hierarchical clustering 
}
\description{
\code{tool.cluster.static} takes dendrogram (clustering tree) and its 
cutting height; then, obtains cluster labels for the nodes of the tree.
}
\usage{
tool.cluster.static(dendro, hlim)
}
\arguments{
\item{dendro}{dendrogram (tree) }
\item{hlim}{cutting height of the dendrogram}
}
\value{
\item{clusters }{cluster labels of the components after static clustering}
}
\examples{
set.seed(1)
## assume that we have a dataset including several samples 
## with distinct features
dataset <- matrix(rnorm(20), ncol=5) ## 4 samples with 5 features
## Find the distances between each sample pair to cluster them
d <- dist(dataset, method = "euclidean", upper=TRUE, diag=TRUE)
tree <- hclust(d)
## Height cutoff.
hlim <- max(tree$height)
## Find clusters.
clusters <- tool.cluster.static(tree, hlim)
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

