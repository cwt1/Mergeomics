\name{kda2cytoscape.colorize}
\alias{kda2cytoscape.colorize}
\title{
Trace module memberships of genes
}
\description{
\code{kda2cytoscape.colorize} assigns color to each node of the given 
module. If a node belongs to more than one module, different colors will 
be assigned to that node, as each color representing one module (shared 
nodes are illustrated as pie charts in the graph).
}
\usage{
kda2cytoscape.colorize(noddata, moddata, modpool, palette)
}
\arguments{
\item{noddata}{
node information of the entire graph
}
\item{moddata}{
module data including node (member gene) list
}
\item{modpool}{
unique module list including significant key drivers 
}
\item{palette}{
assigned unique color map for all modules 
}
}
\value{
\item{res }{data frame including the assigned color labels for the nodes of
the given module. If a node is concurrently member of many modules, multiple
colors will be assigned to that node (one color for each of these modules)
}
}
\references{
Shu L, Zhao Y, Kurt Z, Byars S, Tukiainen T, Kettunen J, Ripatti S, Zhang B,
Inouye M, Makinen VP, Yang X. Mergeomics: integration of diverse genomics
resources to identify pathogenic perturbations to biological systems.
bioRxiv doi: http://dx.doi.org/10.1101/036012
}
\examples{
## Trace module memberships for each KD and its neighbors
## If a KD (and its neigbor nodes) is member of multiple modules, assign 
## multiple colors to these multi-member nodes.
## We need to know data of all possible modules and all possible module ids
## to assign multiple colors to a shared node (between modules) when needed
if(exists("valdata"))
cat("Marker pvalues will be used to determine node sizes 
in the network illustration")
# noddata <- kda2cytoscape.colorize(neighs, job.kda$moddata, modpool, palette)
}
\author{
Zeyneb Kurt
}
\seealso{
\code{\link{kda2cytoscape}}
}

