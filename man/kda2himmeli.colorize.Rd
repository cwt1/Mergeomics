\name{kda2himmeli.colorize}
\alias{kda2himmeli.colorize}
\title{
Trace module memberships of genes
}
\description{
\code{kda2himmeli.colorize} assigns color to each node of the given 
module. If a node belongs to more than one module, different colors will 
be assigned to that node, as each color representing one module (shared 
nodes are illustrated as pie charts in the graph).
}
\usage{
kda2himmeli.colorize(noddata, moddata, modpool, palette)
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
\item{res }{data frame including the assigned color labels for the nodes 
of the given module. If a node is concurrently member of many modules, 
many colors will be assigned to that node (one color for each of these 
modules)
}
}
\references{
Shu L, Zhao Y, Kurt Z, Byars SG, Tukiainen T, Kettunen J, Orozco LD, 
Pellegrini M, Lusis AJ, Ripatti S, Zhang B, Inouye M, Makinen V-P, Yang X.
Mergeomics: multidimensional data integration to identify pathogenic 
perturbations to biological systems. BMC genomics. 2016;17(1):874.
}
\examples{
## Trace module memberships for each KD
## If a KD is member of multiple modules, assign multiple colors to it
## Also consider the locus pval of the top locus of each KD (by valdata)
## We need to know data of all possible modules and all possible module ids
## to assign multiple colors(sectors) to a KD when needed
if(exists("valdata"))
cat("Marker pvalues will be used to determine node sizes 
in the network illustration")
# noddata <- kda2himmeli.colorize(valdata, job.kda$moddata, modpool, palette)
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{kda2himmeli}}
}

