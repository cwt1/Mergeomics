\name{kda2himmeli.edges}
\alias{kda2himmeli.edges}
\title{
Find edges of a given node with a specified depth
}
\description{
\code{kda2himmeli.edges} finds the sub-graph (node and edge lists) of a
central node and its neighborhood at a particular search depth. The 
central node is a member of a module, which is defined at 
\code{\link{kda2himmeli.exec}}. 
}
\usage{
kda2himmeli.edges(graph, center, depth, direction)
}
\arguments{
\item{graph}{
entire graph
}
\item{center}{
the node, whose interactions with neighbors will be searched within
\code{graph}.
}
\item{depth}{
search depth for graph neighborhood 
}
\item{direction}{
edge direction. 0 for undirected, negative for downstream and 
positive for upstream
}
}
\value{
\item{g }{the sub-graph including TAIL, HEAD, WEIGHT information of the
central node, which belongs to the specified module.}
}
\examples{
## get the prepared and KDA applied dataset:(see kda.analyze for details)
data(job_kda_analyze)
## set the relevant parameters:
job.kda$label<-"HDLC"
## parent folder for results
job.kda$folder<-"Results"
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-system.file("extdata","network.mouseliver.mouse.txt", 
package="Mergeomics")
## Gene sets derived from ModuleMerge, containing two columns, MODULE, 
## NODE, delimited by tab 
job.kda$modfile<- system.file("extdata","mergedModules.txt", 
package="Mergeomics")
job.kda$nodfile <- system.file("extdata","msea2kda.nodes.txt", 
package="Mergeomics")
## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-0.0
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 1

## Finish the KDA process
job.kda <- kda.finish(job.kda)
## Select a center node to seek its neighbors in the graph:
edges.of.center.node <- kda2himmeli.edges(job.kda$graph, 1, 
job.kda$depth, job.kda$direction)

## remove the results folder
unlink("Results", recursive = TRUE)
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
\seealso{
\code{\link{kda2himmeli}}
}

