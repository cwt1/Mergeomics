\name{kda2cytoscape}
\alias{kda2cytoscape}
\title{
Generate input files for Cytoscape
}
\description{
\code{kda2cytoscape} generates input files for Cytoscape to visualize the 
graph and hubnets after the wKDA process finished. The network 
visualization is a streamlined depiction of the module enrichment in hub 
neighborhoods.
}
\usage{
kda2cytoscape(job, node.list = NULL, modules = NULL, ndrivers = 5, 
depth = 1)
}
\arguments{
\item{job}{
wKDA result data list as returned by \code{\link{kda.finish}}
}
\item{node.list}{
array of node/gene names to be visualized with their neighbor node. 
if this is not specified top \code{ndrivers} of each module and their
neighborhoods will be illustrated.
}
\item{modules}{
array of module names to be visualized
}
\item{ndrivers}{
maximum number of drivers per module
}
\item{depth}{
depth for neighborhood search in the graph
}
}
\details{
\code{\link{kda2cytoscape}} first, selects top scoring key drivers for 
each module; then, assigns a colormap to modules, processes each module
separately, finds key nodes' neighborhoods within a particular search 
depth, and saves the edge and node lists of the modules to the specified 
output folder. Besides, it returns this configuration data to the user. 
Created file list for Cytoscape are given below: 
\preformatted{
kda2cytoscape.top.kds.txt: top key drivers of the modules are
listed in this file. Number of the key drivers can be set by
user with ndrivers parameter. 
kda2cytoscape.edges.txt: edge lists of the integrated graph 
that includes the subnetworks of all modules.
kda2cytoscape.nodes.txt: node lists of the integrated graph 
that includes the subnetworks of all modules.
module.color.mapping.txt: color mapping for the  modules, 
i.e. one color is assigned to each module.
}
}
\value{
\item{job }{updated data list including the node and edge information of
the modules converted to Cytoscape format}
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
## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-0.0
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 1

## finish the KDA process
job.kda <- kda.finish(job.kda)

## prepare the cytoscape-ready files:
job.kda <- kda2cytoscape(job.kda)

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
Zeyneb Kurt
}
\seealso{
\code{\link{kda.analyze}}, \code{\link{kda.finish}}
}
