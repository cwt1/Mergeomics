\name{tool.graph}
\alias{tool.graph}
\title{
Convert an edge list to a graph representation
}
\description{
\code{tool.graph} translates an edge list including TAIL, HEAD and WEIGHT 
information into a graph representation-adapted data list. It also 
provides in-degree and out-degree statistics for nodes.
}
\usage{
tool.graph(edges)
}
\arguments{
\item{edges}{a data frame with three columns TAIL, HEAD and WEIGHT}
}
\value{a datalist including following components:
\item{nodes }{N-element array of node names}
\item{tails }{K-element array of node indices}
\item{heads}{K-element array of node indices}
\item{weights}{K-element array of edge weights}
\item{tail2edge}{N-element list of adjacent edge indices}
\item{head2edge}{N-element list of adjacent edge indices}
\item{outstats}{N-row data frame of out-degree node statistics}
\item{instats}{N-row data frame of in-degree node statistics}
\item{stats}{N-row data frame of node statistics}
}
\examples{
job.kda <- list()
job.kda$label<-"HDLC"
## parent folder for results
job.kda$folder<-"Results"
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-system.file("extdata","network.mouseliver.mouse.txt", 
package="Mergeomics")
## module file:
job.kda$modfile<- system.file("extdata","mergedModules.txt", 
package="Mergeomics")
## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-0.0
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 1
job.kda$nperm <- 20 # the default value is 2000, use 20 for unit tests

## kda.start() process takes long time while seeking hubs in the given net
## Here, we used a very small subset of the module list (1st 10 mods
## from the original module file):
moddata <- tool.read(job.kda$modfile)
mod.names <- unique(moddata$MODULE)[1:min(length(unique(moddata$MODULE)),
10)]
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
## save this to a temporary file and set its path as new job.kda$modfile:
tool.save(moddata, "subsetof.supersets.txt")
job.kda$modfile <- "subsetof.supersets.txt"

job.kda <- kda.configure(job.kda)
## Import data for weighted key driver analysis:
## Import topology.
edgdata <- kda.start.edges(job.kda)
## Create an indexed graph structure.
job.kda$graph <- tool.graph(edgdata)

## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.subgraph}}
}
