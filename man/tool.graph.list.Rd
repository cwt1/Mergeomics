\name{tool.graph.list}
\alias{tool.graph.list}
\title{
Return edge list for each node
}
\description{
\code{tool.graph.list} finds and returns the edge list of each node for 
both tail and head node lists.
}
\usage{
tool.graph.list(entries, nnodes)
}
\arguments{
\item{entries}{either tail nodes list or head nodes list}
\item{nnodes}{total number of all nodes including both tails and heads}
}
\value{
\item{groups }{a data list including edge list of each node}
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
edges <- kda.start.edges(job.kda)
## Create an indexed graph structure.
tails <- as.character(edges$TAIL)
heads <- as.character(edges$HEAD)
wdata <- as.double(edges$WEIGHT)

nedges <- length(tails)
# Create factorized representation.
labels <- as.character(c(tails, heads))
labels <- as.factor(labels)
labelsT <- as.integer(labels[1:nedges])
labelsH <- as.integer(labels[(nedges+1):(2*nedges)])
# Create edge lists.
nodnames <- levels(labels)
nnodes <- length(nodnames)
elistT <- tool.graph.list(labelsT, nnodes)
elistH <- tool.graph.list(labelsH, nnodes)
## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.graph}}
}
