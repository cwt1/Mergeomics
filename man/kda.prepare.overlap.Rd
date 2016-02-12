\name{kda.prepare.overlap}
\alias{kda.prepare.overlap}
\title{
Extract overlapping co-hubs
}
\description{
\code{\link{kda.prepare.overlap}} finds overlapping co-hubs of the given 
graph.
}
\usage{
kda.prepare.overlap(graph, direction, rmax)
}
\arguments{
\item{graph}{
entire graph, whose overlapping co-hubs will be found
}
\item{direction}{
the direction of the interactions among graph components. 0 for 
undirected, negative for downstream, and positive for upstream
}
\item{rmax}{
maximum allowed overlap between two key driver neighborhoods
}
}
\value{
\item{graph}{
Updated graph including overlapping co-hubs: \preformatted{
cohubsetsco-hubs of the given graph}
}
}
\examples{
job.kda <- list()
job.kda$label<-"HDLC"
## parent folder for results
job.kda$folder<- "Results"
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-"network.mouseliver.mouse.txt"
## Gene sets derived from ModuleMerge, containing two columns, MODULE, 
## NODE, delimited by tab 
job.kda$modfile<- "mergedModules.txt"
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 1

## Configure the parameters for KDA:
# job.kda <- kda.configure(job.kda)
## Create the object properly
# job.kda <- kda.start(job.kda)

## Find the hubs, co-hubs, and hub neighborhoods (hubnets) by kda.prepare()
## and its auxiliary functions kda.prepare.screen and kda.prepare.overlap
## First, determine the minimum and maximum hub degrees:
# nnodes <- length(job.kda$graph$nodes)
# if (job.kda$mindegree == "automatic") {
#   dmin <- as.numeric(quantile(job.kda$graph$stats$DEGREE,0.75))
#   job.kda$mindegree <- dmin
# }
# if (job.kda$maxdegree == "automatic") {
#   dmax <- as.numeric(quantile(job.kda$graph$stats$DEGREE,1))
#   job.kda$maxdegree <- dmax
# } 
## Collect neighbors.
# job.kda$graph <- kda.prepare.screen(job.kda$graph, job.kda$depth, 
# job.kda$direction, job.kda$edgefactor, job.kda$mindegree, job.kda$maxdegree)

## Then, extract overlapping co-hubs by kda.prepare.overlap():
## Collect overlapping co-hubs.
# job.kda$graph <- kda.prepare.overlap(job.kda$graph, job.kda$direction,
# job.kda$maxoverlap)

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
\code{\link{kda.analyze}}, \code{\link{kda.prepare}}
}
