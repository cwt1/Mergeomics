\name{kda.prepare.screen}
\alias{kda.prepare.screen}
\title{
Prepare hubs and hubnets
}
\description{
\code{\link{kda.prepare.screen}} finds hubs and their neighborhoods
(hubnets) from the given graph.
}
\usage{
kda.prepare.screen(graph, depth, direction, efactor, dmin, dmax)
}
\arguments{
\item{graph}{
entire graph, whose hubs and hubnets will be obtained
}
\item{depth}{
search depth for subgraph search
}
\item{direction}{
the direction of the interactions among graph components. 0 for 
undirected, negative for downstream, and positive for upstream
}
\item{efactor}{
influence of node strengths (weights): 0.0 no influence, 1.0 full 
influence
}
\item{dmin}{
minimum hub degree to include
}
\item{dmax}{
maximum hub degree to include
}
}
\value{
\item{graph}{Updated graph including obtained hubs and hubnets:
\preformatted{
hubs: hub nodes list
hubnets: neighborhoods of hubs (hubnets)
}
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

## Then, extract overlapping co-hubs by kda.prepare.overlap()
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
\seealso{
\code{\link{kda.analyze}}, \code{\link{kda.prepare}}
}
