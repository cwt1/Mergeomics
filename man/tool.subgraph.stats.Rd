\name{tool.subgraph.stats}
\alias{tool.subgraph.stats}
\title{
Calculate node degrees and strengths
}
\description{
\code{tool.subgraph.stats} graph statistics (degrees and strengths) of the 
seed nodes obtained from their neighborhoods.
}
\usage{
tool.subgraph.stats(frame, edgemap, heads, weights)
}
\arguments{
\item{frame}{a data frame including following components: \preformatted{
RANK: indices of neighboring nodes (including seeds)
LEVEL: number of edges away from seed 
STRENG: sum of adjacent edge weights within neighborhood
DEGREE: number of adjacent edges within neighborhood
}
}
\item{edgemap}{list of adjacent edge information for detected neighborhoods 
of seed nodes. \code{edgemap} can belong to either tails or heads.}
\item{heads}{list of either head (destination) or tail (source) nodes for 
neighborhoods of the seed nodes}
\item{weights}{weights of the edges in the entire graph}
}
\value{a data list including seed nodes neighborhood information with 
following components:
\item{RANK}{indices of neighboring nodes (including seeds)}
\item{LEVEL}{number of edges away from seed }
\item{STRENG}{sum of adjacent edge weights within neighborhood}
\item{DEGREE}{number of adjacent edges within neighborhood}
}
\examples{
data(job_kda_analyze)
depth <- 1
direction <- 0
## Take one or multiple center nodes (seeds) to search the neighborhoods:
## e.g. take the first node in the graph as the seed, find its neighborhood:
center.node = job.kda$graph$nodes[1]
## Convert center node (seed) names to indices:
nodes <- job.kda$graph$nodes
ranks <- match(center.node, nodes)
ranks <- ranks[which(ranks > 0)]
## we already know that rank is 1, since we took the first node in the graph
## as an example:
ranks <- as.integer(ranks) 
## Find edges to adjacent nodes. (both up- and down-stream searches)
visited <- ranks
levels <- 0*ranks
for(i in 1:depth) {
## Find edges to adjacent nodes.
foundT <- tool.subgraph.find(ranks, job.kda$graph$tail2edge, 
job.kda$graph$heads, visited)
foundH <- tool.subgraph.find(ranks, job.kda$graph$head2edge, 
job.kda$graph$tails, visited)        
## Expand neighborhood for the further depths of the neighborhood search
ranks <- unique(c(foundT, foundH))
visited <- c(visited, ranks)
levels <- c(levels, (0*ranks + i)) ## level shows the depth
if(length(ranks) < 1) break
}
## Calculate node degrees and strengths.
res <- data.frame(RANK=visited, LEVEL=levels, DEGREE=0,
STRENG=0.0, stringsAsFactors=FALSE)
res <- tool.subgraph.stats(res, job.kda$graph$tail2edge, 
job.kda$graph$heads, job.kda$graph$weights)
res <- tool.subgraph.stats(res, job.kda$graph$head2edge,
job.kda$graph$tails, job.kda$graph$weights)
}
\author{
Ville-Petteri Makinen 
}
