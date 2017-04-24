\name{tool.subgraph.find}
\alias{tool.subgraph.find}
\title{
Find edges to adjacent nodes
}
\description{
\code{tool.subgraph.find} finds the edge lists between given seed nodes 
and their neighbors
}
\usage{
tool.subgraph.find(seeds, edgemap, heads, visited)
}
\arguments{
\item{seeds}{seed nodes' indices }
\item{edgemap}{list of adjacent edge information for entire graph. 
\code{edgemap} can belong to either tails or heads.}
\item{heads}{list of either head (destination) or tail (source) nodes of 
the entire graph}
\item{visited}{flag holding already visited node indices during
neighborhood searching}
}
\value{
\item{neighbors }{neighbor edge lists of seed nodes (for either tails or 
heads)}
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
foundT <- tool.subgraph.find(ranks, job.kda$graph$tail2edge, 
job.kda$graph$heads, visited)
foundH <- tool.subgraph.find(ranks, job.kda$graph$head2edge, 
job.kda$graph$tails, visited)
}
\author{
Ville-Petteri Makinen 
}
\note{
Neighbor edge lists of the seed nodes should be obtained separately for
tail and head nodes.
}
