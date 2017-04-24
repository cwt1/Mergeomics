\name{tool.subgraph.search}
\alias{tool.subgraph.search}
\title{
Search neighborhoods for given nodes
}
\description{
\code{tool.subgraph.search} looks for both upstream and downstream 
neighborhoods of given seed node list for a given depth, gets the directed 
edge information among seed nodes and their neighbors, obtains statistics 
(degrees and strengths) for seed nodes.
}
\usage{
tool.subgraph.search(graph, seeds, depth, direction)
}
\arguments{
\item{graph}{a datalist including following components: \preformatted{
nodes: N-element array of node names
tails: K-element array of node indices
heads: K-element array of node indices
weights:  K-element array of edge weights
tail2edge: N-element list of adjacent edge indices
head2edge: N-element list of adjacent edge indices
outstats: N-row data frame of out-degree node statistics
instats: N-row data frame of in-degree node statistics
stats: N-row data frame of node statistics
}
}
\item{seeds}{seed nodes' indices }
\item{depth}{the maximum number of links to connect neighbors}
\item{direction}{sets the directionality: use a negative value for dowstream,
positive for upstream or zero for undirected}
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
## Find neighbors.
res <- tool.subgraph.search(job.kda$graph, ranks, depth, direction)
}
\author{
Ville-Petteri Makinen 
}
