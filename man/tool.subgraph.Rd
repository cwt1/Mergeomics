\name{tool.subgraph}
\alias{tool.subgraph}
\title{
Determine network neighbors for a set of nodes
}
\description{
\code{tool.subgraph} finds the sub-network, i.e. neighborhood, for a given 
seed node list with a specified depth. It also provides graph statistics
(degrees and strengths) for seed nodes.
}
\usage{
tool.subgraph(graph, seeds, depth = 1, direction = 0)
}
\arguments{
\item{graph}{a datalist including following components: \preformatted{
nodes: N-element array of node names
tails: K-element array of node indices
heads: K-element array of node indices
weights: K-element array of edge weights
tail2edge: N-element list of adjacent edge indices
head2edge: N-element list of adjacent edge indices
outstats: N-row data frame of out-degree node statistics
instats: N-row data frame of in-degree node statistics
stats: N-row data frame of node statistics
}
}

\item{seeds}{list of seed node names }
\item{depth}{the maximum number of links to connect neighbors}
\item{direction}{sets the directionality: use a negative value for 
dowstream, positive for upstream or zero for undirected}
}
\value{a data list including following components:
\item{RANK}{indices of neighboring nodes (including seeds)}
\item{LEVEL}{number of edges away from seed }
\item{STRENG}{sum of adjacent edge weights within neighborhood}
\item{DEGREE}{number of adjacent edges within neighborhood}
}
\examples{
data(job_kda_analyze)
## take the first node in the graph as the seed, find its neighborhood:
center.node = job.kda$graph$nodes[1]
subnet = tool.subgraph(job.kda$graph, center.node, depth=1, direction=0)
}

\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{tool.subgraph.search}}
}
