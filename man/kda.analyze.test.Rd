\name{kda.analyze.test}
\alias{kda.analyze.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate enrichment score for wKDA
}
\description{
Obtains the enrichment score of a given hub (center node) belonging to a 
specified module. Enrichment score of a center node depends on the shared 
node number between the neighbor nodes of this center node (derived from 
the provided graph topology) and member nodes of this center node's module.
The more a center node has neighbors in the graph among the member genes 
belonging to the module of this center node, the greater enrichment score 
it has.
}
\usage{
kda.analyze.test(neigh, w, members, nnodes)
}

\arguments{
\item{neigh}{
Neighbor nodes of the given hub node (i.e. nodes in the hubnet)
}
\item{w}{
Weigths of the given hub node based on its in-degree and out-degree 
edge density in the hubnet
}
\item{members}{
Node indices -within the entire graph- of the member genes of given 
hub's module.
}
\item{nnodes}{
Number of the nodes in the entire graph of the dataset.
}
}
\details{
\code{\link{kda.analyze.test}} takes a hub node's neigbor list and weight 
list; additionally, it takes the member node list of relevant module. It 
searches the masses of the shared nodes between hubnet and the given module 
(gene set). 
The shared edge mass is normalized with respect to number of the expected 
match ratio between hubnet and the given node list. This normalized ratio 
is assigned as the observed enrichment score of the hubnet according to 
the given member node list. 
}
\value{
\item{z}{Calculated enrichment score}
}
\examples{
job.kda <- list()
job.kda$label<-"HDLC"
## parent folder for results
job.kda$folder<- "Results"
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
job.kda$direction<-1
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

## Let's prepare KDA object for KDA:
job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
set.seed(job.kda$seed)
i = 1 ## index of the module, whose p-val is calculated:
memb <- job.kda$module2nodes[[i]]
graph <- job.kda$graph  ## we need to import a network
nsim <- job.kda$nperm   ## number of simulations
## This auxiliary function is called by kda.analyze.exec(), which is called
## by kda.analyze() main function, see this main function for more details

hubs <- graph$hubs
hubnets <- graph$hubnets
nhubs <- length(hubs)
nnodes <- length(graph$nodes)
nmemb <- length(memb)
    
## Observed enrichment scores for the hubs of the given module.
obs <- rep(NA, nhubs)
k <- 1 ## actual using: for(k in 1:nhubs){}, for test, use only the 1st hub
g <- hubnets[[hubs[k]]]
obs[k] <- kda.analyze.test(g$RANK, g$STRENG, memb, nnodes)

## Then, estimate preliminary and final P-values by kda.analyze.simulate()
## See kda.analyze() for more details

## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
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
\code{\link{kda.analyze}}, \code{\link{kda.analyze.exec}}, 
\code{\link{kda.analyze.simulate}}
}

