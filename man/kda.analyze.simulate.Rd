\name{kda.analyze.simulate}
\alias{kda.analyze.simulate}
\title{
Weighted key driver analysis (wKDA) simulation 
}
\description{
Generates simulations for permutation test, which is performed to obtain 
the p-value for the enrichment score of a given hub for a specified module
during the wKDA process.
}
\usage{
kda.analyze.simulate(o, g, nmemb, nnodes, nsim)
}
\arguments{ 
\item{o}{
Observed enrichment score of a hub node assigned for a given module.
}
\item{g}{
Sub-graph of a given hub and its neighbors (hubnet).
}
\item{nmemb}{
Number of the members included in a given module.
}
\item{nnodes}{
Number of the nodes in the whole graph (network) of the dataset.
}
\item{nsim}{
Number of the iterations (simulations) performed for the permutation 
test.
}
}
\details{
\code{\link{kda.analyze.simulate}} performs permutation tests to obtain 
p-values for the enrichment score of a given hub node for a given module.
It takes the observed enrichment score of the given hub, hubnet (subgraph 
of the hub and its neighbors), number of the members of the given module, 
total number of the nodes in the entire graph of the dataset, and number of
the simulations for the permutation test. 
In each iteration (simulation), it samples \code{nmemb} nodes randomly 
among the entire nodes of the graph. Then, it tests the overlapped nodes 
among the randomly chosen nodes and the given node's neigborhood. At the 
end, it obtains an enrichment score for each simulation and evaluates these
permuted enrichment scores with respect to the observed enrichment score of
the hub. Among \code{nsim} random simulations; maximally, enrichment scores
of 10 iterations are allowed to be greater than the observed (actual)
enrichment score of the hub. If this limitation is exceeded, simulation 
will be finalized at that point and the enrichment score list of the 
iterations will be returned. 
}
\value{
\item{x}{A list containing enrichment scores of the simulation's 
iterations}
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
    
## Observed enrichment scores.
# obs <- rep(NA, nhubs)
# k <- 1 ## actual using: for(k in 1:nhubs){}, for unit test, use the 1st hub
# g <- hubnets[[hubs[k]]]
# obs[k] <- kda.analyze.test(g$RANK, g$STRENG, memb, nnodes)

## Estimate P-values.
# pvals <- rep(NA, nhubs)
# for(k in which(obs > 0)) {
# g <- hubnets[[hubs[k]]]
## First pass:
# x <- kda.analyze.simulate(obs[k], g, nmemb, nnodes, 200)
## Then, use x to estimate preliminary and final P-values. 
## See kda.analyze() for more detail

## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
## remove the results folder
unlink("Results", recursive = TRUE)
# } ## finishing for loop
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
\code{\link{kda.analyze}}, \code{\link{kda.analyze.exec}},
\code{\link{kda.analyze.test}}
}
