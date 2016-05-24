\name{kda.analyze}
\alias{kda.analyze}
\title{
Weighted key driver analysis (wKDA) main function
}
\description{
Finds the statistics (enrichment score, p-value, FDR, etc.) of the key 
driver (hub) genes belonging to the specified modules based on the graph
topology. The enrichment score of a hub node based on the shared nodes 
between this hub's neighbor nodes in the graph and the member nodes of the 
hub's module. The hub node enrichment P-values reflect the degree of 
observed enrichment of the hub, when compared to the null distribution of
randomly expected enrichment of this hub within graph's nodes. Permutation
test is used to obtain these statistics. 
}
\usage{
kda.analyze(job)
}
\arguments{
\item{job}{
The data list that will be subjected to KDA. It involves the modules, 
member genes belonging to each module, graph (network) topology, hubs 
of the graph, and sub-graph around each hub (hubnets of the graph).
}
}
\details{
\code{\link{kda.analyze}} analyzes each module individually and determines
the p-values and FDRs of hub nodes of each module via permutation test. 
It returns the hit hub (key driver) gene name and member list of each module.
}
\value{
\item{job}{
The KDA applied data list. It involves the modules, hub gene and
member genes belonging to each module, and False Discovery Rate (FDR) 
adjusted p-values of hub nodes for each module.
}
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

## Let's run KDA!
\donttest{
job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)
}

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
\code{\link{kda.analyze.exec}}, \code{\link{kda.analyze.simulate}}, 
\code{\link{kda.analyze.test}}
}
