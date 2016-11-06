\name{kda.analyze.exec}
\alias{kda.analyze.exec}
\title{
Auxiliary function for weight key driver analysis (wKDA) 
}
\description{
Obtains the enrichment scores (and p-values of these scores) of the hub 
nodes by the module member genes for a given module.
The hub node enrichment P-values reflect the degree of enrichment of hub's 
neighbor nodes within the member genes of the module, to whom this hub 
belongs to, when compared to the null distribution of randomly expected 
enrichment of hub within graph's nodes.
}
\usage{
kda.analyze.exec(memb, graph, nsim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{memb}{Member nodes of the given module.}
\item{graph}{Entire graph (network) of the dataset.}
\item{nsim}{Number of the simulations for the permutation test to obtain 
p-values of the enrichment scores belonging to the hub nodes for a given
module.}
}
\details{
\code{\link{kda.analyze.exec}} obtains the p-values of the enrichment scores 
belonging to the hub nodes for a given module. Enrichment score of a hub node 
for a given module is obtained by the overlapped (shared) nodes between this 
hub's neighbor nodes and the member nodes of the given module. If a hub node 
does not have at least a particular number of neighbors, its enrichment score 
is assigned as 0.0. 
}
\value{
\item{pvals }{
P-values of the enrichment scores belonging to the hub nodes for the given
module.}
}
\examples{
## This auxiliary function is called by kda.analyze(), 
## see this main function for more details
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
## calculate p-vals of KDs for the specified module:
# p <- kda.analyze.exec(memb, graph, nsim) ## see kda.analyze() for details

## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
## remove the results folder
unlink("Results", recursive = TRUE)
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
\code{\link{kda.analyze}}, \code{\link{kda.analyze.simulate}},
\code{\link{kda.analyze.test}}
}
