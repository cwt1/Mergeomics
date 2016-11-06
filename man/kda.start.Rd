\name{kda.start}
\alias{kda.start}
\title{
Import data for weighted key driver analysis
}
\description{
\code{\link{kda.start}} converts identities (such as module descriptions, 
module identifiers, and module nodes) to indices. It prepares graph 
topology and module information for wKDA process.
}
\usage{
kda.start(job)
}
\arguments{
\item{job}{
a data frame including fields for edges and nodes information of the 
graph (TAIL, HEAD, WEIGHT). It also involves path of input files 
including module descriptions and module-gene lists.
}
}
\details{
\code{\link{kda.start}} imports graph and relevant module descriptor input
files, creates an indexed graph structure, and converts identities to 
indices from module descriptions and module-gene lists. Hence, it concludes
with a graph structure and a module set involving member gene IDs for
each module.
}
\value{
\item{job }{Updated data frame including indexed graph topology, modules,
and nodes information: \preformatted{
graph:     indexed topology
modules:   module identities
modinfo:   module descriptions (indexed)
moddata:   module data (indexed)
module2nodes:  lists of node indices for each module
modulesizes:   module sizes
}
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

job.kda <- kda.configure(job.kda)
## Import data for weighted key driver analysis:
job.kda <- kda.start(job.kda) 

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
\code{\link{kda.analyze}}, \code{\link{kda.finish}}, 
\code{\link{kda.prepare}}, \code{\link{kda.start.edges}}, 
\code{\link{kda.start.identify}}, \code{\link{kda.start.modules}}
}

