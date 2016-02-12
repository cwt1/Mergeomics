\name{kda.prepare}
\alias{kda.prepare}
\title{
Prepare graph topology for weighted key driver analysis
}
\description{
\code{\link{kda.prepare}} gets graph topology required by wKDA process, 
then provides the information including hub list, hubnets, and overlapping
co-hubs. 
}
\usage{
kda.prepare(job)
}
\arguments{
\item{job}{
a parameter including restirictions while determining the graph topology
information (such as hubs, hubnets, co-hubs, etc.), which is required by 
the wKDA process: \preformatted{
graph: graph of the dataset
depth: search depth for subgraph search
direction: use 0 for undirected, negative for downstream 
and positive for upstream
maxoverlap: maximum allowed overlap between two key driver
neighborhoods
mindegree: minimum hub degree to include
edgefactor: influence of node strengths; 0.0 no influence,
1.0 full influence
}
}
}
\details{
\code{\link{kda.prepare}} determines minimum hub degree if it is not 
specified by the user, finds hubs and their neighborhoods (hubnets), 
extracts overlapping co-hubs, returns this information to user, and prints 
it to the screen.
}
\value{
\item{job }{Updated data frame including information about the graph 
topology in terms of hubs, hubnets, and overlapping co-hubs: 
\preformatted{
hubs: hub nodes list
hubnets: neighborhoods of hubs (hubnets)
cohubsets: overlapping hubs (co-hubs)
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

## Configure the parameters for KDA:
job.kda <- kda.configure(job.kda)
## Create the object properly
job.kda <- kda.start(job.kda)
## Find the hubs, co-hubs, and hub neighborhoods (hubnets), etc.:
job.kda <- kda.prepare(job.kda)
## After that, we need to call kda.analyze() and kda.finish()

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
\code{\link{kda.analyze}}, \code{\link{kda.prepare.overlap}}, 
\code{\link{kda.prepare.screen}}
}
