\name{kda.start.edges}
\alias{kda.start.edges}
\title{
Import nodes and edges of graph topology
}
\description{
\code{\link{kda.start.edges}} imports network file, gets edge data 
(in TAIL, HEAD, WEIGHT format), eliminates the nodes -whose degree is 
smaller than the maximum allowed node degree-, and returns the edges of 
remaining nodes. 
}
\usage{
kda.start.edges(job)
}
\arguments{
\item{job}{a data frame including information such as network file name, 
maximum allowed node degree, edge direction (\code{job$netfile}, 
\code{job$maxdegree}, 
\code{job$direction}, and so on.) }
}
\value{
\item{edgdata }{filtered edge list, i.e. edges of the nodes, whose degree 
is smaller than the maximum allowed node degree}
}
\references{
Shu L, Zhao Y, Kurt Z, Byars S, Tukiainen T, Kettunen J, Ripatti S, Zhang B,
Inouye M, Makinen VP, Yang X. Mergeomics: integration of diverse genomics
resources to identify pathogenic perturbations to biological systems.
bioRxiv doi: http://dx.doi.org/10.1101/036012
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
job.kda <- kda.configure(job.kda)

## Import topology of the graph for KDA
## This is already had been done in the kda.start() main function, due to
## the time constraint while running examples, we did not run it again.
# edgdata <- kda.start.edges(job.kda)

## remove the results folder
unlink("Results", recursive = TRUE)
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{kda.analyze}}, \code{\link{kda.finish}}, 
\code{\link{kda.prepare}}, \code{\link{kda.start}} 
}
