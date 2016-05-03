\name{kda.start.modules}
\alias{kda.start.modules}
\title{
Import module descriptions
}
\description{
\code{\link{kda.start.modules}} searches the whole nodes of the modules 
within the nodes of \code{edgdata} edgelist, filters out the nodes that 
does not exist in the nodes of \code{edgdata}, and deletes the modules,
which does not have enough nodes.
}
\usage{
kda.start.modules(job, edgdata)
}

\arguments{
\item{job}{a data frame including information such as module data file 
name, edge direction, minimum acceptable module size (\code{job$modfile}, 
\code{job$direction}, \code{job$minsize}, and so on.) }
\item{edgdata}{edge list data obtained from \code{\link{kda.start.edges}} }
}
\value{
\item{moddata }{module descriptions and their member node lists}
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

## Import topology of the graph for KDA, then find the module statistics
## This is already had been done in the kda.start() main function, due to
## the time constraint while running examples, we did not run it again.
# edgdata <- kda.start.edges(job.kda)
## Find module memberships of the graph nodes and obtain module statistics:
# moddata <- kda.start.modules(job.kda, edgdata)

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
\code{\link{kda.analyze}}, \code{\link{kda.finish}}, 
\code{\link{kda.prepare}}, \code{\link{kda.start}} 
}
