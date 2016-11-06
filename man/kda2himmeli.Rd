\name{kda2himmeli}
\alias{kda2himmeli}
\title{
Generate input files for Himmeli
}
\description{
\code{kda2himmeli} generates input files for Himmeli to visualize the 
graph and hubnets after the wKDA process finished. The network visualization
is a streamlined depiction of the module enrichment in hub neighborhoods.
}
\usage{
kda2himmeli(job, modules = NULL, ndrivers = 5)
}
\arguments{
\item{job}{
KDA result data list as returned by \code{\link{kda.finish}}
}
\item{modules}{
array of module names to be visualized
}
\item{ndrivers}{
maximum number of drivers per module
}
}
\details{
\code{\link{kda2himmeli}} first, selects top scoring key drivers for each 
module; then, assigns a colormap to modules, processes each module 
separately, finds key nodes' neighborhoods, and saves the edge and node 
lists of the modules to the specified output folder. Besides, it returns 
this configuration data to the user. 
}
\value{
\item{job }{updated data list including the node and edge information of 
the modules converted to Himmeli format}
}
\examples{
## get the prepared and KDA applied dataset:(see kda.analyze for details)
data(job_kda_analyze)
## set the relevant parameters:
job.kda$label<-"HDLC"
## parent folder for results
job.kda$folder<-"Results"
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-system.file("extdata","network.mouseliver.mouse.txt", 
package="Mergeomics")
job.kda$nodfile <- system.file("extdata","msea2kda.nodes.txt", 
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

## finish the KDA process
job.kda <- kda.finish(job.kda)

## prepare the cytoscape-ready files:
job.kda <- kda2himmeli(job.kda)

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
\code{\link{kda.analyze}}, \code{\link{kda.finish}}
}

