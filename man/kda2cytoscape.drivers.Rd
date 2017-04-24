\name{kda2cytoscape.drivers}
\alias{kda2cytoscape.drivers}
\title{
Select top key drivers for each module
}
\description{
\code{kda2cytoscape.drivers} finds maximally top \code{ndriv} key drivers 
for each module with respect to the significance level of the drivers. 
}
\usage{
kda2cytoscape.drivers(data, modules, ndriv)
}
\arguments{
\item{data}{
data frame including information of the modules (key driver list, 
p-values, node list, false discovery rates (fdr), and so on.)
}
\item{modules}{
top scoring modules among KDA results
}
\item{ndriv}{
maximum number of drivers that can be chosen for per module
}
}
\value{
\item{data }{top key drivers (maximally \code{ndriv} drivers for each
module) for top modules}
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

## Finish the KDA process
job.kda <- kda.finish(job.kda)
## Select top key drivers from each module.
## First, take module names from kda results
modules <- unique(job.kda$results$MODULE)
## Take top 2 KDs:
drivers <- kda2cytoscape.drivers(job.kda$results, modules, ndriv=2)

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
Zeyneb Kurt
}
\seealso{
\code{\link{kda2cytoscape}}
}

