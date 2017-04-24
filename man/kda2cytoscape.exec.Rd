\name{kda2cytoscape.exec}
\alias{kda2cytoscape.exec}
\title{
Evaluate each module separately for visualization
}
\description{
\code{kda2cytoscape.exec} deals with the modules individually; takes a 
particular amount of top key drivers of the given module in company with
the top key driver lists and colormap of all modules; traces module 
memberships and produces colormap, it finds the edge and node lists for 
the top key drivers and their neighborhood for a given module. 
}
\usage{
kda2cytoscape.exec(job, drivers, modpool, palette, graph.depth = 1)
}
\arguments{
\item{job}{
data list including entire graph, nodes, modules information
}
\item{drivers}{
top key drivers of the specified module
}
\item{modpool}{
unique key driver list for all modules 
}
\item{palette}{
assigned unique color map for all modules 
}
\item{graph.depth}{
search depth for graph neighborhood 
}
}
\value{
\item{res }{uniquely identified node and edge lists of the members
belonging to the given module}
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
drivers <- as.data.frame(drivers)
colnames(drivers) <- c("MODULE" , "NODE")
    
mods <- unique(drivers$MODULE)
modnames <- job.kda$modules[mods] 
modnames[which(mods == 0)] <- "NON.MODULE"
palette <- kda2cytoscape.colormap(length(mods)) 
palette[,which(mods == 0)] <- c(90,90,90)
drivers$MODNAMES <- modnames[match(drivers$MODULE, mods)]
drivers$NODNAMES <- job.kda$graph$nodes[drivers$NODE] 
for(i in 1:nrow(drivers))
drivers$COLOR[i] <- paste(palette[1, match(drivers$MODULE[i], mods)],
palette[2, match(drivers$MODULE[i], mods)],
palette[3, match(drivers$MODULE[i], mods)], sep=" ")
## Process each module separately. Just perform for the 1st module:
i <- 1
rows <- which(drivers$MODULE == mods[i])
if(length(rows) > 0)
tmp <- kda2cytoscape.exec(job.kda, drivers[rows,], mods, palette,
job.kda$depth)

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

