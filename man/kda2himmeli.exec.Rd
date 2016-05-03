\name{kda2himmeli.exec}
\alias{kda2himmeli.exec}
\title{
Evaluate each module separately for visualization
}
\description{
\code{kda2himmeli.exec} deals with the modules individually; takes a 
particular amount of top key drivers of the given module in company with 
the top key driver lists and colormap of all modules; traces module
memberships and produces colormap, it finds the edge and node lists for 
the top key drivers and their neighborhood for a given module. 
}
\usage{
kda2himmeli.exec(job, valdata, drivers, modpool, palette)
}
\arguments{
\item{job}{
data list including entire graph, nodes, modules information
}
\item{valdata}{
GWAS pvalues of top loci of the nodes - if this information is available,
sizes of the nodes in the figure will be correlated with the p-value of 
the top loci of the nodes - 
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
job.kda$nodfile <- system.file("extdata","msea2kda.nodes.txt", 
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

## Get valdata including marker pvals
valdata <- tool.read(job.kda$nodfile)
z <- as.double(valdata$VALUE)
z <- (z/quantile(z, 0.95) + rank(z)/length(z))
valdata$SIZE <- pmin(4.0, z)
## Select subset of genes.
valdata <- kda2himmeli.identify(valdata, "NODE", job.kda$graph$nodes)

## Select top key drivers from each module.
## First, take module names from kda results
modules <- unique(job.kda$results$MODULE)
## Take top 2 KDs:
drivers <- kda2himmeli.drivers(job.kda$results, modules, ndriv=2)
drivers <- as.data.frame(drivers)
colnames(drivers) <- c("MODULE" , "NODE")

mods <- unique(drivers$MODULE)
modnames <- job.kda$modules[mods] 
modnames[which(mods == 0)] <- "NON.MODULE"
palette <- kda2himmeli.colormap(length(mods)) 
palette[,which(mods == 0)] <- c(90,90,90)
drivers$MODNAMES <- modnames[match(drivers$MODULE, mods)]
drivers$NODNAMES <- job.kda$graph$nodes[drivers$NODE] 
for(i in 1:nrow(drivers))
drivers$COLOR[i] <- paste(palette[1, match(drivers$MODULE[i], mods)],
palette[2, match(drivers$MODULE[i], mods)],
palette[3, match(drivers$MODULE[i], mods)], collapse=" ")
## Process each module separately. Just perform for the 1st module:
i <- 1
rows <- which(drivers$MODULE == mods[i])
if(length(rows) > 0)
tmp <- kda2himmeli.exec(job.kda, valdata, drivers[rows,], mods, palette)

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
\code{\link{kda2himmeli}}
}

