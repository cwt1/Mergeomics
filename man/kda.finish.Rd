\name{kda.finish}
\alias{kda.finish}
\title{
Organize and save results 
}
\description{
After wKDA process is accomplished, \code{\link{kda.finish.estimate}} sums 
up the results and log them to the relevant files and folders. Besides, 
return them within the given job parameter.  
}
\usage{
kda.finish(job)
}
\arguments{
\item{job}{
the data list including label and folder fields to specify a unique 
identifier for the wKDA process and the output folder for the obtained
results, respectively.  
}
}
\details{
\code{\link{kda.finish.estimate}} estimates additional measures if needed,
saves results into relevant files, trims numbers to provide a simpler file 
for viewing, and stores a summary file of top hits after the wKDA prcess is 
accomplished. It also obtains the overlaps of the modules with hub 
neighborhoods, finds co-hubs information, determines the top key driver 
for each module and saves the updated and sorted p-values belonging to 
them. 
}
\value{
\item{job}{
updated information including the overlapping hub neighborhoods, 
co-hubs information, top driver of each module, and their updated and 
sorted p-values. 
}
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

## finish the KDA process
job.kda <- kda.finish(job.kda)

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
\code{\link{kda.finish.estimate}}, \code{\link{kda.finish.save}},
\code{\link{kda.finish.summarize}}, \code{\link{kda.finish.trim}}
}
