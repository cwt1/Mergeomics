\name{kda.finish.estimate}
\alias{kda.finish.estimate}
\title{
Estimate measures for accomplished wKDA results
}
\description{
Estimates additional measures based on overlapping of module member nodes 
with hub neighbor nodes in the graph. 
}
\usage{
kda.finish.estimate(job)
}
\arguments{
\item{job}{
The data list that was subjected to wKDA. It involves the modules, 
member genes belonging to each module, graph information of the 
dataset, hubs and hubnets of the graph. 
}
}
\details{
\code{\link{kda.finish.save}} determines the overlaps of modules with hub 
neighborhoods, obtains graph measures based on the ratio of the observed 
overlap amounts to the expected overlap amount, and returns the values of 
this measure. 
}
\value{
\item{res }{Returns the overlapping ratio of the modules with hubnets. }
}
\examples{
## get the prepared and KDA applied dataset:(see kda.analyze for details)
data(job_kda_analyze)
## finish the KDA process by estimating additional measures for the modules
## such as module sizes, overlaps with hub neighborhoods, etc.
# job.kda <- kda.finish(job.kda)
# if (nrow(job.kda$results)==0){
# cat("No Key Driver Found!!!!")
# } else{
##  Estimate additional measures - see kda.analyze and kda.finish for details
#   res <- kda.finish.estimate(job.kda)
# }
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
\code{\link{kda.finish}}, \code{\link{kda.finish.save}}, 
\code{\link{kda.finish.summarize}}, \code{\link{kda.finish.trim}}
}
