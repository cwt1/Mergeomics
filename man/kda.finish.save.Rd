\name{kda.finish.save}
\alias{kda.finish.save}
\title{
Save full wKDA results
}
\description{
\code{\link{kda.finish.save}} sorts (according to KD p-values) and saves 
the wKDA results into specified files and folders. 
}
\usage{
kda.finish.save(res, job)
}
\arguments{
\item{res}{
the results obtained from \code{\link{kda.finish.estimate}}. They
will be stored into specified folder.
}
\item{job}{
information including the entire graph, nodes, modules, co-hubs, 
top key driver of each module, and their updated and sorted p-values. 
All the information included \code{job}  will be stored into relevant 
files.
}
}
\value{
\item{res }{
the results obtained from \code{\link{kda.finish.estimate}} are merged 
with the module and graph nodes information gained from \code{job} 
data frame. At the end, merged information is both written to file and 
returned to the user.
}
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
##  Save full results about modules such as co-hub, nodes, P-values info etc.
#   res <- kda.finish.save(res, job.kda)
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
\code{\link{kda.finish}}, \code{\link{kda.finish.estimate}}, 
\code{\link{kda.finish.summarize}}, \code{\link{kda.finish.trim}}
}

