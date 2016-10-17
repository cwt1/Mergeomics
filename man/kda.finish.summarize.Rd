\name{kda.finish.summarize}
\alias{kda.finish.summarize}
\title{
Summarize the wKDA results
}
\description{
Create a summary file of top key drivers. The file includes the key driver 
of each block of the dataset and their p-values.}
\usage{
kda.finish.summarize(res, job)
}
\arguments{
\item{res}{
the data frame including the p-values, false discovery rates, and fold 
scores of the nodes obtained from \code{\link{kda.finish.trim}} 
}
\item{job}{
the data frame including the path of output file which will briefly 
contain top key drivers of the blocks and ranked p-values of those top
key drivers
}
}
\details{
\code{\link{kda.finish.summarize}} determines the ranking scores of blocks, 
finds the top node for each block, selects and saves top key drivers, and 
stores P-values into file. top drovers of the blocks are also returned to 
the user.
}
\value{
\item{res}{data frame including top node for each block}
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
##  Create a simpler file for viewing by trimming floating numbers
#   res <- kda.finish.trim(res, job.kda)
##  Create a summary file of top hit KDs.
#   res <- kda.finish.summarize(res, job.kda)
# }
## See kda.analyze() and kda.finish() for details
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
\code{\link{kda.finish}}, \code{\link{kda.finish.estimate}}, 
\code{\link{kda.finish.save}}, \code{\link{kda.finish.trim}}
}

