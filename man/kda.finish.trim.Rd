\name{kda.finish.trim}
\alias{kda.finish.trim}
\title{
Trim numbers before save
}
\description{
\code{\link{kda.finish.trim}} trims p-values, false discovery rates, and 
fold scores to make them nicer to look at before saving the file. It also 
returns trimmed results to the user.
}
\usage{
kda.finish.trim(res, job)
}
\arguments{
\item{res}{
includes p-values, false discovery rates, and fold scores of the nodes
}
\item{job}{
data frame including output folder path to store trimmed results
}
}
\value{
\item{res }{
Trimmed and formatted p-values, false discovery rates, and fold 
scores of the nodes
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
##  Create a simpler file for viewing by trimming floating numbers
#   res <- kda.finish.trim(res, job.kda)
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
\code{\link{kda.finish.save}}, \code{\link{kda.finish.summarize}}
}
