\name{tool.coalesce.merge}
\alias{tool.coalesce.merge}
\title{
Merge overlapping clusters
}
\description{
\code{tool.coalesce.merge} determines combinable groups and trims clusters 
by removing rarest items.
}
\usage{
tool.coalesce.merge(data, ncore)
}
\arguments{
\item{data }{data list including following components: \preformatted{
CLUSTER: cluster label
NODE: item (node) name
}
}
\item{ncore}{
minimum number of items required for trimming
}
}
\value{
\item{res }{data list including GROUPS, ITEMs, and their hit COUNTs}
}
\examples{
## Generate item and group labels for 100 items:
## Assume that unique gene number (items) is 60:
members <- 1:100 ## will be updated
modules <- 1:100 ## will be updated
set.seed(1)
for (i in 1:10){
## each time pick 10 items (genes) from 60 unique item labels
members[(i*10-9):(i*10)] <- sample(60,10) 
}
## Assume that unique group labels is 30:
for (i in 1:10){
## each time pick 10 items (genes) from 30 unique group labels
modules[(i*10-9):(i*10)] <- sample(30, 10)
}
rcutoff <- 0.33
ncore <- length(members)
## Default output.
res <- data.frame(CLUSTER=modules, GROUPS=modules, ITEM=members,
stringsAsFactors=FALSE)
## Iterative merging and trimming.
res$COUNT <- 0.0
while(TRUE) {
clust <- tool.coalesce.find(res, rcutoff)    
if(is.null(clust)) break
res <- tool.coalesce.merge(clust, ncore)
}
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
