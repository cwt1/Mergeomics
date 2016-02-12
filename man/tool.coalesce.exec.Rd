\name{tool.coalesce.exec}
\alias{tool.coalesce.exec}
\title{
Find, merge, and trim overlapping clusters
}
\description{
\code{tool.coalesce.exec} searchs overlaps, iteratively merges and trims 
overlapping clusters (by using \code{\link{tool.coalesce.find}} and 
\code{\link{tool.coalesce.merge}}, respectively) until no more overlap is
available, and assigns representative label for the merged clusters.
}
\usage{
tool.coalesce.exec(items, groups, rcutoff, ncore)
}
\arguments{
\item{items}{array of item identities}
\item{groups}{array of group identities for items}
\item{rcutoff}{maximum overlap not coalesced}
\item{ncore}{minimum number of items required for trimming}
}
\value{a data list with the following components:
\item{CLUSTER }{cluster identities after merging and triming (a subset of 
group identities)}
\item{GROUPS }{comma separated overlapping group identities}
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
## Find and trim clusters after iteratively merging the overlapping ones:
res <- tool.coalesce.exec(members, modules, rcutoff, ncore)
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

