\name{tool.coalesce}
\alias{tool.coalesce}
\title{
Calculate overlaps between groups (main function)
}
\description{
\code{tool.coalesce} is utilized to merge and trim either overlapping 
modules (containing shared genes) or overlapping genes (containing 
shared markers)
}
\usage{
tool.coalesce(items, groups, rcutoff = 0, ncore = NULL)
}
\arguments{
\item{items}{
array of item identities
}
\item{groups}{
array of group identities for items
}
\item{rcutoff}{
maximum overlap not coalesced
}
\item{ncore}{
minimum number of items required for trimming
}
}
\value{a data list with the following components:
\item{CLUSTER }{cluster identities after merging and triming (a subset of 
group identities)}
\item{ITEM }{item identities}
\item{GROUPS }{comma separated overlapping group identities}  
}
\examples{
## read the coexpr module file as an example:
moddata <- tool.read(system.file("extdata",
"modules.mousecoexpr.liver.human.txt", package="Mergeomics"))

## let us find the overlapping ratio between first 10 modules in the file:
## to merge overlapping modules first collect member genes:
mod.names <- unique(moddata$MODULE)[1:min(length(unique(moddata$MODULE)),
10)]
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]

## Merge and trim overlapping modules.(max allowed overlap ratio is 0.33)
rmax <- 0.33
moddata$OVERLAP <- moddata$MODULE
moddata <- tool.coalesce(items=moddata$GENE, groups=moddata$MODULE,
rcutoff=rmax)
moddata$MODULE <- moddata$CLUSTER
moddata$GENE <- moddata$ITEM
moddata$OVERLAP <- moddata$GROUPS
moddata <- moddata[,c("MODULE", "GENE", "OVERLAP")]
moddata <- unique(moddata)
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
