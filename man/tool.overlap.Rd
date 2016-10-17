\name{tool.overlap}
\alias{tool.overlap}
\title{
Calculate overlaps between groups of specified items
}
\description{
\code{tool.overlap} checks each pair of blocks, finds number of shared 
items, and obtains significance values of the sharings for block pairs.
}
\usage{
tool.overlap(items, groups, nbackground = NULL)
}
\arguments{
\item{items}{array of item identities}
\item{groups}{array of group identities for items}
\item{nbackground}{total number of items}
}
\value{a data list including following components
\item{A }{group name}
\item{B }{group name}
\item{POSa }{group name rank}
\item{POSb }{group name rank}
\item{Na }{group A size}
\item{Nb }{group B size}
\item{Nab }{shared items}
\item{R }{overlap ratio}
\item{F }{fold change to null expectation}
\item{P }{overlap P-value (Fisher's test)}
}
\examples{
## read the coexpr module file as an example:
moddata <- tool.read(system.file("extdata",
"modules.mousecoexpr.liver.human.txt", package="Mergeomics"))
## let us find the overlapping ratio between two modules:
## pick the first and last modules:
mod.names <- unique(moddata$MODULE)[c(1,length(unique(moddata$MODULE)))]
if(length(mod.names) > 0){
modA.members <- moddata[which(moddata$MODULE == mod.names[1]),]
modB.members <- moddata[which(moddata$MODULE == mod.names[2]),]
}
mod.pool <- rbind(modA.members, modB.members)
overlap.stats <- tool.overlap(mod.pool[,2], mod.pool[,1])
}
\author{
Ville-Petteri Makinen 
}
