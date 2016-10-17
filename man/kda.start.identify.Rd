\name{kda.start.identify}
\alias{kda.start.identify}
\title{
Convert identities to indices for wKDA
}
\description{
\code{\link{kda.start.identify}} searches the members of \code{dat} among
the members of \code{labels} with respect to the \code{varname} attribute, 
returns the matching rows of the \code{dat}.   
}
\usage{
kda.start.identify(dat, varname, labels)
}
\arguments{
\item{dat}{
data list of the identities that will be searched
}
\item{varname}{
search will be performed with respect to which attribute (MODULE or NODE)
}
\item{labels}{
the place, where data list (i.e. \code{dat}) will be searched
}
}
\value{
\item{res}{matched rows of \code{dat} among the members of \code{labels} 
list according to the \code{varname} attribute}
}
\examples{
## Converts identities (either module names or gene names) to the indices
aa<- data.frame(MODULE=c("Mod1", "Mod1", "Mod2", "Mod2", "Mod3"), 
NODE=c("GeneA", "GeneC", "GeneB", "GeneC", "GeneA"))
aa
bb <- kda.start.identify(aa, "MODULE", c("Mod1"))
bb
cc <- kda.start.identify(aa, "MODULE", c("Mod1", "Mod3"))
cc
dd <- kda.start.identify(aa, "NODE", c("GeneA"))
dd
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
\code{\link{kda.analyze}}, \code{\link{kda.finish}}, 
\code{\link{kda.prepare}}, \code{\link{kda.start}} 
}

