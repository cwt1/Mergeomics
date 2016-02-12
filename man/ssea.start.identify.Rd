\name{ssea.start.identify}
\alias{ssea.start.identify}
\title{
Convert identities to indices for MSEA
}
\description{
\code{\link{ssea.start.identify}} finds matching identities for the given
variable name. It searches the members of \code{dat} among the members of
\code{labels} with respect to the \code{varname} attribute, returns the 
matching rows of the \code{dat}.
}
\usage{
ssea.start.identify(dat, varname, labels)
}
\arguments{
\item{dat}{data list (source) of the identities that will be searched. 
e.g. the information after merging of overlapped genes (containing shared 
markers) }
\item{varname}{search and match will be performed with respect to 
which attribute (MODULE or NODE or MARKER) }
\item{labels}{the place, where the identities of \code{dat} will be 
searched and matched. }
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
bb <- ssea.start.identify(aa, "MODULE", c("Mod1"))
bb
cc <- ssea.start.identify(aa, "MODULE", c("Mod1", "Mod3"))
cc
dd <- ssea.start.identify(aa, "NODE", c("GeneA"))
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
\code{\link{ssea.start}}
}
