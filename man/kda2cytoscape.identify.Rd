\name{kda2cytoscape.identify}
\alias{kda2cytoscape.identify}
\title{
Match identities with respect to given variable name
}
\description{
\code{kda2cytoscape.identify} searches the given data list \code{dat}
within the \code{labels} according to the specified attribute (variable 
name). It returns the matched rows. Hence, it finds identifier numbers 
for the searched data list \code{dat}.
}
\usage{
kda2cytoscape.identify(dat, varname, labels)
}
\arguments{
\item{dat}{
node ID list whose symbols or names will be collected from network 
node name (or symbol) list.
}
\item{varname}{
specifies that \code{dat} will be searched among \code{labels} 
according to which variable (attribute). Here, gene symbols whose IDs 
are given, will be searched in the causal network node list according 
to the NODE attribute.
}
\item{labels}{
the data list possibly including names or symbols corresponding to 
the given IDs in the \code{dat} data list. 
}
}
\value{
\item{res }{the matching rows of \code{labels} with the identifiers of 
given data list \code{dat} }
}
\examples{
## Converts identities (either module names or gene names) to the indices
aa<- data.frame(MODULE=c("Mod1", "Mod1", "Mod2", "Mod2", "Mod3"), 
NODE=c("GeneA", "GeneC", "GeneB", "GeneC", "GeneA"))
aa
bb <- kda2cytoscape.identify(aa, "MODULE", c("Mod1"))
bb
cc <- kda2cytoscape.identify(aa, "MODULE", c("Mod1", "Mod3"))
cc
dd <- kda2cytoscape.identify(aa, "NODE", c("GeneA"))
dd
}
\references{
Shu L, Zhao Y, Kurt Z, Byars S, Tukiainen T, Kettunen J, Ripatti S, Zhang B,
Inouye M, Makinen VP, Yang X. Mergeomics: integration of diverse genomics
resources to identify pathogenic perturbations to biological systems.
bioRxiv doi: http://dx.doi.org/10.1101/036012
}
\author{
Zeyneb Kurt
}
\seealso{
\code{\link{kda2cytoscape}}
}

