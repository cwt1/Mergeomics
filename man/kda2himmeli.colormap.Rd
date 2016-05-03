\name{kda2himmeli.colormap}
\alias{kda2himmeli.colormap}
\title{
Assign one color to each unique module
}
\description{
\code{kda2himmeli.colormap} takes number of the modules and assigns a 
particular color to each module. Returns the color list (palette).
}
\usage{
kda2himmeli.colormap(ncolors)
}
\arguments{
\item{ncolors}{
number of the unique modules 
}
}
\value{
\item{palette }{color list: one color is assigned to each module}
}
\examples{
color.number = 5 ## let us assume we have 5 modules, assign 1 color to each:
palette <- kda2himmeli.colormap(color.number)
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
\code{\link{kda2himmeli}}
}

