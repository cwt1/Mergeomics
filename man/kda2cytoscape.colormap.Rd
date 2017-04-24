\name{kda2cytoscape.colormap}
\alias{kda2cytoscape.colormap}
\title{
Assign one color to each unique module
}
\description{
\code{kda2cytoscape.colormap} takes number of the modules and assigns a 
particular color to each module. Returns the color list (palette).
}
\usage{
kda2cytoscape.colormap(ncolors)
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
palette <- kda2cytoscape.colormap(color.number)
}
\references{
Shu L, Zhao Y, Kurt Z, Byars SG, Tukiainen T, Kettunen J, Orozco LD, 
Pellegrini M, Lusis AJ, Ripatti S, Zhang B, Inouye M, Makinen V-P, Yang X.
Mergeomics: multidimensional data integration to identify pathogenic 
perturbations to biological systems. BMC genomics. 2016;17(1):874.
}
\author{
Zeyneb Kurt
}
\seealso{
\code{\link{kda2cytoscape}}
}

