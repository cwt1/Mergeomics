\name{tool.read}
\alias{tool.read}
\title{
Read a data frame from a file
}
\description{
\code{tool.read} reads contents of given input file.
}
\usage{
tool.read(file, vars = NULL)
}
\arguments{
\item{file}{file name to be read}
\item{vars}{if we want to read particular attributes (columns) from the 
input file, we need to specify names of these attributes within list 
\code{vars} (attribute names can be e.g. MODULE, GENE, LOCUS, etc.)}
}
\details{
All lines with NAs are excluded.
}
\value{
\item{dat }{data frame including content of the given file. If \code{vars} 
is specified, only the listed columns inside the \code{vars} list will be
returned.}
}
\examples{
## read the network file as an example:
net.info <- tool.read(system.file("extdata","network.mouseliver.mouse.txt", 
package="Mergeomics"))
dim(net.info)
names(net.info)
}
\author{
Ville-Petteri Makinen 
}
