\name{tool.translate}
\alias{tool.translate}
\title{
Translate gene symbols
}
\description{
\code{tool.translate} converts the symbols given in the list \code{from}
into the list \code{to}. e.g. we can translate human gene symbols into the 
mouse orthologs (or vice versa) if the symbol mapping file is provided.
}
\usage{
tool.translate(words, from, to)
}
\arguments{
\item{words}{translation table including words (i.e. gene symbols) that 
will be translated}
\item{from}{a list denoting the words will be translated from which symbols}
\item{to}{a list denoting the words will be translated to which symbols}
}
\value{
\item{words }{translated table (words)}
}
\examples{
syms <- tool.read(system.file("extdata", "symbols.txt", 
package="Mergeomics"))
syms <- syms[,c("HUMAN", "MOUSE")]
names(syms) <- c("FROM", "TO")
moddata <- tool.read(system.file("extdata", 
"modules.mousecoexpr.liver.human.txt", package="Mergeomics"))
moddata$NODE <- moddata$GENE
moddata$NODE <- tool.translate(words=moddata$NODE, from=syms$FROM, 
to=syms$TO)
}

\author{
Ville-Petteri Makinen 
}
