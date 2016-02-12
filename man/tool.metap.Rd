\name{tool.metap}
\alias{tool.metap}
\title{
Estimate meta P-values
}
\description{
\code{tool.metap} returns the meta p-values of given datasets with 
multiple p-values.
}
\usage{
tool.metap(datasets, idcolumn, pcolumn, weights = NULL)
}
\arguments{
\item{datasets}{data list, whose meta p-values will be obtained}
\item{idcolumn}{column number of the datasets that includes identities}
\item{pcolumn}{column number of the datasets that includes p-values}
\item{weights}{weight list of the data list}
}
\value{
\item{res }{data list including identities and meta p-values of the 
given \code{datasets}}
}
\examples{
set.seed(1)
## let us assume we have p-values for the coexpr modules obtained from 
## distinct analyses by using different gene-marker mapping sets (e.g. eQTLs
## from diff tissues) and we would like to make a meta-analysis for 
## these multiple Pvalues of the modules:
datasets=list()
## we have 3 datasets and 3 diff result sets
datasets[[1]] <- data.frame(MODULE=c("Mod1", "Mod2", "Mod3", "Mod4"), 
P=c(rnorm(4)))
datasets[[2]] <- data.frame(MODULE=c("Mod1", "Mod2", "Mod3", "Mod4"), 
P=c(rnorm(4)))
datasets[[3]] <- data.frame(MODULE=c("Mod1", "Mod2", "Mod3", "Mod4"),
P=c(rnorm(4)))
idcolumn <- "MODULE" ## identifiers of the modules are in the 1st col
pcolumn <- "P" ## p values of the modules are in the 2nd col
tool.metap(datasets, idcolumn, pcolumn)
}
\author{
Ville-Petteri Makinen 
}
