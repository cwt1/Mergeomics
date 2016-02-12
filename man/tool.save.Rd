\name{tool.save}
\alias{tool.save}
\title{
Save a data frame in tab-delimited file
}
\description{
\code{tool.save} saves a given data frame into a specified file within a 
given directory.
}
\usage{
tool.save(frame, file, directory = NULL, verbose = TRUE, compression = FALSE)
}
\arguments{
\item{frame}{data frame to be saved into file}
\item{file}{name of the output file to be written}
\item{directory}{path of the directory for the file}
\item{verbose}{specifies whether the information about file saving process
will be displayed to user}
\item{compression}{specifies whether the file is compressed while saving.
Applicable for only UNIX-family systems with gzip.}
}
\value{
\item{fname }{returns file name with full path}
}
\examples{
aa<- data.frame(MODULE=c("Mod1", "Mod1", "Mod2", "Mod2", "Mod3"), 
NODE=c("GeneA", "GeneC", "GeneB", "GeneC", "GeneA"))
tool.save(aa, "aa.save.txt")
file.remove("aa.save.txt") ## delete the saved file!
}
\author{
Ville-Petteri Makinen 
}
\note{Compression only works on UNIX-family systems with gzip.}
