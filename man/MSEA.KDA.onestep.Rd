\name{MSEA.KDA.onestep}
\alias{MSEA.KDA.onestep}
\title{
Run MSEA and/or KDA in one step
}

\description{
\code{\link{MSEA.KDA.onestep}} performs Marker Set Enrichment Analysis (MSEA)
and/or Key Driver Anlaysis (KDA) processes in one step.
}
\usage{
MSEA.KDA.onestep(plan, apply.MSEA=TRUE, apply.KDA=FALSE, 
maxoverlap.genesets=0.33, symbol.transfer.needed=FALSE, 
sym.from=c("HUMAN", "MOUSE"), sym.to=c("HUMAN", "MOUSE"))
}
\arguments{
\item{plan}{a data list including file and parameter settings for MSEA and/or
KDA processes: \preformatted{
label: unique identifier for the analysis
folder: output folder for results
modfile: path to module file (cols: MODULE GENE) 
genfile: path to gene file (cols: GENE LOCUS) (MSEA-specific)
marfile: path to marker file (cols: MARKER VALUE) (MSEA-specific)
inffile: path to module info file (cols: MODULE DESCR)
seed: seed for random number generator
permtype: gene for gene-level, locus for marker-level
nperm: max number of random permutations
mingenes: min number of genes per module (after merging)
maxgenes: max number of genes per module
quantiles: cutoffs for test statistic
maxoverlap: max overlap allowed between genes
netfile:  path to network file (TAIL HEAD WEIGHT) (KDA-specific)
}
}
\item{apply.MSEA}{determines whether MSEA will be performed to the given
set. Default value is TRUE.}
\item{apply.KDA}{determines whether KDA will be performed to the given
set. Default value is FALSE.}
\item{maxoverlap.genesets}{maximum overlapping ratio for the genesets. This is
applicable if KDA is performed following the MSEA process in one-step running.
Default value is 0.33.}
\item{symbol.transfer.needed}{determines whether gene symbols in the 
gene sets are needed to be transformed between different species.
Default value is FALSE.}
\item{sym.from}{defines the species, whose gene symbols will be converted
to the gene symbols of \code{sym.to} species. It can be either HUMAN or MOUSE.
It is applicable if \code{symbol.transfer.needed} is TRUE.}
\item{sym.to}{defines the species, whose gene symbols will be converted
from the gene symbols of \code{sym.from} species. It can be either HUMAN or
MOUSE. It is applicable if \code{symbol.transfer.needed} is TRUE.}
}
\details{
\code{\link{MSEA.KDA.onestep}} performs MSEA and/or KDA operations in one
built-in function. Users can run both MSEA and KDA sequentially, or they can
run either MSEA or KDA in one step with the same function. 
If MSEA and KDA will be applied sequentially, significantly enriched gene sets
(having FDR < 0.25), coming from MSEA results, will be merged if their 
overlapping ratios are larger than a given threshold, i.e. 
\code{maxoverlap.genesets}, to proceed the next step with relatively indepent
gene sets. Then, KDA is applied to this relatively independent gene sets.
}
\value{
\item{plan }{the updated data frame after performing MSEA and/or KDA.
If MSEA is performed, results will include standard MSEA results (see 
\code{\link{ssea.analyze}} for details); if KDA is applied, results will
include standard KDA results (see \code{\link{kda.analyze}} for details).}
}
\examples{
plan <- list()
plan$label <- "hdlc"
plan$folder <- "Results"
plan$genfile <- system.file("extdata", 
"genes.hdlc_040kb_ld70.human_eliminated.txt", package="Mergeomics")
plan$marfile <- system.file("extdata", 
"marker.hdlc_040kb_ld70.human_eliminated.txt", package="Mergeomics")
plan$modfile <- system.file("extdata", 
"modules.mousecoexpr.liver.human.txt", package="Mergeomics")
plan$inffile <- system.file("extdata", 
"coexpr.info.txt", package="Mergeomics")
plan$nperm <- 100 ## default value is 20000
\donttest{plan <- MSEA.KDA.onestep(plan, apply.MSEA=TRUE)}
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
\code{\link{ssea.analyze}}, \code{\link{kda.analyze}}
}
