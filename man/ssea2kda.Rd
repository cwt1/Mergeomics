\name{ssea2kda}
\alias{ssea2kda}
\title{
Generate inputs for wKDA
}
\description{
\code{ssea2kda} forwards MSEA results to weighted key driver analysis (wKDA) 
from the first MSEA results, merges the overlapped modules according to 
a given overlapping ratio to obtain a relatively independent module set,
apply a second MSEA on the merged modules (supersets), updates and saves 
the second MSEA results properly for wKDA process.
}
\usage{
ssea2kda(job, symbols = NULL, rmax = NULL, min.module.count=NULL)
}
\arguments{
\item{job}{data list including the organized results of MSEA process. It 
has following components:\preformatted{
results: updated information of modules including: 
number of distinct member genes (NGENES), 
number of distinct member markers (NLOCI), 
ratio of distinct to non-distinct markers (DENSITY), 
false discovery rates (FDR).
generesults: updated gene-specific information including: 
indexed gene identity (GENE), 
gene size (NLOCI), 
unadjusted enrichment score (SCORE),  
marker with max value (LOCUS), 
marker value (VALUE).
}
}
\item{symbols}{
dataframe for translating gene symbols
}
\item{rmax}{
maximum allowed overlap ratio between gene sets
}
\item{min.module.count}{
minimum number of the pathways to be taken from the MSEA results to merge.
Default value is NULL. If it is not specified, all the pathways having MSEA-FDR 
value less than 0.25 will be considered for merging if they are overlapping 
with the given ratio \code{rmax}. 
}
}
\details{
\code{ssea2kda} gets genes and top markers from input files, selects 
significant modules with respect to ordered p-values, gets identities of 
modules and genes, merges and trims the overlapping modules (either having FDR
less than 0.25 or top \code{min.module.count} modules when ranked up to the
P-values), obtains enrichment scores for merged modules, translates the gene 
symbols (between species) if needed, and finally saves the module, gene, node, 
and marker information into relevant output files. 
}
\value{
\item{plan }{an updated data list with the following components: 
\preformatted{
label: unique identifier for the analysis.
parent: parent folder for results.
modfile: path of module file (columns: MODULE NODE).
inffile: path of module information file 
(columns: MODULE DESCR).
nodfile: path of node selection file (columns: NODE).
}
}
}
\examples{
job.msea <- list()
job.msea$label <- "hdlc"
job.msea$folder <- "Results"
job.msea$genfile <- system.file("extdata", 
"genes.hdlc_040kb_ld70.human_eliminated.txt", package="Mergeomics")
job.msea$marfile <- system.file("extdata", 
"marker.hdlc_040kb_ld70.human_eliminated.txt", package="Mergeomics")
job.msea$modfile <- system.file("extdata", 
"modules.mousecoexpr.liver.human.txt", package="Mergeomics")
job.msea$inffile <- system.file("extdata", 
"coexpr.info.txt", package="Mergeomics")
job.msea$nperm <- 100 ## default value is 20000

## ssea.start() process takes long time while merging the genes sharing high
## amounts of markers (e.g. loci). it is performed with full module list in
## the vignettes. Here, we used a very subset of the module list (1st 10 mods
## from the original module file) and we collected the corresponding genes
## and markers belonging to these modules:
moddata <- tool.read(job.msea$modfile)
gendata <- tool.read(job.msea$genfile)
mardata <- tool.read(job.msea$marfile)
mod.names <- unique(moddata$MODULE)[1:min(length(unique(moddata$MODULE)),
10)]
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
gendata <- gendata[which(!is.na(match(gendata$GENE, 
unique(moddata$GENE)))),]
mardata <- mardata[which(!is.na(match(mardata$MARKER, 
unique(gendata$MARKER)))),]

## save this to a temporary file and set its path as new job.msea$modfile:
tool.save(moddata, "subsetof.coexpr.modules.txt")
tool.save(gendata, "subsetof.genfile.txt")
tool.save(mardata, "subsetof.marfile.txt")
job.msea$modfile <- "subsetof.coexpr.modules.txt"
job.msea$genfile <- "subsetof.genfile.txt"
job.msea$marfile <- "subsetof.marfile.txt"
## run ssea.start() and prepare for this small set: (due to the huge runtime)
job.msea <- ssea.start(job.msea)
job.msea <- ssea.prepare(job.msea)
job.msea <- ssea.control(job.msea)
job.msea <- ssea.analyze(job.msea)
job.msea <- ssea.finish(job.msea)

###############  Create intermediary datasets for KDA ##################
syms <- tool.read(system.file("extdata", "symbols.txt", 
package="Mergeomics"))
syms <- syms[,c("HUMAN", "MOUSE")]
names(syms) <- c("FROM", "TO")
job.ssea2kda <- ssea2kda(job.msea, symbols=syms)

## Remove the temporary files used for the test:
file.remove("subsetof.coexpr.modules.txt")
file.remove("subsetof.genfile.txt")
file.remove("subsetof.marfile.txt")
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
\code{\link{ssea.analyze}}, \code{\link{kda.analyze}}
}
