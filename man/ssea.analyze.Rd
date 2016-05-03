\name{ssea.analyze}
\alias{ssea.analyze}
\title{
Marker set enrichment analysis (MSEA) 
}
\description{
\code{\link{ssea.analyze}} finds the enrichment of the pathways or 
co-expression modules by a marker set (e.g. associated risk variants 
-loci- of a relevant disease). 
Association study by mapping markers (e.g. SNPs) to genes (e.g. via 
expression QTLs). 
Enrichment P-values obtained by the MSEA denote the degree of enrichment 
of significantly disease-associated (high ranking) markers (e.g. eSNPs)
within these pathways when compared to the null distribution of expected 
uniform distribution of all ranks of the markers.
MSEA is performed with either gene-level or marker-level permutations 
based on Gaussian distribution. 
}
\usage{
ssea.analyze(job)
}
\arguments{
\item{job}{
the data list including fields: seed for random number generator 
\code{(job$seed)} (to obtain the same results from permutation when
the same input set is given), random permutation level 
\code{(job$permtype)} (either gene- or marker-level), maximum number 
of permutations \code{(job$nperm)} for the permutation test, and the
database that uses indexed identities for modules, genes, and markers
(e.g. loci) \code{(job$database)}.  
}
}
\details{
\code{\link{ssea.analyze}} associates the gene sets (pathways or 
co-expression modules) with relevant disease (e.g. Coronary Artery Disease)
association data by mapping markers (e.g. SNPs) to genes (e.g. via expression
QTLs). 
It performs the MSEA by using observed and estimated enrichment scores. 
First, the observed enrichment scores of the pathways by markers (e.g. 
loci) are calculated. Then, a Gaussian distribution based simulation is 
performed, by using the statistics of the observed scores (mean, std.dev.,
etc.), to obtain the estimated enrichment scores, enrichment frequencies,
and other statistics e.g. p-values for the pathways.
}
\value{
\item{job }{the updated data frame including results: indexed module 
identity, enrichment P-values, raw frequencies (raw frequency of a gene 
set defines the number of the estimated enrichment scores that are larger 
than this gene set's enrichment score under the null distribution based 
on Gaussian function)}.
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
\code{\link{ssea.control}}, \code{\link{ssea.finish}}, 
\code{\link{ssea.prepare}}, \code{\link{ssea.start}}, 
\code{\link{ssea2kda}}
}
