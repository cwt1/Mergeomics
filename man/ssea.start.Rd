\name{ssea.start}
\alias{ssea.start}
\title{
Create a job for MSEA
}
\description{
Creates identities (for modules, member genes, and loci) to start MSEA 
process.
}
\usage{
ssea.start(plan)
}
\arguments{
\item{plan}{a data list with the following components: \preformatted{
label: unique identifier for the analysis
folder: output folder for results
modfile: path to module file (cols: MODULE GENE)
marfile: path to marker file (cols: MARKER VALUE)
genfile: path to gene file (cols: GENE LOCUS)
inffile: path to module info file (cols: MODULE DESCR)
seed: seed for random number generator
permtype: gene for gene-level, locus for marker-level
nperm: max number of random permutations
mingenes: min number of genes per module (after merging)
maxgenes: max number of genes per module
quantiles: cutoffs for test statistic
maxoverlap: max overlap allowed between genes
}
}
}
\details{
\code{ssea.start} imports modules, genes-locus mapping, and locus values; 
removes the genes with no locus values from the list, find identities for
modules, genes, loci components, and excludes missing data and factorize 
identities for these components.
}
\value{
\item{job }{a data list with the following components: \preformatted{
modules: module identities as characters.
genes: gene identities as characters.
loci: marker identities as characters.
moddata: preprocessed module data (indexed identities)
modinfo: description of the modules. 
gendata: preprocessed mapping data between genes and 
markers (indexed identities).
locdata: preprocessed marker data (indexed identities)
geneclusters: genes with shared markers.
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
## run ssea.start() for this small set:(due to the huge runtime we did not use
## full sets of modules, genes, and markers)
job.msea <- ssea.start(job.msea)

## Remove the temporary files used for the test:
file.remove("subsetof.coexpr.modules.txt")
file.remove("subsetof.genfile.txt")
file.remove("subsetof.marfile.txt")
}
\references{
Shu L, Zhao Y, Kurt Z, Byars SG, Tukiainen T, Kettunen J, Orozco LD, 
Pellegrini M, Lusis AJ, Ripatti S, Zhang B, Inouye M, Makinen V-P, Yang X.
Mergeomics: multidimensional data integration to identify pathogenic 
perturbations to biological systems. BMC genomics. 2016;17(1):874.
}
\author{
Ville-Petteri Makinen 
}
\seealso{
\code{\link{ssea.analyze}}, \code{\link{ssea.control}}, 
\code{\link{ssea.finish}}, \code{\link{ssea.prepare}},
\code{\link{ssea2kda}}
}
