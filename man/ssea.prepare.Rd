\name{ssea.prepare}
\alias{ssea.prepare}
\title{
Prepare an indexed database for MSEA
}
\description{
\code{ssea.prepare} prepares a database that includes hierarchical for 
modules, i.e. it collects gene list and unique marker list of the modules
for MSEA process
}
\usage{
ssea.prepare(job)
}
\arguments{
\item{job}{a data list with the following components: \preformatted{
modules: module identities as characters.
genesgene: identities as characters.
loci: marker identities as characters.
moddata: preprocessed module data (indexed identities).
gendata: preprocessed mapping data (indexed identities).
locdata: preprocessed marker data (indexed identities).
mingenes: minimum module size allowed.
maxgenes: maximum module size allowed.
maxoverlap: maximum module overlap allowed (1.0 to skip).
quantiles: quantile points for test statistic.
}
}
}
\details{
\code{ssea.prepare} removes extreme-sized modules, constructs a
hierarchical representation of genes and modules, obtains hit counts for
markers, and returns the finalized module, genes, markers, database 
information.
}
\value{
\item{job }{an updated data list with the following components:
\preformatted{
modules: finalized module names.
moddata: finalized module data.
gendata: finalized mapping data.
locdata: finalized marker data.
quantiles: verified quantile points.
database$modulesizes: gene counts for modules.
database$modulelengths: distinct markers counts for
modules.
database$moduledensities: ratio between distinct and
non-distinct markers.
database$genesizeslocus: count for each gene.
database$module2genes: gene lists for each module.
database$gene2locilocus: lists for each gene.
database$locus2row: indices in the marker data frame
for each marker.
database$observed: matrix of observed counts of values
that exceed each quantile point for each marker.
database$expected: 1.0 - quantile points.
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
\code{\link{ssea.finish}},  \code{\link{ssea.start}}, 
\code{\link{ssea2kda}}
}
