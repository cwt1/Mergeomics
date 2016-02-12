\name{ssea.prepare.structure}
\alias{ssea.prepare.structure}
\title{
Construct hierarchical representation of components 
}
\description{
\code{ssea.prepare.structure} represents modules, genes, and markers in a 
hierarchical structure.  
}
\usage{
ssea.prepare.structure(moddata, gendata, nmods, ngens)
}
\arguments{
\item{moddata}{module data (indexed identities)}
\item{gendata}{mapping data (indexed identities)}
\item{nmods}{number of modules}
\item{ngens}{number of all genes}
}
\details{
\code{ssea.prepare.structure} finds member genes of modules and marker 
lists of genes; counts distinct markers within each module and obtains 
module's density from this count; at the end, it returns hierarchically 
structured results.
}
\value{
\item{res }{a data list with the following components: \preformatted{
modulesizes: module size
modulelengths: module length
moduledensities: module densities
genesizes: gene sizes of module
module2genesgene: list of module
gene2loci: markers lists of genes
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

## Remove extremely big or small modules:
st <- tool.aggregate(job.msea$moddata$MODULE)
mask <- which((st$lengths >= job.msea$mingenes) & 
(st$lengths <= job.msea$maxgenes))
pos <- match(job.msea$moddata$MODULE, st$labels[mask])
job.msea$moddata <- job.msea$moddata[which(pos > 0),]

## Construct hierarchical representation for modules, genes, and markers:
ngens <- length(job.msea$genes) 
nmods <- length(job.msea$modules)
db <- ssea.prepare.structure(job.msea$moddata, job.msea$gendata, 
nmods, ngens)

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
\code{\link{ssea.prepare}}
}
