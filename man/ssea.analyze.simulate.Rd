\name{ssea.analyze.simulate}
\alias{ssea.analyze.simulate}
\title{
Simulate scores for MSEA 
}
\description{
\code{\link{ssea.analyze.simulate}} simulates enrichment scores by
randomly permuting database with respect to the specified permutation 
type (either gene-level or marker-level).
}
\usage{
ssea.analyze.simulate(db, observ, nperm, permtype)
}
\arguments{
\item{db}{
database including the indexed identities for modules, genes and markers 
(e.g. loci): \preformatted{
modulesizes: gene counts for modules.
modulelengths: distinct marker counts for modules.
moduledensities: ratio between distinct and 
non-distinct markers.
genesizes: marker count for each gene.
module2genes: gene lists for each module.
gene2loci: marker lists for each gene.
locus2row: row indices in the marker data frame for
each marker.
observed: matrix of observed counts of values that
exceed each quantile point for each marker.
expected: 1.0 - quantile points.
}
} 
\item{observ}{
observed enrichment scores
}
\item{nperm}{
maximum nubmer of permutations (for simulation)
}
\item{permtype}{
permutation type (either gene or locus)
}
}
\value{
\item{scoresets }{simulated score lists for the statistically significant 
modules}
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

## Observed enrichment scores.
db <- job.msea$database
scores <- ssea.analyze.observe(db)
nmods <- length(scores)

## Simulated scores.
nperm <- job.msea$nperm
nullsets <- ssea.analyze.simulate(db, scores, nperm, job.msea$permtype)

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
\code{\link{ssea.analyze}}
}
