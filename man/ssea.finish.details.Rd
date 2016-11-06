\name{ssea.finish.details}
\alias{ssea.finish.details}
\title{
Organize and save module, gene, top locus, Ps of MSEA results
}
\description{
\code{ssea.finish.details} finds significant modules and their gene lists, 
and top marker (with GWAS -log10 transformed p-vals) of these genes, 
merge results of markers, genes and module statistics, sort results 
according to first, module enrichment score, then marker P-value, 
and saves these sorted results into the relevant files.
}
\usage{
ssea.finish.details(job)
}
\arguments{
\item{job}{
data list including the results of MSEA process: 
\preformatted{
label: unique identifier for the analysis.
folder: output folder for results.
modinfo: descriptions of the modules.
resultsdata: frame including indexed module identities
(MODULE) and enrichment P-values (P).
database: database including indexed identities for 
modules, genes, and markers.
}
}
}
\value{
None.
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

## Estimate mod FDR values, sort according to significance, save full results:
job.msea <- ssea.finish.fdr(job.msea)
## Collect top markers(e.g.loci) within genes, save genes with top marker Pval
job.msea <- ssea.finish.genes(job.msea)
## Find signficant modules, collect gene members of top modules,
## Merge gene results (with top marker info), 
## Sort and save details according to enrichment and marker value:
job.msea <- ssea.finish.details(job.msea)

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
\code{\link{ssea.finish}}
}
