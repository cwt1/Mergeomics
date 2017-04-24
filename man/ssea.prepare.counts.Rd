\name{ssea.prepare.counts}
\alias{ssea.prepare.counts}
\title{
Calculate hit counts up to a given quantile
}
\description{
Counts unique loci in a module, maps the marker data of a module
to the all available markers by creating a bit matrix for values above the 
given quantiles. Created bit matrix contains either TRUE (above quantiles)
or FALSE (below or equals to quantiles) values as a resuls of these 
comparisons. It returns the results (marker mapping and bit matrix)
}
\usage{
ssea.prepare.counts(locdata, nloci, quantiles)
}
\arguments{
\item{locdata}{marker data}
\item{nloci}{number of elements in markers list}
\item{quantiles}{quantile points for test statistic}
}
\value{
\item{res }{a data list with the following components: \preformatted{
locus2row: mapped marker information
observed: bit matrix that involves TRUEs and FALSEs 
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
## Determine test cutoffs:
if(is.null(job.msea$quantiles)) {
lengths <- db$modulelengths
mu <- median(lengths[which(lengths > 0)])
job.msea$quantiles <- seq(0.5, (1.0 - 1.0/mu), length.out=10)
}
## Calculate hit counts:
nloci <- length(job.msea$loci)
hits <- ssea.prepare.counts(job.msea$locdata, nloci, job.msea$quantiles)
db <- c(db, hits)

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
\code{\link{ssea.prepare}}
}
