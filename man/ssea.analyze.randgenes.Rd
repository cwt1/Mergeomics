\name{ssea.analyze.randgenes}
\alias{ssea.analyze.randgenes}
\title{
Estimate enrichment from randomized genes
}
\description{
\code{\link{ssea.analyze.randgenes}} simulates enrichment scores by
randomizing the genes from all modules (from database - db)
}
\usage{
ssea.analyze.randgenes(db, targets, gene_sel)
}
\arguments{
\item{db}{
database including the indexed identities for modules, genes and markers: 
\preformatted{
modulesizes: gene counts for modules.
modulelengths: distinct marker counts for modules.
moduledensities: ratio between distinct and non-distinct
markers.
genesizes: marker count for each gene.
module2genes: gene lists for each module.
gene2loci: marker lists for each gene .
locus2row: row indices in the marker data frame for each
marker.
observed: matrix of observed counts of values that exceed
each quantile point for each marker.
expected: 1.0 - quantile points.
}
} 
\item{targets}{all modules}
\item{gene_sel}{selected genes to be trimmed away to avoid signal inflation
 of null background in gene permutation.
}
}
\value{
\item{scores }{randomly simulated enrichment scores}
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
gene2loci <- db$gene2loci
locus2row <- db$locus2row
observed <- db$observed
#Calcuate individual gene enrichment score
trim_scores <- rep(NA, length(gene2loci))
for(k in 1:length(trim_scores)) {
  genes <- k
  # Collect markers.
  loci <- integer()
  for(i in genes) 
    loci <- c(loci, gene2loci[[i]])
        
  # Determine data rows.
  loci <- unique(loci)
  rows <- locus2row[loci]
  nloci <- length(rows)
        
  # Calculate total counts.
  e <- (nloci/length(locus2row))*colSums(observed)
  o <- observed[rows,]
  if(nloci > 1) o <- colSums(o)
        
  # Estimate enrichment.
  trim_scores[k] <- ssea.analyze.statistic(o, e)
}
trim_start=0.002 # default
trim_end=1-trim_start
cutoff=as.numeric(quantile(trim_scores,probs=c(trim_start,trim_end)))
gene_sel=which(trim_scores>cutoff[1]&trim_scores<cutoff[2])

scores <- ssea.analyze.observe(db)
nmods <- length(scores)

## Simulated scores.
nperm <- job.msea$nperm
observ <- scores
## Include only non-empty modules for simulation.
nmods <- length(db$modulesizes)
targets <- which(db$modulesizes > 0)
hits <- rep(NA, nmods)
hits[targets] <- 0
    
## Prepare data structures to hold null samples.
keys <- rep(0, nperm)
scores <- rep(NA, nperm)
scoresets <- list()
for(i in 1:nmods) scoresets[[i]] <- double()
## Simulate random scores.
## within a for loop: check capacity, find new statistics, update snull  
## distribution (simulated null distr.) by permuting genes
snull <- ssea.analyze.randgenes(db, targets, gene_sel)

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
\code{\link{ssea.analyze}}
}
