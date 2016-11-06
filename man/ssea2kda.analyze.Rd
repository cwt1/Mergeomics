\name{ssea2kda.analyze}
\alias{ssea2kda.analyze}
\title{
Apply second MSEA after merging the modules
}
\description{
\code{ssea2kda.analyze} performs a second MSEA for the updated modules 
after merging the highly overlapped modules (according to a specified 
overlapping ratio)
}
\usage{
ssea2kda.analyze(job, moddata)
}
\arguments{
\item{job}{the data list including the information of modules, genes, and 
markers, and also involving the database that uses indexed identities for 
modules, genes, and markers \code{(job$database)}.  
}
\item{moddata}{
merged modules including MODULE, GENE, and OVERLAP information
}
}
\details{
\code{ssea2kda.analyze} constructs new gene lists for merged modules and 
updates module database including module sizes, lengths, densities (based 
on marker sizes), and gene list. Then, it runs  a second MSEA and returns 
the enrichment scores of the updated module database. 
}
\value{
\item{res }{data list including updated information (after merge) such as, 
enrichment scores of merged modules}
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
## Collect genes and top markers from original files.
noddata <- ssea2kda.import(job.msea$genfile, job.msea$locfile)

## Select candidate modules (significant ones according to FDRs)
res <- job.msea$results
res <- res[order(res$P),]
rows <- which(res$FDR < 0.25)
res <- res[rows,]
    
## Collect member genes.
moddata <- job.msea$moddata
pos <- match(moddata$MODULE, res$MODULE)
moddata <- moddata[which(pos > 0),]

## Restore original identities.
modinfo <- job.msea$modinfo
modinfo$MODULE <- job.msea$modules[modinfo$MODULE]
moddata$MODULE <- job.msea$modules[moddata$MODULE]
moddata$GENE <- job.msea$genes[moddata$GENE]

## Merge and trim overlapping modules.
moddata$OVERLAP <- moddata$MODULE
rmax <- 0.33
moddata <- tool.coalesce(items=moddata$GENE, groups=moddata$MODULE,
rcutoff=rmax)
moddata$MODULE <- moddata$CLUSTER
moddata$GENE <- moddata$ITEM
moddata$OVERLAP <- moddata$GROUPS
moddata <- moddata[,c("MODULE", "GENE", "OVERLAP")]
moddata <- unique(moddata)
    
## Calculate enrichment scores for merged modules.
tmp <- unique(moddata[,c("MODULE","OVERLAP")])
res <- ssea2kda.analyze(job.msea, moddata)

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
\code{\link{ssea2kda}}
}
