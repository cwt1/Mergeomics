\name{ssea.start.relabel}
\alias{ssea.start.relabel}
\title{
Update gene symbols after merging overlapped markers
}
\description{
\code{ssea.start.relabel} updates gene symbols within the modules after 
merging overlapping genes that contain shared markers
}
\usage{
ssea.start.relabel(dat, grp)
}
\arguments{
\item{dat}{
module data corresponding gene sets
}
\item{grp}{ 
gene data that is needed to be relabeled after the merging process 
of the overlapping markers 
}
}
\value{
\item{dat}{relabeled module data of \code{grp}}
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
job.msea <- ssea.start.configure(job.msea)

## Import moddata:
moddata <- tool.read(job.msea$modfile, c("MODULE", "GENE"))
moddata <- unique(na.omit(moddata))

## Import marker (e.g. locus) values:
locdata <- tool.read(job.msea$locfile, c("LOCUS", "VALUE"))
locdata$VALUE <- as.double(locdata$VALUE)
rows <- which(0*(locdata$VALUE) == 0)
locdata <- unique(na.omit(locdata[rows,]))
locdata_ex <- locdata
names(locdata_ex) <- c("MARKER","VALUE")
## Import mapping data between genes and markers: 
gendata <- tool.read(job.msea$genfile, c("GENE", "LOCUS"))
gendata <- unique(na.omit(gendata))
gendata_ex <- gendata
names(gendata_ex) <- c("GENE","MARKER")

## Remove genes with no marker values:
pos <- match(gendata$LOCUS, locdata$LOCUS)
gendata <- gendata[which(pos > 0),]

## Merge overlapping genes:
gendata <- tool.coalesce(items=gendata$LOCUS, groups=gendata$GENE,
rcutoff=job.msea$maxoverlap)
job.msea$geneclusters <- gendata[,c("CLUSTER","GROUPS")]
job.msea$geneclusters <- unique(job.msea$geneclusters)

## Update gene symbols after merging the overlapping ones:
moddata <- ssea.start.relabel(moddata, gendata)  
gendata <- unique(gendata[,c("GROUPS", "ITEM")])
names(gendata) <- c("GENE", "LOCUS")

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
\code{\link{ssea.start}}
}
