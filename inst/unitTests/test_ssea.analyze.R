test_ssea.analyze <- function() {
    library(RUnit)
    library(Mergeomics)
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
    set.seed(1)
    ## ssea.start() process takes long time while merging the genes 
    ## sharing high amounts of markers (e.g. loci). it is performed 
    ## with full module list in the vignettes. Here, we used a very 
    ## subset of the module list (1st 10 mods from the original module file)
    ## and we collected the corresponding genes and markers belonging to 
    ## these modules:
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
    ## run ssea.start() and prepare for this small set: 
    ## (due to the huge runtime)
    job.msea <- ssea.start(job.msea)
    job.msea <- ssea.prepare(job.msea)
    job.msea <- ssea.control(job.msea)
    job.msea <- ssea.analyze(job.msea)
    
    ## compare the pvals with the expected ones:
    ## since we set the seed for random # generation, we know the exact
    ## results for our input sets:
    checkEqualsNumeric(sort(job.msea$results$P)[1], 1.35e-20,
    tolerance=1.0e-4)

    
    ## Remove the temporary files used for the test:
    file.remove("subsetof.coexpr.modules.txt")
    file.remove("subsetof.genfile.txt")
    file.remove("subsetof.marfile.txt")
}
