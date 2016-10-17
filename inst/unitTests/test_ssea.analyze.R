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
    
    ## run ssea.start() 
    job.msea <- ssea.start(job.msea)
    job.msea <- ssea.prepare(job.msea)
    job.msea <- ssea.control(job.msea)
    job.msea <- ssea.analyze(job.msea)
    
    ## compare the pvals with the expected ones:
    ## since we set the seed for random # generation, we know the exact
    ## results for our input sets:
    checkEqualsNumeric(sort(as.numeric(job.msea$results$P))[1], 
    3.67e-33, tolerance=1.0e-4)

}
