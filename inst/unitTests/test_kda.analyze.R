test_kda.analyze <- function() {
    library(RUnit)
    library(Mergeomics)

    job.kda <- list()
    job.kda$label<-"HDLC"
    ## parent folder for results
    job.kda$folder<-"Results"
    ## Input a network
    ## columns: TAIL HEAD WEIGHT
    job.kda$netfile<-system.file("extdata","network.mouseliver.mouse.txt", 
    package="Mergeomics")
    ## Gene sets derived from ModuleMerge, containing two columns, MODULE, 
    ## NODE, delimited by tab 
    job.kda$modfile<- system.file("extdata","mergedModules.txt", 
    package="Mergeomics")
    ## "0" means we do not consider edge weights while 1 is opposite.
    job.kda$edgefactor<-0.0
    ## The searching depth for the KDA
    job.kda$depth<-1
    ## 0 means we do not consider the directions of the regulatory 
    ## interactions; while 1 is opposite.
    job.kda$direction <- 1
    job.kda$nperm <- 100 # the default value is 2000, use 100 for unit tests
    ## Let's run KDA!
    set.seed(1)
    job.kda <- kda.configure(job.kda)
    job.kda <- kda.start(job.kda)
    job.kda <- kda.prepare(job.kda)
    job.kda <- kda.analyze(job.kda)
    
    ## compare the pvals with the expected ones:
    ## since we set the seed for random # generation, we know the exact
    ## results for our input sets:
    ## check FDRs of the top KD
    
    checkEqualsNumeric(sort(as.numeric(job.kda$results$FDR))[1], 
    3.1492e-12, tolerance=1.0e-4)
    
    ## check Pvalues of the top KD
    checkEqualsNumeric(as.numeric(sort(job.kda$results$P))[1], 
    1.2907e-14, tolerance=1.0e-4)
}
