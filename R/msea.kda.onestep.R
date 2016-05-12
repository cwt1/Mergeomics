MSEA.KDA.onestep <- function(plan, apply.MSEA=TRUE, apply.KDA=FALSE, 
maxoverlap.genesets=0.33, symbol.transfer.needed=FALSE, 
sym.from=c("HUMAN", "MOUSE"), sym.to=c("HUMAN", "MOUSE")){

    if((!apply.MSEA) & (!apply.KDA)){
        stop("You need to set either apply.MSEA argument or apply.KDA argument
        to TRUE.")
    }

    cat("\tIf you want to perform only MSEA, set apply.MSEA=T and
    apply.KDA=F.\n")
    cat("\tIf you want to perform only KDA, set apply.MSEA=F and 
    apply.KDA=T.\n")

    if(apply.MSEA)
        plan.check <- ssea.start.configure(plan)   
    if(apply.KDA)
        plan.check <- kda.configure(plan)
   
    ################ MSEA (Marker set enrichment analysis)  ###
    if(apply.MSEA){
        plan <- ssea.start(plan)
        plan <- ssea.prepare(plan)
        plan <- ssea.control(plan)
        plan <- ssea.analyze(plan)
        plan <- ssea.finish(plan)
    }

    ## if you want to first run MSEA, then run KDA:
    if((apply.MSEA) & (apply.KDA)){
        netfile.backup <- plan$netfile
        
        ## delete control groups if there are any:
        indA <- which(plan$modules == "_ctrlA")
        if (length(indA) > 0)  
            indA <- which(plan$results$MODULE == indA)
        if (length(indA) > 0) plan$results <- plan$results[-indA,]
        indB <- which(plan$modules == "_ctrlB")
        if (length(indB) > 0)  
            indB <- which(plan$results$MODULE == indB)
        if (length(indB) > 0) plan$results <- plan$results[-indB,]
        
        #########  Create intermediary datasets for KDA ##########
        if(symbol.transfer.needed){
            syms <- tool.read(system.file("extdata", "symbols.txt", 
            package="Mergeomics"))
            syms <- syms[,c(as.character(sym.from), as.character(sym.to))]
            names(syms) <- c("FROM", "TO")
            plan <- ssea2kda(plan, symbols=syms, rmax=maxoverlap.genesets)
        }
        else
            plan <- ssea2kda(plan, rmax=maxoverlap.genesets)
        plan$netfile <- netfile.backup 
    }

    ## if you want to both directly run KDA or right after running the MSEA:
    if (apply.KDA){
        #######   wKDA (Weighted key driver analysis)    ##########
        plan <- kda.configure(plan)
        plan <- kda.start(plan)
        plan <- kda.prepare(plan)
        plan <- kda.analyze(plan)
        plan <- kda.finish(plan)
        ######  Prepare network files for visualization   #########
        ## Creates the input files for Cytoscape (http://www.cytoscape.org/)
        plan <- kda2cytoscape(plan)
    }
    return(plan)
}