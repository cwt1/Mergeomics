\name{kda.configure}
\alias{kda.configure}
\title{
Set parameters for weighted key driver analysis (wKDA) 
}
\description{
takes the configuration (plan) parameter for wKDA process as input and 
assigns default values if needed.
The fields of this parameter are listed in the arguements section in 
detail.
}
\usage{
kda.configure(plan)
}
\arguments{
\item{plan}{
a parameter including fields about the details of the wKDA process: 
\preformatted{
label:    unique identifier for the analysis
folder:   parent folder for results
netfile:  path to network file (TAIL HEAD WEIGHT)
modfile:  path to module file (MODULE GENE)
inffile:  path to module info file
nodfile:  path to node selection file
depthsearch: depth for subgraph search
direction: 0 for undirected, negative for downstream and
positive for upstream
maxoverlap: maximum allowed overlap between two key driver
neighborhoods
minsize:  minimum module size
mindegreeminimum: node degree to qualify as a hub
maxdegreemaximum: node degree to include
edgefactor: influence of node strengths: 0.0 no influence,
1.0 full influence
seed: seed for random number generator
}
}
}
\details{
\code{\link{kda.configure}} prepares the environment for wKDA process, 
checks the fields of the input plan parameter (that includes paths of 
required input files and output folder, min module size, etc.), and 
assigns the default values to these fields if they are not specified. 
}
\value{
\item{plan}{configured and -if needed updated- plan parameter to be used
in wKDA process. }
}
\examples{
## for KDA the essential parameters should be assigned by user is as follows:
plan <- list()
## assign job label:
plan$label<-"HDLC"
## specify parent folder for results:
plan$folder<-"Results"
## Get an input network (columns: TAIL HEAD WEIGHT)
plan$netfile <-"network.mouseliver.mouse.txt"
## Get the gene sets derived from ModuleMerge, containing two columns, 
## MODULE and NODE, delimited by tab 
plan$modfile<- "moddata.txt"
## If above parameters are not assigned by users, code will stop with error:
if(is.null(plan$folder)) stop("No parent folder.")
if(is.null(plan$label)) stop("No job label.")
if(is.null(plan$netfile)) stop("No network file.")
if(is.null(plan$modfile)) stop("No module file.")

## other parameters are optional, if they are not specified by user,
## kda.configure assigns their default values:
## graph search depth parameter:
if(is.null(plan$depth)) plan$depth <- 1
## edge directionality in the network: O means undirected
if(is.null(plan$direction)) plan$direction <- 0
## max overlap allowed between two modules
if(is.null(plan$maxoverlap)) plan$maxoverlap <- 0.33
## min size of the modules 
if(is.null(plan$minsize)) plan$minsize <- 20
## min and max hub degree to be included:
if(is.null(plan$mindegree)) plan$mindegree <- "automatic"
if(is.null(plan$maxdegree)) plan$maxdegree <- "automatic"
## number of simulations for permutation test:
if(is.null(plan$nperm)) plan$nperm <- 2000
## seed for random number generator:
if(is.null(plan$seed)) plan$seed <- 1
## these are the main parameters needed to be assigned default values.
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
\code{\link{kda.analyze}}
}

