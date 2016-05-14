\name{Mergeomics-package}
\alias{Mergeomics}
\alias{Mergeomics-package}
\docType{package}
\title{
Integrative network analysis of omics data
}
\description{
The Mergeomics pipeline serves as a flexible framework for integrating 
multidimensional omics-disease associations, functional genomics, 
canonical pathways and gene-gene interaction networks to generate 
mechanistic hypotheses. It includes two main parts,
1) Marker set enrichment analysis (MSEA);
2) Weighted Key Driver Analysis (wKDA).

}
\details{
\tabular{ll}{
Package: \tab Mergeomics\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2016-01-04\cr
License: \tab GPL (>= 2)\cr
Depends: R (>= 3.0.1)\cr
URL: \tab http://mergeomics.research.idre.ucla.edu/\cr
}

Mergeomics amalgamates disease association information derived from 
multidimensional omics data (e.g., genome, epigenome, transcriptome, 
metablome) with functional genomics (e.g., eQTLs, ENCODE), canonical 
pathways (e.g., KEGG, Reactome), and molecular networks (e.g., gene 
regulatory networks, protein-protein interaction networks). 
Two main steps of the pipeline are: Marker set enrichment analysis 
(MSEA) and weighted key driver analysis (wKDA). MSEA takes the following 
data as input: i) disease association data (GWAS, EWAS, TWAS...), 
ii) functional genomics (eQTLs and/or ENCODE information), and iii) 
functionally related genes information extracted from knowledge-based 
biological pathways or data-driven network modules (e.g., coexpressed 
genes in a given tissue relevant to a disease of interest). These datasets 
are integrated via MSEA to return gene sets that are significantly enriched 
for markers showing low p value associations with a given disease. Then, 
the disease related gene sets are examined to detect the key drivers by 
using the wKDA step of the pipeline, which requires pre-defined directional 
networks such as tissue-specific Bayesian networks, protein-protein 
interaction networks, etc. wKDA maps the disease related gene sets to the 
pre-defined directional networks to identify key driver genes that are more 
likely regulators of the disease gene sets based on their central positions 
in the gene networks. The key drivers and their local network topology can 
be viewed and downloaded after the completion of the analysis via 
Visualization step. Our pipeline provides users to perform MSEA and wKDA 
together or separately using either their own input data or selecting 
preloaded sample datasets.
The details of the functions and parameter settings are described in the 
Manual of the package.  
}
\author{
Ville-Petteri Makinen, Le Shu, Yuqi Zhao, Zeyneb Kurt, Bin Zhang, Xia Yang
Maintainer: <zeyneb@ucla.edu>
}
\references{
Shu L, Zhao Y, Kurt Z, Byars S, Tukiainen T, Kettunen J, Ripatti S, Zhang B,
Inouye M, Makinen VP, Yang X. Mergeomics: integration of diverse genomics
resources to identify pathogenic perturbations to biological systems.
bioRxiv doi: http://dx.doi.org/10.1101/036012
}

\keyword{ 
Integrative Genomics; Multidimensional Data Integration; Gene Networks; 
Key Drivers
}
