Right now this is just an incomplete rough draft of the methods that I was involved in.

# Methods

## Differential Expression

The raw reads were mapped to the _Drosophila melaogaster_ reference genome (Release 6) using STAR (info from Seth ?). Those aligned reads were mapped to genes using featureCounts (info from Seth about version?).

The RNAseq count data of heat shock and cold shock were compared to the control were analyzed using the DESeq2 package (version 1.24.0; @Love2014) in R (version 3.6.1) using default parameters. Differential gene expression was called at an FDR of 1% (Benjamini-Hochberg corrected p-value < 0.01). 

## Integration with the GWAS

The top SNPs (AvgMixedPval < $10^{-4}$)​ in the GWAS were annotated to Flybase 5.57. Each associated gene was matched with the corresponding gene in the expression data set

# Results

## Sequencing

Do we need anything talking about the sequencing results?

## Differential Gene Expression

Almost a third of the genes in the _Drosophila melanogaster_ genome were down-regulated compared to control in each of the temperature shock conditions (5,126 in the cold shock; 6,241 in the heat shock). A smaller set of genes were up-regulated in each condition (1826 in the cold; 2314 in the heat). The overall transcriptional response was largely shared between the two extreme temperature treatments. A majority of genes fall along the 1:1 line (Figure X; 6000 shared) with only a small subset of genes with a unique response in different directions depending on the temperature treatment.

## Integration with the GWAS

59 of the 99 unique genes assocated with $CT_{max}$, belong to differentially expressed genes under heat-shock. And of the 151 unique gene associated with $CT_{min}$, 72 of them are differentially expressed in the cold-shock condition. 



