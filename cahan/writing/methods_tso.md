# Frontiers Paper

# Methods

## Differential Expression

The raw reads were mapped to DM6...

RNAseq count data of heat shock and cold shock versus the control were analyzed using the DESeq2 package [version X ;@Love2014]. The default normalization was used. Differential gene expression was called at an FDR of 1% (padj < 0.01) with no log fold-change cut off.  

## Integration with the GWAS

The top SNPs (AvgMixedPval < $10^{-4}$)​ in the GWAS were annotated to Flybase 5.57? Each associated gene was matched with the corresponding gene in the expression data set

# Results

## Differential Gene Expression

Almost a third of the genes in the _Drosophila melanogaster_ genome were down-regulated compared to control in each of the temperature shock conditions (5,126 in the cold shock; 6,241 in the heat shock). A smaller set of genes were up-regulated in each condition (1826 in the cold; 2314 in the heat). The overall transcriptional response was largely shared between the two extreme temperature treatments. A majority of genes fall along the 1:1 line (Figure X; 6000 shared) with only a small subset of genes with a unique response in different directions depending on the temperature treatment.

## Integration with the GWAS

 59 of the 99 unique genes assocated with $CT_{max}$, belong to differentially expressed genes under heat-shock. And of the 151 unique gene associated with $CT_{min}$, 72 of them are differentially expressed in the cold-shock condition. 



