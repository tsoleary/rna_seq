# R script to compare several conditions with the SARTools and DESeq2 packages 
# designed to be executed with SARTools 1.6.9
# Hugo Varet
# March 20th, 2018
rm(list=ls())

# parameters: to be modified by the user ---------------------------------------
workDir <- here::here("cahan")

projectName <- "whole_body_heat_cold_shock" 
author <- "TSO"                              

# path to metadata file
targetFile <- "samples.txt"
# path to the directory containing raw counts files
rawDir <- "counts" 

# names of the features to be removed
# (specific HTSeq-count information and rRNA for example)
# NULL if no feature to remove
featuresToRemove <- c("alignment_not_unique",        
                      "ambiguous", "no_feature",     
                      "not_aligned", "too_low_aQual")

# factor of interest
varInt <- "treat"                                    
# reference biological condition
condRef <- "ctrl"                                    
# blocking factor: NULL (default) or "batch" for example
batch <- NULL                                        


# mean-variance relationship: "parametric" (default), "local" or "mean"
fitType <- "parametric"                              
# TRUE/FALSE to perform the outliers detection (default is TRUE)
cooksCutoff <- TRUE                                  
# TRUE/FALSE to perform independent filtering (default is TRUE)
independentFiltering <- TRUE                        
# threshold of statistical significance
alpha <- 0.05                                        
# p-value adjustment method: "BH" (default) or "BY"
pAdjustMethod <- "BH"                                

# transformation for PCA/clustering: "VST" or "rlog"
typeTrans <- "VST"                                   
# "median" (default) or "shorth" to estimate the size factors
locfunc <- "median"                                  
# vector of colors of each biological condition on the plots
colors <- c("springgreen", "dodgerblue", "firebrick1")             

forceCairoGraph <- FALSE

# run script -------------------------------------------------------------------
setwd(workDir)
require(SARTools)
if (forceCairoGraph) options(bitmapType = "cairo")

# checking parameters
checkParameters.DESeq2(projectName = projectName,
                       author = author,
                       targetFile = targetFile,
                       rawDir = rawDir,
                       featuresToRemove = featuresToRemove,
                       varInt = varInt,
                       condRef = condRef,
                       batch = batch,
                       fitType = fitType,
                       cooksCutoff = cooksCutoff,
                       independentFiltering = independentFiltering,
                       alpha = alpha,
                       pAdjustMethod = pAdjustMethod,
                       typeTrans = typeTrans,
                       locfunc = locfunc,
                       colors = colors)

# loading target file
target <- loadTargetFile(targetFile = targetFile, 
                         varInt = varInt, 
                         condRef = condRef, 
                         batch = batch)

# loading counts
counts <- loadCountData(target = target, 
                        rawDir = rawDir, 
                        featuresToRemove = featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts = counts, 
                                 group = target[, varInt], 
                                 col = colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts = counts, 
                         target = target, 
                         varInt = varInt, 
                         batch = batch,
                         locfunc = locfunc, 
                         fitType = fitType, 
                         pAdjustMethod = pAdjustMethod,
                         cooksCutoff = cooksCutoff, 
                         independentFiltering = independentFiltering, 
                         alpha = alpha)

# PCA + clustering
exploreCounts(object = out.DESeq2$dds, 
              group = target[ , varInt], 
              typeTrans = typeTrans, 
              col = colors)

# summary of the analysis (boxplots, dispersions, diag size factors, 
# export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, 
                                          group = target[ , varInt], 
                                          col = colors,
                                          independentFiltering = 
                                            independentFiltering,
                                          cooksCutoff = cooksCutoff, 
                                          alpha = alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target = target, 
                   counts = counts, 
                   out.DESeq2 = out.DESeq2, 
                   summaryResults = summaryResults,
                   majSequences = majSequences, 
                   workDir = workDir, 
                   projectName = projectName, 
                   author = author,
                   targetFile = targetFile, 
                   rawDir = rawDir, 
                   featuresToRemove = featuresToRemove, 
                   varInt = varInt,
                   condRef = condRef, 
                   batch = batch, 
                   fitType = fitType, 
                   cooksCutoff = cooksCutoff,
                   independentFiltering = independentFiltering, 
                   alpha = alpha, 
                   pAdjustMethod = pAdjustMethod,
                   typeTrans = typeTrans, 
                   locfunc = locfunc, 
                   colors = colors)

