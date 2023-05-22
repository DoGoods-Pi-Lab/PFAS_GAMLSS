##PFAS script analysis part 3b

#Author: Rebekah Petroff
#Date of Creation: 4/14/2022


###Introduction
#The goal of this script is to assess the link between PFAS exposure in first 
#trimester pregnancy and DNA methylation/hydroxymethylation at birth. To do so, 
#we are using preprocessed beta matrices that describe the amount of methylation
#and hydroxymethylation. We also have a phenotype file with all relevant covariates
#already described. This script is designed to complete the following actions. 

#Setup - load files, R packages, and define the main sample and subset sample. Randomly select 5000 CpGs to fit models (making sure that they are CpGs in OX-BS dataset)
#1. Determine best and second best fit models on these CpGs across all PFAS
#Model selection will use GAMLSS 
#Y distributions on beta and normal distributions
#add/drop terms on smoking, sex, and cell type to see which ones produce best fit
#calculate average AIC and genomic lambdas across each model
#2. Fit best and second best model to all BS methylation data (n=141)
#Extract CpGs meeting the  top 500 + any q<0.2
#3. Fit models with interaction term with OX and BS data (n=70), allowing for variance to be different between groups
#4. Fit models for ONLY OX and BS on select sites meeting q-value <0.2. 

#Based on these results, can consider sex-interaction in the future. 




###########################################################################################################################
#####Load libraries and files
library(tidyverse)
library(QCEWAS)

setwd(getwd())

#models have been run in 53 parts to speed computing 
load("../Part2/expAll.Rdata")

annot <- read.csv("../annotation.csv")


GL <- matrix(nrow=length(expAll),ncol=1)
rownames(GL) <- names(expAll)
colnames(GL) <- c("Beta Regression: No Sex-Interxn")

####collect the results from the array

for(i in 1:ncol(expAll)){
  #run model for the ith exposure
  #read results from 53 runs
    
   maindat <- matrix(nrow=0, ncol= 6)

	
    for(k in seq(1:53)){
  	tmp <- read.csv(paste0("./array_job/output/",names(expAll[i]),"_", k,"_results.csv"))
	maindat <- rbind(maindat, tmp)
	    }

print(paste0("Done pulling all files for ", names(expAll[i]),"."))

  ### Multiple testing correction using false discovery rate
  # Pull out p-values from modeling results data frame:
  # Adjust p-values using Benjamini-Hochberg FDR adjustment method in p.adjust function
  # Append the FDR values to the dataset.

  maindat <- unique(maindat)
  maindat$fdr <- p.adjust(maindat$pval, method = "BH", n = nrow(maindat))
  maindat <- left_join(maindat, annot, by="cpg") #annotate
  
  #caluclate genomic lambda
  GL[i,1] <- P_lambda(maindat$pval)
  
  FDR <- maindat[(maindat$fdr < 0.2), ]
  
  ##assign it a name and save it
  write.csv(maindat, paste0("./FinalAll/",names(expAll[i]),"_all.csv"), row.names=FALSE)
  
  write.csv(FDR, paste0("./FDR_0.2/",names(expAll[i]),"_fdr.csv"), row.names=FALSE)
  
print(paste0("Done writing main model for ", names(expAll[i]),"."))
}
      
write.csv(GL, "./genomic_lambdas.csv")

print(paste0("Done with all"))
###########################################################################################################################

