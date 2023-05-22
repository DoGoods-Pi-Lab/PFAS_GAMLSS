##PFAS script analysis part 4

#Author: Rebekah Petroff
#Date of Creation: 4/28/2022


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

###########################################################################################################################
#this is ONLY modeling the indepented OX and BS data with those sites that had a significant interaction term. 
#####Load libraries and files
library(tidyverse)
library(QCEWAS)
library(gamlss)

setwd(getwd())
load("./session04282022.Rdata")

annot <- read.csv("../annotation.csv")



fit.function <- function(mycpg){
#my model fit and design
fit <- gamlss(data=pheno, mycpg ~  expn +  PC1 + PC2 + Parity + race_B + smoke_Preg + sex + CD4T + CD8T + Gran + nRBC, family = BE, trace = FALSE)
# Extract the estimate, se, t-value and p-value for beta coefficient (interaction term)
#Note: change the [6,c(1:4)], depending on how many terms are in the model.
sum <- summary(fit, save=TRUE, trace=FALSE)
coef <- sum$coef.table[1:12,] #take all but leave the coef for sigma var 

return(coef)
}

fitNames <- c("intercept", "exp", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC")

pheno <- pheno %>% dplyr::select(-Treatment)

print(paste("Starting Modeling..."))

###need to run by hand b/c only some things have more than 1 CpG to look at.

for(i in 1:ncol(expAll)){
  interxnCpGs <- read.csv(paste0("../Part3/FDR_0.2/",names(expAll[i]),"_fdr.csv"))  

  ##select HMC subset
  pheno$expn <- as.numeric(expHMC[,i])
  pheno <- as.data.frame(apply(pheno,2,as.numeric))		
	
  hmc_tmp <- hmc[,which(colnames(hmc) %in% interxnCpGs$cpg)]
  hmc_tmp <- as.data.frame(apply(hmc_tmp, 2, as.numeric))
  
  #fit the model
  tmpout <- lapply(hmc_tmp, function(x) fit.function(x)) 
  
  # Create data frames of modeling results and add variable name:
  tmpout<- plyr::ldply(tmpout, data.frame) 
  names(tmpout) <- c("cpg", "estimate", "se", "tval", "pval")
  

  tmpout$terms <- rep(fitNames, ncol(hmc_tmp)) 
  tmpout <- tmpout%>% subset(terms=="exp")
  names(tmpout) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  tmpout$fdr <- p.adjust(tmpout$pval, method = "BH", n = nrow(tmpout))
  tmpout <- left_join(tmpout, annot, by= "cpg") #annotate
  
  ##assign it a name and save it
  write.csv(tmpout, paste0("./output/HMC/",names(expHMC[i]),"_results.csv"))
  remove(tmpout, hmc_tmp)
  
  print(paste0("Done modeling hmc for ", names(expHMC[i]),"."))
  
  ##############################MC ONLY

  mc_tmp <- mc[,which(colnames(mc) %in% interxnCpGs$cpg)]
  mc_tmp <- as.data.frame(apply(mc_tmp, 2, as.numeric))
  

  #fit the model
  tmpout <- lapply(mc_tmp, function(x) fit.function(x)) 
  
  # Create data frames of modeling results and add variable name:
  tmpout<- plyr::ldply(tmpout, data.frame) 
  names(tmpout) <- c("cpg", "estimate", "se", "tval", "pval")
  
  
  tmpout$terms <- rep(fitNames, ncol(mc_tmp)) 
  tmpout <- tmpout%>% subset(terms=="exp")
  names(tmpout) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  tmpout$fdr <- p.adjust(tmpout$pval, method = "BH", n = nrow(tmpout))
  tmpout <- left_join(tmpout, annot, by= "cpg") #annotate
  
  ##assign it a name and save it
  write.csv(tmpout, paste0("./output/MC/",names(expMC[i]),"_results.csv"))
  
  print(paste0("Done modeling mc for ", names(expMC[i]),"."))
  remove(tmpout, mc_tmp)
}