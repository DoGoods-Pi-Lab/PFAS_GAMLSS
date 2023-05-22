##PFAS script analysis part 3a

#Author: Rebekah Petroff
#Date of Creation: 3/31/2022


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
#1. Setup

#####1a. Load libraries and files
library(tidyverse)
library(gamlss)
library(QCEWAS)

array_index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
step=10000
begin <- array_index*step-(step-1)
end <- array_index*step


setwd(getwd())

#most of the setup was already done, so no need to rehash things. 
#we are selecting two models: a restricted term one and an exposure only one. 
load("../session04192022.Rdata")

#for these data, n=70 subjects with ~520,000 probes

#phenoSelectOXBS <- phenoOXBS %>% dplyr::select(c(PC1, PC2, Parity, race_B, smoke_Preg, sex, CD4T, CD8T, Gran, nRBC, Treatment))


##write the function to fit the model      
#fit.function <- function(mycpg, myexp){
  #my model fit and design
  #fit <- gamlss(data=phenoSelectOXBS, mycpg ~  myexp + Treatment + myexp*Treatment + PC1 + PC2 + Parity + race_B + smoke_Preg +
                         #sex + CD4T + CD8T + Gran + nRBC, random= ~1 | randeff, 
			 #sigma.formula = ~ Treatment,
			  #family = BE, trace = FALSE)
  # Extract the estimate, se, t-value and p-value for beta coefficient (interaction term)
  #Note: change the [6,c(1:4)], depending on how many terms are in the model.
  #sum <- summary(fit, save=TRUE, trace=FALSE)
  #coef <- sum$coef.table[1:14,] #take all but leave the coef for sigma var 

  #return(coef)
#}


####2. run the model!
  i <- as.numeric(Sys.getenv("exposure"))

  betaOXBSt<- betaOXBSt[,begin:end]
  #for testing: betaOXBSt<- betaOXBSt[,1:250]
  
  #run model for the ith exposure
  expn <- expOXBS[,i]
  tmp <- lapply(betaOXBSt, function(x) fit.function(x, expn)) 

  
  # Create data frames of modeling results and add variable name:
  tmp<- plyr::ldply(tmp, data.frame) 
  names(tmp) <- c("cpg", "estimate", "se", "tval", "pval")
  
  fitNames <- c("intercept", "exp", "OX-BS", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC", "expOXBS_interxn")
  tmp$terms <- rep(fitNames, ncol(betaOXBSt)) 
  intrxn <- tmp%>% subset(terms=="expOXBS_interxn")
  names(intrxn) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  
  ##assign it a name and save it
  write.csv(intrxn, paste0("./output/",names(expOXBS[i]),"_", array_index, "_results.csv"))
    


###########################################################################################################################

#end script - finish the concat of array job files and look for sig interactions


