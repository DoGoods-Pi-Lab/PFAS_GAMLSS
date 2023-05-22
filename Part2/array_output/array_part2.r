##PFAS script analysis part 2a

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
step=12000
begin <- array_index*step-(step-1)
end <- array_index*step


setwd(getwd())

#most of the setup was already done, so no need to rehash things. 
#we are selecting two models: a restricted term one and an exposure only one. 
load("betaAll.Rdata")
load("phenoAll.Rdata")
load("expAll.Rdata")

print(paste0("Starting model fitting process. After filtering, ", ncol(betaAllt)," probes are included for bisulfite convereted data over ",nrow(phenoSelect)," individuals are modeled."))



##write the function to fit the model      
fit.function1 <- function(mycpg, myexp){
  #my model fit and design
  simple.fit <- gamlss(data=phenoSelect, mycpg ~  myexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                         sex + CD4T + CD8T + Gran + nRBC, family = BE, trace = FALSE)
  # Extract the estimate, se, t-value and p-value for beta coefficient (interaction term)
  #Note: change the [6,c(1:4)], depending on how many terms are in the model.
  simple.sum <- summary(simple.fit, save=TRUE, trace=FALSE)
  simple.coef <- simple.sum$coef.table[-13,] #take all but leave the coef for sigma var 

  return(simple.coef)
}


##write the function to fit the model      
fit.function2 <- function(mycpg, myexp){
  #fit complex beta distributions
  complex.fit <- gamlss(data=phenoSelect, mycpg ~  myexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                                        sex + CD4T + CD8T + Gran + nRBC + sex*myexp, family = BE, trace = FALSE)
  complex.sum <- summary(complex.fit, save=TRUE, trace=FALSE)
  complex.coef <- complex.sum$coef.table[-14,] #take all but leave the coef for sigma var 
  return(complex.coef)
}
    

fit1Names <- c("intercept", "exp", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC")
fit2Names <- c("intercept", "exp", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC", "Sex-Interaction")

####run the model!
#for(i in 1:ncol(expAll)){

#i =1 
i <- as.numeric(Sys.getenv("exposure"))


  # Optional: Test on a small number of probes

  betaAllt<- betaAllt[,begin:end]

  #run model for the ith exposure
  expn <- expAll[,i]
  tmp <- lapply(betaAllt, function(x) fit.function1(x, expn)) 

  
  # Create data frames of modeling results and add variable name:
  tmp<- plyr::ldply(tmp, data.frame) 
  names(tmp) <- c("cpg", "estimate", "se", "tval", "pval")
  
  tmp$terms <- rep(fit1Names, ncol(betaAllt)) 
  tmpexp <- tmp %>% subset(terms=="exp")
  names(tmpexp) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  
  ##assign it a name and save it
  write.csv(tmpexp, paste0("./output/BS/all/",names(expAll[i]),"_", array_index, "_results.csv"))
    
  print(paste0("Finished main model for", names(expAll[i]),".")) 

  ##final go
  tmpI <-lapply(betaAllt, function(x) fit.function2(x, expn)) 
  tmpI<- plyr::ldply(tmpI, data.frame)
  
  tmpI$terms <- rep(fit2Names, ncol(betaAllt))  #look at sex interaction
  tmpSI <- tmpI %>% subset(terms=="Sex-Interaction")
  names(tmpSI) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")

  write.csv(tmpSI, paste0("./output/BS/Sex-InteractionCoef/",names(expAll[i]),"_", array_index, "_SEXinterxn_results.csv"))

  print(paste0("Done modeling ", names(expAll[i]),"."))

  remove(tmp, tmpexp, tmpI, tmpSI)
#}
print(paste0("Finished test data fiting: BS only, full model & sex-interaction model."))


###########################################################################################################################


#end script 


