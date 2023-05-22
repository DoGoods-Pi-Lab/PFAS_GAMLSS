##PFAS script analysis parts 4-6

#Author: Rebekah Petroff
#Date of Creation: 3/31/2022


###Introduction
#The goal of this script is to assess the link between PFAS exposure in first 
#trimester pregnancy and DNA methylation/hydroxymethylation at birth. To do so, 
#we are using preprocessed beta matrices that describe the amount of methylation
#and hydroxymethylation. We also have a phenotype file with all relevant covariates
#already described. This script is designed to complete the following actions. 

#1. Setup - load files, R packages, and define the main sample and subset sample.
#2. Randomly select 5000 CpGs to fit models (making sure that they are CpGs in OX-BS dataset)
#3. Determine best and second best fit models on these CpGs across all PFAS
#Model selection will use GAMLSS 
#Y distributions on beta and normal distributions
#add/drop terms on smoking, sex, and cell type to see which ones produce best fit
#calculate average AIC and genomic lambdas across each model
#4. Fit best and second best model to all BS methylation data (n=141)
#Extract CpGs meeting the  top 500 + any q<0.2
#5. Fit models with interaction term with OX and BS data (n=70), allowing for variance to be different between groups
#6. Fit models for ONLY OX and BS on select sites meeting q-value <0.2. 

#Based on these results, can consider sex-interaction in the future. 




###########################################################################################################################
#1. Setup

#####1a. Load libraries and files
library(tidyverse)
library(gamlss)
library(QCEWAS)


setwd(getwd())

#most of the setup was already done, so no need to rehash things. 
#we are selecting two models: a restricted term one and an exposure only one. 
load("betaAll.Rdata")
load("annotAll.Rdata")
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
    
# Optional: Test on a small number of probes
betaAllt<- betaAllt[,1:1000]
annotAll<- annotAll[1:1000,]


GL <- matrix(nrow=length(expAll), ncol=2)
rownames(GL) <- names(expAll)
colnames(GL) <- c("Beta Regression: No Sex-Interxn", "Beta Regression: Sex-Interxn")
fit1Names <- c("intercept", "exp", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC")
fit2Names <- c("intercept", "exp", "PC1", "PC2", "Parity", "race_B", "smoke_Preg", "sex", "CD4T", "CD8T", "Gran", "nRBC", "Sex-Interaction")

####run the model!
#for(i in 1:ncol(expAll)){

i =1 
  #run model for the ith exposure
  expn <- expAll[,i]
  tmp <- mclapply(betaAllt, function(x) fit.function1(x, expn),mc.cores = 10, mc.preschedule = FALSE, mc.silent=TRUE) 

  
  # Create data frames of modeling results and add variable name:
  tmp<- plyr::ldply(tmp, data.frame) 
  names(tmp) <- c("cpg", "estimate", "se", "tval", "pval")
  
  tmp$terms <- rep(fit1Names, ncol(betaAllt))
  ### Multiple testing correction using false discovery rate
  # Pull out p-values from modeling results data frame:
  # Adjust p-values using Benjamini-Hochberg FDR adjustment method in p.adjust function
  # Append the FDR values to the dataset.
 
  tmpexp <- tmp %>% subset(terms=="exp")
  names(tmpexp) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  tmpexp$fdr <- p.adjust(tmpexp$pval, method = "BH", n = nrow(tmpexp))
  tmpexp <- cbind(tmpexp, annotAll) #annotate
  
  #caluclate genomic lambda
  GL[i,1] <- P_lambda(tmpexp$pval)
  
  FDR <- tmpexp[(tmpexp$fdr < 0.2), ]
  top500 <- tmpexp%>%slice_min(pval, n=500)
  
  ##assign it a name and save it
  myname <-  paste0("all_", names(expAll[i]))
  assign(myname, tmpexp)
  write.csv(tmpexp, paste0("./output/BS/all/",names(expAll[i]),"_results.csv"))
  
  write.csv(FDR, paste0("./output/BS/FDR_0.2/",names(expAll[i]),"_results.csv"))
  write.csv(top500, paste0("./output/BS/top500/",names(expAll[i]),"_results.csv"))
  
  print(paste0("Finished main model for", names(expAll[i]),".")) 

  ##final go
  tmpI <-mclapply(betaAllt, function(x) fit.function2(x, expn), mc.preschedule = FALSE, mc.cores = 10, mc.silent=TRUE) 
  tmpI<- plyr::ldply(tmpI, data.frame)
  
  tmpI$terms <- rep(fit2Names, ncol(betaAllt))
  ### Multiple testing correction using false discovery rate
  # Pull out p-values from modeling results data frame:
  # Adjust p-values using Benjamini-Hochberg FDR adjustment method in p.adjust function
  # Append the FDR values to the dataset.
  
  #caluclate genomic lambda
  GL[i,2] <- P_lambda(tmpexp$pval)
    
  #look at sex interaction
  tmpSI <- tmpI %>% subset(terms=="Sex-Interaction")
  names(tmpSI) <- c("cpg", "estimate", "se", "tval", "pval", "coef type")
  tmpSI$fdr <- p.adjust(tmpSI$pval, method = "BH", n = nrow(tmpexp))
  tmpSI <- cbind(tmpSI, annotAll) #annotate 
  myname <-  paste0("allComp", names(expAll[i]))
  assign(myname, tmpSI)
  write.csv(tmpSI, paste0("./output/BS/Sex-InteractionCoef/",names(expAll[i]),"_SEXinterxn_results.csv"))

  print(paste0("Done modeling ", names(expAll[i]),"."))

  remove(tmp, tmpexp, FDR, top500, tmpI, tmpSI)
#}
      
write.csv(GL, "./output/BS/genomic_lambdas.csv")

print(paste0("Finished test data fiting: BS only, full model & sex-interaction model."))


###########################################################################################################################

#save.image(file='session04122022.RData')
#end script - use files to visualize results and determine if you need to run sex-stratified models


