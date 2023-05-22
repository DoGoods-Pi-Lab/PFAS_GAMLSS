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


