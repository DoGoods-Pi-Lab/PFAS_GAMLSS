##PFAS script analsysis parts 1-3 

#Author: Rebekah Petroff
#Date of Creation: 3/16/2022


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

#Based on these results, can consider sex-intereaction in the future. 

###########################################################################################################################
#1. Setup

#####1a. Load libraries and files
library(tidyverse)
library(gamlss)
library(Hmisc) #for the nin function


#load("//sph-gaydojac-win.turbo.storage.umich.edu/sph-gaydojac/MMIP EPIC/oxBS EPIC for analysis/mlml_mmip74.Rda")
#load("L:/MMIP EPIC/oxBS EPIC for analysis/mlml_mmip74.Rda")
#load("C:/Users/rlpet/Dropbox (University of Michigan)/Petroff - postdoc/Projects/PFAS/analysis/PFAS_bs/preprocess_output/beta_fornosex.rda")
#load("C:/Users/rlpet/Dropbox (University of Michigan)/Petroff - postdoc/Projects/PFAS/analysis/PFAS_bs/preprocess_output/finalpheno.rda")
#load("C:/Users/rlpet/Dropbox (University of Michigan)/Petroff - postdoc/Projects/PFAS/analysis/PFAS_bs/preprocess_output/annot_fornosex.rda")

setwd(getwd())

load("mlml_mmip74.Rda")
load("beta_fornosex.Rda")
load("finalpheno.Rda")
load("annot_fornosex.Rda")

betaAll <- betaXY
phenoAll <- orderedfinalgroup
annotAll <- annotXY


mc <- as.data.frame(mlml_result$mC)
hmc <- as.data.frame(mlml_result$hmC)

mcPheno <- pheno %>% subset(type == "bs")
hmcPheno <- pheno %>% subset(type == "oxbs")

#take out families without PFAS in ox-bs data
mcPheno <- mcPheno[which(mcPheno$FamilyID %in% phenoAll$FamilyID),]
hmcPheno <- hmcPheno[which(hmcPheno$FamilyID %in% phenoAll$FamilyID),]

hmc <- hmc[,colnames(hmc) %in% paste0(mcPheno$biosample2,"_hmC")]
mc <- mc[,colnames(mc) %in% paste0(mcPheno$biosample2,"_mC")]


###1b. Subsetting data into what we will use.
#we can use this drop xys already and any probes that maybe snuck through processing 
hmc <- hmc[rownames(hmc) %in% rownames(betaAll),]
mc <- mc[rownames(mc) %in% rownames(betaAll),]

#drop probes with less than 5% methylation in at least 5% of samples in only the ox-bs combined dataset
minN <- length(hmc)*0.05
lowMethProbes <- row.names(betaAll[rowSums(betaAll <= 0.05) >= minN, ]) 

#number of probes that have less than 5% methylation in at least 5% of people
#length(lowMethProbes)

#drop these probes
hmc <- hmc[!rownames(hmc) %in% lowMethProbes,]
mc <- mc[!rownames(mc) %in% lowMethProbes,]

print(paste0("Finished loading all data."))

######1b - Processing our datasets 

print(paste0("Working on imputing pheno data..."))


#first we need to drop the duplicate sample (randomly) in the ox-bs combined dataset
#then we will impute missing vars in the main dataset and then grab those imputed values for the ox-bs data
#finally we will extract the exposures to allow for easy looping over them

#flip a coin to select the duplication
set.seed(123)


select <- c(TRUE, FALSE)
sample1 <- sample(select, size=1) # this produces a 50/50 chance of selecting either True or False. We can then apply that to our one dup. 
print(sample1)

fam132 <- subset(pheno, FamilyID == 132) #pull out only this family
fam132Discard <- ifelse(sample1 == "FALSE", fam132$biosample2[1], fam132$biosample2[2]) #grab the barcode matching the coin flip

#filter out the sample to not use
mcPheno <- mcPheno %>% subset(biosample2 != paste(fam132Discard))
hmcPheno <- hmcPheno %>% subset(biosample2 != paste(fam132Discard))

#filter out in the dataset as well 
mc <- mc %>% dplyr::select(- paste0(fam132Discard,"_mC"))
hmc <- hmc %>% dplyr::select(- paste0(fam132Discard,"_hmC"))

remove(fam132, sample1, lowMethProbes, minN, fam132Discard, select, pheno, annotXY, betaXY, orderedfinalgroup)


#####1c. Imputing data
phenoAll$race3_imputed <- phenoAll$race_3cat #duplicate instead of writing over

#to imput, we #write our own sample based on the proportions that exist in our population. before using it on our data, we can test on a large imaginary sample population with no data, then replace just our missing values if that looks good
test <- sample(c(0,1,2), size=1000, prob = c((sum(which(phenoAll$race_3cat==0))/length(phenoAll$race_3cat)),
                                             (sum(which(phenoAll$race_3cat==1))/length(phenoAll$race_3cat)), 
                                             (sum(which(phenoAll$race_3cat==2))/length(phenoAll$race_3cat))), 
               replace = TRUE) 
print(test) #looks great!

phenoAll$race3_imputed <- ifelse(phenoAll$race3_imputed == 99, #if value ==99, then: 
                                #replace by a representative coin flip of our other proportion
                                sample(c(0,1,2), size=1, prob =
                                         c((sum(which(phenoAll$race_3cat==0))/length(phenoAll$race_3cat)),
                                           (sum(which(phenoAll$race_3cat==1))/length(phenoAll$race_3cat)), 
                                           (sum(which(phenoAll$race_3cat==2))/length(phenoAll$race_3cat)))), 
                                phenoAll$race3_imputed) #if not 99, keep the same.







sum(which(phenoAll$race3_imputed==99)) #we should have no 99s

#now seperate out the categories for model use.
#1 = race = black
#2 = race = nonwhite, nonblack

phenoAll$race_B <- ifelse(phenoAll$race3_imputed ==1, 1, 0) #category for this var
phenoAll$race_NBNW <- ifelse(phenoAll$race3_imputed ==2, 1, 0) #category for this var

#for smoking_stat: 
#0 = never
#1 =none
#2 = quit before preg
#3 = quit during preg (some smoking in preg)
#4 = smoked during preg
#6 = not reported
phenoAll$smoke_Preg <- ifelse(phenoAll$smoking_stat == 3 | phenoAll$smoking_stat == 4, 1, 0) #if =3 or 4, 1
phenoAll$smoke_unknown <- ifelse(phenoAll$smoking_stat == 6, 1, 0) #if 6, then 1

#checking the work
length(which(phenoAll$smoking_stat==6)) #9 folks have 6 beofre
sum(phenoAll$smoke_unknown) #9 folks check in with that cat

print(phenoAll$smoke_Preg)
sum(phenoAll$smoke_Preg) #7 have 3 or 4 beofre
length(which(phenoAll$smoking_stat==3 | phenoAll$smoking_stat==4)) #and 7 still have the same!


#copy the vars we created to the other pheno s
mcPheno<- left_join(as.data.frame(mcPheno), 
                    phenoAll %>% dplyr::select(FamilyID, race3_imputed, race_NBNW, race_B, smoke_Preg, smoke_unknown), 
                    by= "FamilyID")

hmcPheno<- left_join(as.data.frame(hmcPheno), 
                    phenoAll %>% dplyr::select(FamilyID, race3_imputed, race_NBNW, race_B, smoke_Preg, smoke_unknown), 
                    by= "FamilyID")


#put our pheno together
phenoOXBS <- rbind(mcPheno, hmcPheno)

#create the paired ID for these data only
pairedID <- nrow(phenoOXBS)/2
randeff <- as.factor(c(1:pairedID,1:pairedID)) 


#####1d. Pull out exposure data and finalize pheno variables for each type of dataset

expAll <- as.data.frame(phenoAll) %>% dplyr::select(starts_with("ln")) #some vars we have numerical, log transformed values

#two vars we have categorical detected or not
expAll$MePFOSA_cat <-phenoAll$MePFOSA_cat
expAll$PFUA_cat <- phenoAll$PFUA_cat
rownames(expAll) <-phenoAll$FamilyID


#let's order things to keep them neat later on: 
expList <- c("PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUnDA", "MeFOSAA") #keep for table labels
expoOrder <- c("lnPFHxS", "lnPFOA", "lnPFOS", "lnPFNA", "lnPFDeA", "PFUA_cat", "MePFOSA_cat")


#create the vars for our other datasets as well
expMC <- expAll %>% subset(rownames(expAll) %in% mcPheno$FamilyID)
expHMC <- expMC
expOXBS <- rbind.data.frame(expMC, expHMC)

remove(mlml_result)

#####1d. Prep the beta data for downstream modeling
# Establish probe variable using rownames
mc$probe <- rownames(mc)
hmc$probe <- rownames(hmc)

# Merge 5-mC and 5-hmC data for remaining probes and transpose dataframe for downstream processing
betaOXBS <- inner_join(mc, hmc, by = "probe") 
rownames(betaOXBS) <- betaOXBS$probe 

# Remove extra columns
betaOXBS <- betaOXBS %>% dplyr::select(-probe)

# Transpose for regression
betaAllt <- as.data.frame(t(betaAll))
betaOXBSt <- t(betaOXBS)
remove(betaAll, betaOXBS)


print(paste0("Finished setting up data for model selection"))

###########################################################################################################################
#2. Random CpG Selection & Splitting Test Data

#select random Cpgs (note that we already set our seed earlier in the doc)
namesRandom <- sample(colnames(betaOXBSt), 500, replace = FALSE)

#pull out these cpgs
randomBetaAll <- betaAllt %>% dplyr::select(all_of(namesRandom))

randomBetaAll <- as.data.frame(randomBetaAll)


#define my folds
rand <- sample(nrow(phenoAll))

K1row <- rand[rand %% 5 + 1 == 1]
K2row <- rand[rand %% 5 + 1 == 2]
K3row <- rand[rand %% 5 + 1 == 3]
K4row <- rand[rand %% 5 + 1 == 4]
K5row <- rand[rand %% 5 + 1 == 5]

#checking our folds use everyone
sum(length(K1row) + length(K2row) + length(K3row) + length(K4row) +length(K5row))


#we initially fit our models with K1-4 and then test with K5. To determine best fit,
#we will compare at the average and ranges of AIC and genomic lambdas from the initial models.
#then we will test them against our K5 fold, and make sure they are the best. Top two models that 
#are overall the best fit will get used in the next steps. 

myks <- list(K1row, K2row, K3row, K4row, K5row)
  


phenoSelect <- phenoAll %>% dplyr::select(PC1, PC2, Parity, race_B, race_NBNW, smoke_Preg, smoke_unknown, sex, Bcell, CD4T, CD8T, Gran, Mono, NK, nRBC)

#deciding the model using 5-fold cross validation, looking at 4 different models (a short and long; beta and normal dists)
#uses all PFAS, but only 500 CpGs on the main datasets

print(paste0("Starting model fitting process."))


for (i in 1:length(expAll)){
  singleExp <- as.data.frame(expAll[,i])
  mse1 <- matrix( , nrow=(length(myks)+1), ncol=length(randomBetaAll)) 
  mse2 <- matrix( , nrow=(length(myks)+1), ncol=length(randomBetaAll)) 
  mse3 <- matrix( , nrow=(length(myks)+1), ncol=length(randomBetaAll)) 
  mse4 <- matrix( , nrow=(length(myks)+1), ncol=length(randomBetaAll)) 

  aics <- as.data.frame(c("fit1", "fit2", "fit3", "fit4"))
  names(aics) <- "fit"
  
  for (j in 1:length(myks)){
      #assign vars
      index <- myks[j]
      index <- unlist(index)
      tmpexp <- singleExp[-index,]
      
      tmpbeta <- randomBetaAll[-index,]
      fitdat <- phenoSelect[-index,]

      
      testbeta <- randomBetaAll[index,]
      tmpexp <- singleExp[index,]
      testdat <- phenoSelect[index,]
      testdat$tmpexp <- tmpexp
      
      fitdat$tmpCpG <- NA
      testdat$tmpCpG<- NA
      
      #fit models on folds (minus the 1 to test), looping over each cpg
      for (k in 1:length(randomBetaAll)){
          newCpG <- tmpbeta[,k]
          fitdat$tmpCpG <-newCpG
          
          newCpG <- testbeta[,k]
          testdat$tmpCpG <-newCpG
          
          #fit complex beta distributions
          fit1 <- gamlss(data= fitdat, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + race_NBNW + smoke_Preg + smoke_unknown +
              sex + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, 
              family = BE, trace=FALSE)
      
          #fit simpler beta distributions
          fit2 <- gamlss(data= fitdat, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                sex + CD4T + CD8T + Gran + nRBC, 
                family = BE, trace=FALSE)
          
          #fit complex normal distributions
          fit3 <- gamlss(data= fitdat, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + race_NBNW + smoke_Preg + smoke_unknown +
                           sex + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, 
                         family = NO, trace=FALSE)
          
          #fit simpler normal distrubutions
          fit4 <- gamlss(data= fitdat, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                          sex + CD4T + CD8T + Gran + nRBC, 
                         family = NO, trace=FALSE)
        
          predict1 <- predict(fit1, newdata = testdat, data = fitdat, type= "response")
          predict2 <- predict(fit2, newdata = testdat, data = fitdat, type= "response")
          predict3 <- predict(fit3, newdata = testdat, data = fitdat, type= "response")   
          predict4 <- predict(fit4, newdata = testdat, data = fitdat, type= "response")   
          
          myaics <- AIC(fit1, fit2, fit3, fit4)
          cpg <- namesRandom[k]
          colnames(myaics) <- c("df", paste0(cpg,"_",k))
          myaics$fit <- rownames(myaics)
          aics <- full_join(aics, myaics[2:3], by="fit")


          MSE1 <- mean((testdat$tmpexp - predict1)^2)
          mse1[j,k] <- MSE1

	  MSE2 <- mean((testdat$tmpexp - predict2)^2)
          mse2[j,k] <- MSE2

	  MSE3 <- mean((testdat$tmpexp - predict3)^2)
          mse3[j,k] <- MSE3

	  MSE4 <- mean((testdat$tmpexp - predict4)^2)
          mse4[j,k] <- MSE4

      } #end cpg for a fold
      } #end all cpgs for a fold
  
  colnames(mse1) <- colnames(randomBetaAll)
  myname <-  paste0("MSE1_", names(expAll[i]))
  assign(myname, mse1)
  write.csv(mse1, paste0("./output/",names(expAll[i]),"mse1.fit.csv"))

  colnames(mse2) <- colnames(randomBetaAll)
  myname <-  paste0("MSE2_", names(expAll[i]))
  assign(myname, mse2)
  write.csv(mse2, paste0("./output/",names(expAll[i]),"mse2.fit.csv"))  

  colnames(mse3) <- colnames(randomBetaAll)
  myname <-  paste0("MSE3_", names(expAll[i]))
  assign(myname, mse3)
  write.csv(mse3, paste0("./output/",names(expAll[i]),"mse3.fit.csv"))

  colnames(mse4) <- colnames(randomBetaAll)
  myname <-  paste0("MSE4_", names(expAll[i]))
  assign(myname, mse4)
  write.csv(mse4, paste0("./output/",names(expAll[i]),"mse4.fit.csv"))

  
  myname <-  paste0("AICs", names(expAll[i]))
  assign(myname, aics)
  write.csv(aics, paste0("./output/",names(expAll[i]),"AIC.fit.csv"))

  } #end all folds for an exposure - this is quick with 500 cpgs

print(paste0("Finished test data fiting."))

###########################################################################################################################

for (i in 1:length(expoOrder)){
    mseAves <- c()
  for (j in 1:4){
    mymse<- eval(parse(text=paste0("MSE",j,"_",expoOrder[[i]])))
    aves <- apply(mymse[1:5,], 2, mean)
    mseAves <-cbind(mseAves, aves)
    }
    tmpnm <- paste0("MSE_",expoOrder[[i]])
    assign(tmpnm, mseAves)
    write.csv(mseAves, paste0("./output/MSE_",expoOrder[[i]],"_aves.csv"))
}


remove(myks, namesRandom, rand, K1row, K2row, K3row, K4row, K4row, K5row, 
	mseAves, mymse, aves, tmpnm, singleExp, mse1, mse2, mse3, mse4, aics, index, tmpexp, tmpbeta, fitdat, 
	myaics, testbeta, newCpG, predict1, predict2, predict3, predict4, cpg, k, i, j, 
	myname, mymse, aves, mseAves, test, testdat,)
remove(list=apropos("MSE"))
remove(list=apropos("AIC"))

###########################################################################################################################
#calculate genomic lambdas

GL <- matrix(nrow=length(expAll), ncol=4)
rownames(GL) <- names(expAll)
colnames(GL) <- c("fit1", "fit2", "fit3", "fit4")

for (i in 1:length(expAll)){
  phenoSelect$tmpexp <-expAll[,i]
  chi <- matrix(nrow=length(randomBetaAll), ncol=4)

for (k in 1:length(randomBetaAll)){
          phenoSelect$tmpCpG <-randomBetaAll[,k]
                    
          #fit complex beta distributions
          fit1 <- gamlss(data= phenoSelect, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + race_NBNW + smoke_Preg + smoke_unknown +
              sex + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, 
              family = BE, trace=FALSE)
      
          #fit simpler beta distributions
          fit2 <- gamlss(data= phenoSelect, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                sex + CD4T + CD8T + Gran + nRBC, 
                family = BE, trace=FALSE)
          
          #fit complex normal distributions
          fit3 <- gamlss(data= phenoSelect, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + race_NBNW + smoke_Preg + smoke_unknown +
                           sex + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, 
                         family = NO, trace=FALSE)
          
          #fit simpler normal distrubutions
          fit4 <- gamlss(data= phenoSelect, tmpCpG ~  tmpexp + PC1 + PC2 + Parity + race_B + smoke_Preg +
                          sex + CD4T + CD8T + Gran + nRBC, 
                         family = NO, trace=FALSE)
        

	  sum1 <- summary(fit1, save=TRUE)
          chi[k,1] <- qchisq(1-sum1$coef.table[2,4],1)

 	  sum2 <- summary(fit2, save=TRUE)
          chi[k,2] <- qchisq(1-sum2$coef.table[2,4],1)

  	  sum3 <- summary(fit3, save=TRUE)
          chi[k,3] <- qchisq(1-sum3$coef.table[2,4],1)

 	  sum4 <- summary(fit4, save=TRUE)
          chi[k,4] <- qchisq(1-sum4$coef.table[2,4],1)

      } #end cpgs

        GL[i,] <- apply(chi, 2, function(x) median(x)/qchisq(0.5,1))
}

write.csv(GL, "genomic_lambdas.csv")

remove(i, k, fit1,fit2,fit3,fit4, sum1,sum2,sum3,sum4, chi)

save.image(file='session03312022.RData')
#end script - use files to visualize results and pick best models