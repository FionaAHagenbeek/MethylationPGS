################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: Run GEE analyses combined for all traits across all three cohorts
# (i.e., mega-analysis)
#
# data: Imputed and residualized top 10% most variable CpGs and PGSs (.RData)
#
# Author: F.A. Hagenbeek (fiona.hagenbeek@helsinki.fi)
#
################################################################################

################################################################################
#
# General set-up script
#
################################################################################

# ensure working dir is empty
rm(list=ls())

#load R packages
#install.packages("gee")
library(gee)
#install.packages("foreach")
library(foreach)
#install.packages("doMC")
library(doMC)
registerDoMC(10)  #number of cpu's to use


################################################################################
#
# Load & prepare data
#
################################################################################

# residualized and imputed DNA methylation betas and PGSs (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/final/2024-01-04_ImputedResidualDNAmPGSandCovs_alldats.RData")

# VERY IMPORTANT! Sort by familynumber! This is required for the GEE analysis!
mega2 <- combi2[order(combi2$FamilyNumber),]

# vector with CpG names
CpGs <- names(mega2)[45:length(mega2)]


################################################################################
#
# GEE functions  
#
################################################################################

#function to run GEE model
geemodel <- function(dat, CpG) { 
  r1 <- gee(as.formula(paste0(CpG, " ~ P_inf_SCORE_Height_MRG16_LDp1 + P_inf_SCORE_BMI_MRG16_LDp1 + 
                              P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 + P_inf_SCORE_Schizophrenia_MRG16_LDp1 + 
                              P_inf_SCORE_EducationalAttainment_MRG16_LDp1 + P_inf_SCORE_SocialDeprivation_MRG16_LDp1 + 
                              Moth_P_inf_SCORE_Height_MRG16_LDp1 + Moth_P_inf_SCORE_BMI_MRG16_LDp1 + 
                              Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 + Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1 + 
                              Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1 + Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1 + 
                              Fath_P_inf_SCORE_Height_MRG16_LDp1 + Fath_P_inf_SCORE_BMI_MRG16_LDp1 + 
                              Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 + Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1 + 
                              Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1 + Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1 + 
                              datasetdummy2")),
         data=dat, id=FamilyNumber, family=gaussian, corstr="independence", maxiter=500, na.action=na.omit,silent=TRUE)
  # extract output
  coeff <- summary(r1)$coefficients
  #add error term
  coeff <- cbind(coeff,r1$error)
  #output
  return(coeff)
}

extract.coeffs <- function(model.output,CpG,dat) {
  # create a list, were each list contains a vector of betas for each of the
  # tested covariates/covariate levels
  Beta <- foreach(i=1:nrow(model.output[[1]])) %dopar% {
    c(unlist(lapply(model.output, function(x) { x[i,1] })))
  }
  # convert list into a data frame
  Beta <- data.frame(Beta)
  # rename the columns to reflect the covariate it was obtained for
  names(Beta) <- c(paste0(row.names(model.output[[1]]),"_beta"))
  
  # create a list, were each list contains a vector of se's for each of the
  # tested covariates/covariate levels
  SE <- foreach(i=1:nrow(model.output[[1]])) %dopar% {
    c(unlist(lapply(model.output, function(x) { x[i,4] })))
  }
  # convert list into a data frame
  SE <- data.frame(SE)
  # rename the columns to reflect the covariate it was obtained for
  names(SE) <- c(paste0(row.names(model.output[[1]]),"_se"))
  
  # create a list, were each list contains a vector of p-values for each of the
  # tested covariates/covariate levels
  P <- foreach(i=1:nrow(model.output[[1]])) %dopar% {
    c(unlist(lapply(model.output, function(x) { 2*pnorm(-abs(x[i,5])) })))
  }
  # convert list into a data frame
  P <- data.frame(P)
  # rename the columns to reflect the covariate it was obtained for
  names(P) <- c(paste0(row.names(model.output[[1]]),"_p"))
  
  # create a list, where each list contains a vector of errors
  Error <- foreach(i=1:nrow(model.output[[1]])) %dopar% {
    c(unlist(lapply(model.output, function(x) { x[i,6] })))
  }
  # convert list into a data frame
  Error <- data.frame(Error)
  # rename columns
  names(Error) <- c(paste0(row.names(model.output[[1]]),"_error"))
  
  # create trait analyzed columns
  trait <- CpG
  
  # calculate the N per analysis
  N <- rep(NA,length(CpG))
  for (i in 1:length(CpG)) {
    N[i] <- nrow(na.omit(dat[,c(CpG[i],"P_inf_SCORE_Height_MRG16_LDp1", "P_inf_SCORE_BMI_MRG16_LDp1",
                                "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1", "P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                                "P_inf_SCORE_EducationalAttainment_MRG16_LDp1", "P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                                "Moth_P_inf_SCORE_Height_MRG16_LDp1", "Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                                "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1", "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                                "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1", "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1", 
                                "Fath_P_inf_SCORE_Height_MRG16_LDp1", "Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                                "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1", "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                                "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1", "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                                "datasetdummy1","datasetdummy2","FamilyNumber")]))
  }
  
  # combine the data frames with all coefficients into a single data frame
  out = cbind(trait, N, Beta, SE, P, Error)
  
  # function outputs the data frame with coefficients for each model in the model output list
  return(out)
}


################################################################################
#
# Run GEE models 
#
################################################################################

#run model and write results to list
res.model2  <- foreach(i=1:length(CpGs)) %dopar% {
  geemodel(dat = mega2, CpG = CpGs[i]) 
}

# extract output
res.table2 <- extract.coeffs(model.output = res.model2,CpG = CpGs, dat = mega2)

# rename columns to remove the unneeded PGS related info
names(res.table2) <- gsub("P_inf_SCORE_","", gsub("_MRG16_LDp1","",names(res.table2)))

# replace the multiple error columns with just one
res.table2 <- res.table2[,c(1:63)]
names(res.table2)[63] <- "error"

# write each of the results for each of the phenotypes to file
write.table(res.table2, file=paste0("/home/fhagenbeek/data/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                             "_EPIC_Buccal_megaanalysis_GEE_results_alltraits.txt"), row.names=F, col.names = T, sep="\t",quote = T) 
