################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: create final dataset for analysis including imputed residualized DNA
# methylation betas after excluding outliers and residualized PGSs
#
# data: 1. ACTION pilot, NTR ACTION and Avera Group 2 PGS + covariate files
# (.RData) 2. NTR ACTION estimates buccal cell counts (.RData) 3. ACTION pilot,
# NTR ACTION and Avera Group 2 DNA methylation betas (imputed) (.RData) 4. top
# 10% most variable DNA Methylation probes (.RData)
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
#install.packages("foreach")
library(foreach)
#install.packages("doMC")
library(doMC)
registerDoMC(15)  #number of cpu's to use
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("remotes")
# library(remotes)
# install_github("perishky/ewaff")
library(ewaff)
#install.packages("dplyr")
library(dplyr)
#install.packages("plyr")
#install.packages("ggplot2")
library(ggplot2)


################################################################################
#
# Load data
#
################################################################################

# DNA methylation covariate files (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/ACTIONpilot_PGSandCovs.RData") #load data ACTION pilot
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/NTRACTION_PGSandCovs.RData") #load data NTR ACTION 
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/AveraGroup2_PGSandCovs.RData") #load data Avera Group 2
#adjust names of sex variable Avera
names(dat.Avera)[14] <- "Sex"

#load annotation files 
epicanno <- get(load("/home/fhagenbeek/data/MethylationPGS/data/anno_epic_072017.RData"))

# Buccal cell counts NTR ACTION
load("/home/fhagenbeek/data/MethylationPGS/data/buccalestimates_EPIDISH_FN5pcs_BIOS.Robj") #load data
EPIdish_WBC <- as.data.frame(EPIdish_WBC) #as data frame

#read in imputed methylation betas
betas.pilot <- get(load("/data/PUBLIC/Methylation/ACTION_EPICarray_pilot/ACTION.EPIC.betas_pilot_nomissings.RData")) #load ACTION pilot
betas.ACTION <- get(load("/data/PUBLIC/Methylation/ACTION_EPICarray/ACTION.EPIC.betas_NTR_nomissings.RData")) #load betas NTR ACTION
betas.Avera <- get(load("/data/PUBLIC/Methylation/Avera_2022_EPICarray/Group2_Buccal/Betas_Avera2022_Group2_Buccal_012023_nomissings.RData")) #load betas Avera Group 2
rm(Beta.FunNorm) # remove duplicate betas 

# SD DNA methylation betas ACTION by batch (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/2023-12-20_EPIC_Buccal_SDprobes_alldats.RData")


################################################################################
#
# Order SD of CpGs by largest to smallest SD
# plot variability
#
################################################################################

#order by SD
CpG <- CpG[with(CpG, order(ressd,decreasing = TRUE)),]

# add column to indicate percentiles. The 'ntile' function is used to calculate
# percentiles. In this case, we are calculating percentiles in 10 groups, which
# represents percentiles from 1 to 10.
CpG <- CpG %>%
  mutate(percentile = ntile(ressd, 10))

# plot
tiff(filename = paste0("/home/fhagenbeek/data/MethylationPGS/data/intermediate/", as.character(Sys.Date()),
                       "_EPIC_Buccal_VariabilityCpGs.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
ggplot(CpG, aes(x=as.factor(percentile), y=ressd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE) + scale_fill_grey() + theme_classic() + 
  theme(legend.position="none") + labs(title="",x="Percentile", y = "Variability")
dev.off()


################################################################################
#
# write top 10% most variable probes to file
#
################################################################################

#select top 10% not imputed raw probes
CpG.top10 <- CpG[which(CpG$percentile==10),c("CpG","ressd")]

#annotate top10% most variable probes
CpG.top10.anno <- merge(CpG.top10,epicanno[,c("name","chromosome","position","gene.symbol","gene.accession","gene.region","cpg.island.name",
                                              "relation.to.island")], by.x = "CpG", by.y = "name", all.x = T)

#write output to file
save(CpG.top10.anno, file=paste("/home/fhagenbeek/data/MethylationPGS/data/intermediate/", as.character(Sys.Date()), "_EPIC_Buccal_Top10VariableCpGs_alldats.RData",sep=""))


################################################################################
#
# Merge PGS + covariate files
#
################################################################################

#add column with sample identifiers to buccal cell counts
EPIdish_WBC$IdatName <- rownames(EPIdish_WBC)

#merge buccal cell counts to covariates + PGSs
dat.ACTION <- merge(dat.ACTION, EPIdish_WBC[,c("IdatName","Epi","NK")], by = "IdatName", all.x = T)

# add dataset dummies 
dat.pilot$datasetdummy1 <- 0
dat.pilot$datasetdummy2 <- 0
dat.ACTION$datasetdummy1 <- 0
dat.ACTION$datasetdummy2 <- 1
dat.Avera$datasetdummy1 <- 1
dat.Avera$datasetdummy2 <- 1

# combine PGS + covariate files
covariates <- plyr::rbind.fill(dat.pilot[,c("IdatName","FISnumber","FamilyNumber","age_at_collection","Sex","Extension","multiple_type","twzyg",
                                      "Sample_Plate","Array_rownum","Epi","NK","PC1_1KG","PC2_1KG","PC3_1KG","PC4_1KG","PC5_1KG","PC6_1KG",
                                      "PC7_1KG","PC8_1KG","PC9_1KG","PC10_1KG","PLD_AXIOM","datasetdummy1","datasetdummy2",
                                      "P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1","P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Moth_P_inf_SCORE_Height_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Fath_P_inf_SCORE_Height_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                         dat.ACTION[,c("IdatName","FISnumber","FamilyNumber","age_at_collection","Sex","Extension","multiple_type","twzyg",
                                       "Sample_Plate","Array_rownum","Epi","NK","PC1_1KG","PC2_1KG","PC3_1KG","PC4_1KG","PC5_1KG","PC6_1KG",
                                       "PC7_1KG","PC8_1KG","PC9_1KG","PC10_1KG","PLD_AXIOM","PLD_GSA","datasetdummy1","datasetdummy2",
                                       "P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1","P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                       "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                       "P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Moth_P_inf_SCORE_Height_MRG16_LDp1",
                                       "Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                       "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                       "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Fath_P_inf_SCORE_Height_MRG16_LDp1",
                                       "Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                       "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                       "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                         dat.Avera[,c("IdatName","FISnumber","FamilyNumber","age_at_collection","Sex","Extension","multiple_type","twzyg",
                                      "Sample_Plate","Array_rownum","Epi","NK","PC1_1KG","PC2_1KG","PC3_1KG","PC4_1KG","PC5_1KG","PC6_1KG",
                                      "PC7_1KG","PC8_1KG","PC9_1KG","PC10_1KG","PLD_AXIOM","PLD_GSA","datasetdummy1","datasetdummy2",
                                      "P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1","P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Moth_P_inf_SCORE_Height_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Fath_P_inf_SCORE_Height_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                                      "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")])


################################################################################
#
# In DNA methylation beta files only retain top 10% most variable probes & those
# samples to be included in the analyses
#
################################################################################

# PILOT #
#keep overlapping samples
OS.pilot <- intersect(colnames(betas.pilot),dat.pilot$IdatName) #obtain overlap
betas.pilot <- betas.pilot[,OS.pilot] #only retain overlapping samples in the file with DNA methylation B-values
dim(betas.pilot) #check whether columns have been reduced

# ACTION #
#keep overlapping samples
OS.ACTION <- intersect(colnames(betas.ACTION),dat.ACTION$IdatName) #obtain overlap
betas.ACTION <- betas.ACTION[,OS.ACTION] #only retain overlapping samples in the file with DNA methylation B-values
dim(betas.ACTION) #check whether columns have been reduced

# AVERA #
#keep overlapping samples
OS.Avera <- intersect(colnames(betas.Avera),dat.Avera$IdatName) #obtain overlap
betas.Avera <- betas.Avera[,OS.Avera] #only retain overlapping samples in the file with DNA methylation B-values
dim(betas.Avera) #check whether columns have been reduced


# retain only overlapping betas
# pilot
OL.pilot <- intersect(rownames(betas.pilot),CpG.top10.anno$CpG)
betas.pilot <- betas.pilot[OL.pilot,] #betas ACTION pilot

# NTR ACTION
OL.ACTION <- intersect(rownames(betas.ACTION),CpG.top10.anno$CpG)
betas.ACTION <- betas.ACTION[OL.ACTION,] #betas NTR ACTION

# Avera
OL.Avera <- intersect(rownames(betas.Avera),CpG.top10.anno$CpG)
betas.Avera <- betas.Avera[OL.Avera,] #betas Avera Group 2

# remove outliers
methylation.pilot <- ewaff.handle.outliers(betas.pilot, method="iqr")[[1]] # pilot
methylation.ACTION <- ewaff.handle.outliers(betas.ACTION, method="iqr")[[1]] # NTR ACTION
methylation.Avera <- ewaff.handle.outliers(betas.Avera, method="iqr")[[1]] # Avera

# combine betas
methylation <- cbind(methylation.pilot,methylation.ACTION,methylation.Avera)


################################################################################
#
# Merge overlapping top 10% most variable residualized DNA methylation betas to 
# PGS and covariates for ACTION pilot, NTR ACTION and Avera Group 2
#
################################################################################

# function to obtain the scaled residualized DNA methylation beta values
resid = function(dat,probe) { 
  reg=lm(as.formula(paste0(probe," ~ age + sex + Sample_Plate + Array_rownum + Epi + NK")), data = dat, na.action=na.exclude)
  out=as.data.frame(reg$residuals)
  colnames(out) <- probe
  out$IdatName <- rownames(out)
  return(out)
}

# select variables
age <- covariates$age_at_collection
sex <- as.numeric(factor(covariates$Sex, levels = c("M","F")))-1
Sample_Plate <- as.factor(covariates$Sample_Plate)
Array_rownum <- covariates$Array_rownum
Epi <- covariates$Epi
NK <- covariates$NK

#create 1 dataframe with betas and covariates
mega <- data.frame(as.data.frame(t(methylation)),age,sex,Sample_Plate,Array_rownum,Epi,NK)

# create empty  matrix to store residualized beta values
CpG.dat <- as.data.frame(matrix(NA,ncol=0,nrow=ncol(methylation))) # create empty dataframe
CpG.dat$IdatName <- colnames(methylation) # add sample identifier  to dataframe

#run function to obtain the residuals of each CpG
resid.dat  <- foreach(i=1:nrow(methylation)) %dopar% {
  resid(dat = mega, probe = CpG.top10.anno$CpG[i])
}

#add residualized values to temporary dataframe
for (i in 1:nrow(CpG.top10.anno)) {
  CpG.dat <- merge(CpG.dat,resid.dat[[i]], by = "IdatName",all.x = T)
}


################################################################################
#
# Obtain residualized PGSs for ACTION pilot, NTR ACTION and Avera Group 2
#
################################################################################

#vector with PGSs
PGSs <- c("P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1","P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
          "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
          "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
          "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
          "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
          "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")

#vector with 1st 10 PCs
PCs <- c("PC1_1KG","PC2_1KG","PC3_1KG","PC4_1KG","PC5_1KG","PC6_1KG","PC7_1KG","PC8_1KG","PC9_1KG","PC10_1KG")

# PILOT #
# select variables
PLD_AXIOM <- dat.pilot$PLD_AXIOM

#create 1 dataframe with betas and covariates
pilot <- data.frame(as.data.frame(dat.pilot[,c(PGSs,PCs)]),PLD_AXIOM)
rownames(pilot) <- dat.pilot$IdatName

# create empty  matrix to store residualized beta values
PGS.pilot <- as.data.frame(matrix(NA,ncol=0,nrow=nrow(pilot))) # create empty dataframe
PGS.pilot$IdatName <- rownames(pilot) # add sample identifier  to dataframe

# function to obtain the residualized PGSs
resid.pgs = function(dat,PGS) { 
  reg=lm(as.formula(paste0(PGS," ~ PLD_AXIOM + PC1_1KG + PC2_1KG + PC3_1KG + PC4_1KG + PC5_1KG + PC6_1KG + PC7_1KG + PC8_1KG + PC9_1KG + 
                           PC10_1KG")), data = dat, na.action=na.exclude)
  out=as.data.frame(scale(reg$residuals))
  colnames(out) <- PGS
  out$IdatName <- rownames(out)
  return(out)
}

#run function to obtain the residuals of each PGS
resid.pilot.pgs  <- foreach(i=1:length(PGSs)) %dopar% {
  resid.pgs(dat = pilot, PGS = PGSs[i])
}

#add residualized values to temporary dataframe
for (i in 1:length(PGSs)) {
  PGS.pilot <- merge(PGS.pilot,resid.pilot.pgs[[i]], by = "IdatName",all.x = T)
}

#merge residualized betas to PGSs + covariates
dat.pilot <- merge(dat.pilot[,!names(dat.pilot) %in% PGSs],PGS.pilot, by="IdatName") #merge
rm(PLD_AXIOM,pilot,resid.pilot.pgs,PGS.pilot) # remove vectors and dataframes that are no longer required

# ACTION #
# select variables
PLD_AXIOM <- dat.ACTION$PLD_AXIOM
PLD_GSA <- dat.ACTION$PLD_GSA

#create 1 dataframe with betas and covariates
ACTION <- data.frame(as.data.frame(dat.ACTION[,c(PGSs,PCs)]),PLD_AXIOM,PLD_GSA)
rownames(ACTION) <- dat.ACTION$IdatName

# create empty  matrix to store residualized beta values
PGS.ACTION <- as.data.frame(matrix(NA,ncol=0,nrow=nrow(ACTION))) # create empty dataframe
PGS.ACTION$IdatName <- rownames(ACTION) # add sample identifier  to dataframe

# function to obtain the residualized PGSs
resid.pgs = function(dat,PGS) { 
  reg=lm(as.formula(paste0(PGS," ~ PLD_AXIOM + PLD_GSA + PC1_1KG + PC2_1KG + PC3_1KG + PC4_1KG + PC5_1KG + PC6_1KG + PC7_1KG + PC8_1KG + 
                           PC9_1KG + PC10_1KG")), data = dat, na.action=na.exclude)
  out=as.data.frame(scale(reg$residuals))
  colnames(out) <- PGS
  out$IdatName <- rownames(out)
  return(out)
}

#run function to obtain the residuals of each PGS
resid.ACTION.pgs  <- foreach(i=1:length(PGSs)) %dopar% {
  resid.pgs(dat = ACTION, PGS = PGSs[i])
}

#add residualized values to temporary dataframe
for (i in 1:length(PGSs)) {
  PGS.ACTION <- merge(PGS.ACTION,resid.ACTION.pgs[[i]], by = "IdatName",all.x = T)
}

#merge residualized betas to PGSs + covariates
dat.ACTION <- merge(dat.ACTION[,!names(dat.ACTION) %in% PGSs],PGS.ACTION, by="IdatName") #merge
rm(PLD_AXIOM,PLD_GSA,ACTION,resid.ACTION.pgs,PGS.ACTION) # remove vectors and dataframes that are no longer required

# AVERA #
# select variables
PLD_AXIOM <- dat.Avera$PLD_AXIOM
PLD_GSA <- dat.Avera$PLD_GSA

#create 1 dataframe with betas and covariates
Avera <- data.frame(as.data.frame(dat.Avera[,c(PGSs,PCs)]),PLD_AXIOM,PLD_GSA)
rownames(Avera) <- dat.Avera$IdatName

# create empty  matrix to store residualized beta values
PGS.Avera <- as.data.frame(matrix(NA,ncol=0,nrow=nrow(Avera))) # create empty dataframe
PGS.Avera$IdatName <- rownames(Avera) # add sample identifier  to dataframe

# function to obtain the residualized PGSs
resid.pgs = function(dat,PGS) { 
  reg=lm(as.formula(paste0(PGS," ~ PLD_AXIOM + PLD_GSA + PC1_1KG + PC2_1KG + PC3_1KG + PC4_1KG + PC5_1KG + PC6_1KG + PC7_1KG + PC8_1KG + 
                           PC9_1KG + PC10_1KG")), data = dat, na.action=na.exclude)
  out=as.data.frame(scale(reg$residuals))
  colnames(out) <- PGS
  out$IdatName <- rownames(out)
  return(out)
}

#run function to obtain the residuals of each PGS
resid.Avera.pgs  <- foreach(i=1:length(PGSs)) %dopar% {
  resid.pgs(dat = Avera, PGS = PGSs[i])
}

#add residualized values to temporary dataframe
for (i in 1:length(PGSs)) {
  PGS.Avera <- merge(PGS.Avera,resid.Avera.pgs[[i]], by = "IdatName",all.x = T)
}

#merge residualized betas to PGSs + covariates
dat.Avera <- merge(dat.Avera[,!names(dat.Avera) %in% PGSs],PGS.Avera, by="IdatName") #merge
rm(PLD_AXIOM,PLD_GSA,Avera,resid.Avera.pgs,PGS.Avera) # remove vectors and dataframes that are no longer required

# combine PGSs 
PGS.dat <- rbind(dat.pilot[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                              "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Moth_P_inf_SCORE_Height_MRG16_LDp1",
                              "Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Fath_P_inf_SCORE_Height_MRG16_LDp1",
                              "Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                 dat.ACTION[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                               "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                               "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                               "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                               "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                               "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                 dat.Avera[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                              "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")])

# combine PGS + covariate files
covariates <- merge(covariates[,!names(covariates) %in% PGSs],PGS.dat, by="IdatName") #merge

#merge residualized betas to PGSs + covariates
combi <- merge(covariates,CpG.dat, by="IdatName") #merge

# save file prior to imputation of parental PGSs
save(combi, file=paste("/home/fhagenbeek/data/MethylationPGS/data/final/", as.character(Sys.Date()), 
                           "_ImputedResidualDNAm_ResidualPGSandCovs_alldats.RData",sep="")) 

################################################################################
#
# Impute missing parental PGss with zero (Kong et al. 2018)
#
################################################################################

# vector with parental PGSs
parental.PGSs <- c("Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                   "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                   "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                   "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")

# obtain N parents without genotype information (e.g., one or both parents have
# no genotype information and should be imputed to 0) - only run this for one
# the PGSs as the missingness is the same across PGSs
## PILOT ##
sum(is.na(dat.pilot$Fath_P_inf_SCORE_BMI_MRG16_LDp1) | is.na(dat.pilot$Moth_P_inf_SCORE_BMI_MRG16_LDp1)) # N = 31
## ACTION ##
sum(is.na(dat.ACTION$Fath_P_inf_SCORE_BMI_MRG16_LDp1) | is.na(dat.ACTION$Moth_P_inf_SCORE_BMI_MRG16_LDp1)) # N = 213
## Avera ##
sum(is.na(dat.Avera$Fath_P_inf_SCORE_BMI_MRG16_LDp1) | is.na(dat.Avera$Moth_P_inf_SCORE_BMI_MRG16_LDp1)) # N = 87

# if parental PGS is.na impute with 0
## PILOT ##
for(i in 1:nrow(dat.pilot)) {
  for(j in 1:length(parental.PGSs)) {
    if(is.na(dat.pilot[i,parental.PGSs[j]])) {
      dat.pilot[i,parental.PGSs[j]] <- 0
    }
  }
}

## NTR ACTION ##
for(i in 1:nrow(dat.ACTION)) {
  for(j in 1:length(parental.PGSs)) {
    if(is.na(dat.ACTION[i,parental.PGSs[j]])) {
      dat.ACTION[i,parental.PGSs[j]] <- 0
    }
  }
}

## Avera ##
for(i in 1:nrow(dat.Avera)) {
  for(j in 1:length(parental.PGSs)) {
    if(is.na(dat.Avera[i,parental.PGSs[j]])) {
      dat.Avera[i,parental.PGSs[j]] <- 0
    }
  }
}

# combine PGSs 
PGS2.dat <- rbind(dat.pilot[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                              "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "P_inf_SCORE_Schizophrenia_MRG16_LDp1","P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Moth_P_inf_SCORE_Height_MRG16_LDp1",
                              "Moth_P_inf_SCORE_BMI_MRG16_LDp1","Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1","Fath_P_inf_SCORE_Height_MRG16_LDp1",
                              "Fath_P_inf_SCORE_BMI_MRG16_LDp1","Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1","Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1",
                              "Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                 dat.ACTION[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                               "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                               "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                               "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                               "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                               "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                               "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")],
                 dat.Avera[,c("IdatName","P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",
                              "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                              "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                              "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                              "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")])

# combine PGS + covariate files
covariates2 <- merge(covariates[,!names(covariates) %in% PGSs],PGS2.dat, by="IdatName") #merge

#merge residualized betas to PGSs + covariates
combi2 <- merge(covariates2,CpG.dat, by="IdatName") #merge


#write final file for analysis to file
save(combi2, file=paste("/home/fhagenbeek/data/MethylationPGS/data/final/", as.character(Sys.Date()), 
                           "_ImputedResidualDNAmPGSandCovs_alldats.RData",sep="")) 
