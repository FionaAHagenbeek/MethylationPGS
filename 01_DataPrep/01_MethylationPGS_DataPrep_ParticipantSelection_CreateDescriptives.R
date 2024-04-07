################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: Combine data types (see data), select participants for the analyses,
# and create participant descriptive table
#
# data: 1. ACTION pilot DNA Methylation, ACTION NTR DNA Methylation, Avera Group
# 2 DNA Methylation covariates (RData) 2. DNA methylation pedigree file (.sav)
# 3. NTR MRG16 pedigree file (.sav) 4. PGSs files (.sav) for Height + BMI (Yengo
# '18), Smoking initiation (Liu '19), Schizophrenia (Trubetskoy '22),
# Educational Attainment (Lee '18), and Social Deprivation (Hill '16) 5. basic
# covariates (.sav file)
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
#install.packages("foreign")
library(foreign)
#install.packages("dplyr)
library(dplyr)

#set working directory
setwd("/home/fhagenbeek/data/MethylationPGS/")


################################################################################
#
# Load data
#
################################################################################

# DNA methylation covariate files (ACTION pilot, NTR ACTION, and Avera Group2; .RData)
# pilot
load("data/Covariates_ACTION_EPICarray_pilot_DSR_3454_21032023.RData") #load data
pilot <- covariates #rename file
rm(covariates) #remove original file
#NTR ACTION
load("data/Covariates_NTR_Curium_V3.RData") #load data
ACTION <- covariates #rename file
rm(covariates) #remove original file
#Avera Group2
load("data/Covariates_EPIC_Avera_Group2_Buccal_NTR_DSR_3454_09032023.RData") #load data
Avera <- covariates #rename file
rm(covariates) #remove original file

# DNA methylation pedigree file (.sav)
ped <- read.spss(file = "data/BuccalEPICarray_NTRACTION_and_AVERA2022_Pedigree_18112022.sav",
                 to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)

# MRG16 (genotype) pedigree file (.sav)
MRG16 <- read.spss(file = "data/MRG16_Pedigree_V1.sav",
                  to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
#strip white space from mother and father identifiers
MRG16$MOTH <- trimws(MRG16$MOTH, which = "both")
MRG16$FATH <- trimws(MRG16$FATH, which = "both")


# PGS files for height, BMI, smoking initiation, schizophrenia, educational
# attainment and social deprivation (.sav)
Height <- read.spss(file = "data/NTR-DSR-3454_Height_PMID30124842_MRG16_PedMergedWithScores.sav",
                    to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
BMI <- read.spss(file = "data/NTR-DSR-3454_BMI_PMID30124842_MRG16_PedMergedWithScores.sav",
                 to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
Smoking <- read.spss(file = "data/NTR-DSR-3454_Smoking_SmokingInitiation_PMID30643251_MRG16_PedMergedWithScores.sav",
                     to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
Schizo <- read.spss(file = "data/NTR-DSR-3454_Schizophrenia_DOI10.1101_FSLASH_2020.09.12.20192922_MRG16_PedMergedWithScores.sav",
                    to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
EA <- read.spss(file = "data/NTR-DSR-3454_EducationalAttainment_PMID30038396_MRG16_PedMergedWithScores.sav",
                to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)
SocDep <- read.spss(file = "data/NTR-DSR-3454_SocialDeprivation_PMID27818178_MRG16_PedMergedWithScores.sav",
                    to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)

# basic covariates (.sav)
Pheno <- read.spss(file = "data/PHE_20220727_3454_TVB_FionaHagenbeek_FIS.sav",
                     to.data.frame = T,trim.factor.names = T,reencode = T, use.missings = T)


################################################################################
#
# For each child create maternal and paternal PGS per data set 
# (ACTION pilot, NTR ACTION, Avera Group2)
#
################################################################################

#to further simplify merging, combine PGss into 1 file while only retaining
#necessary variables combine PGSs in list
list_data <- list(Height,BMI,Smoking,Schizo,EA,SocDep)
# in PGS files only retain person identifier, Ldpred PGS with causal fraction = Inf
list_data = lapply(list_data, function(x) select(x, c(FISnumber, starts_with("P_inf_"))))
# merge PGSs into 1 file 
PGS <- Reduce(function(x,y) merge(x = x, y = y, by = "FISnumber",all=T), list_data)
# add mother & father identifiers to PGS overview
PGS <- merge(PGS,Height[,c("FISnumber","FATH","MOTH","EUR_1KG_Outlier")])
#strip white space from mother and father identifiers
PGS$MOTH <- trimws(PGS$MOTH, which = "both")
PGS$FATH <- trimws(PGS$FATH, which = "both")

## for each child extract PGS + ancestry outlier information for mother ##

#PGS height mother
PGS$Moth_P_inf_SCORE_Height_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_Height_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Height_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
# PGS BMI mother
PGS$Moth_P_inf_SCORE_BMI_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_BMI_MRG16_LDp1[i] <- PGS$P_inf_SCORE_BMI_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
# PGS Smoking mother
PGS$Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
#PGS Schizo mother
PGS$Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Schizophrenia_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
#PGS EA mother
PGS$Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1[i] <- PGS$P_inf_SCORE_EducationalAttainment_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
# PGS SocDep mother
PGS$Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1 <- NA #create column with NA's to hold maternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for mother if mother identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1[i] <- PGS$P_inf_SCORE_SocialDeprivation_MRG16_LDp1[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}
# outlier information for mother
PGS$Moth_EUR_1KG_Outlier <- NA #create column with NA's to hold maternal anstry outlier dummy
for (i in 1:nrow(PGS)) { #extract ancestry outlier dummy for father if father identifier occurs in the list with person identifiers
  if(PGS$MOTH[i] %in% PGS$FISnumber) {
    PGS$Moth_EUR_1KG_Outlier[i] <- PGS$EUR_1KG_Outlier[which(PGS$FISnumber==PGS$MOTH[i])]
  }
}

## for each child extract PGS + ancestry outlier information for father ##

#PGS height father
PGS$Fath_P_inf_SCORE_Height_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_Height_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Height_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
# PGS BMI father
PGS$Fath_P_inf_SCORE_BMI_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_BMI_MRG16_LDp1[i] <- PGS$P_inf_SCORE_BMI_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
# PGS Smoking father
PGS$Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
#PGS Schizo father
PGS$Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1[i] <- PGS$P_inf_SCORE_Schizophrenia_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
#PGS EA father
PGS$Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1[i] <- PGS$P_inf_SCORE_EducationalAttainment_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
#PGS SocDep father
PGS$Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1 <- NA #create column with NA's to hold paternal PGS
for (i in 1:nrow(PGS)) { #extract PGS for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1[i] <- PGS$P_inf_SCORE_SocialDeprivation_MRG16_LDp1[which(PGS$FISnumber==PGS$FATH[i])]
  }
}
#outlier information for father
PGS$Fath_EUR_1KG_Outlier <- NA #create column with NA's to hold paternal ancestry outlier dummy
for (i in 1:nrow(PGS)) { #extract ancestry outlier dummy for father if father identifier occurs in the list with person identifiers
  if(PGS$FATH[i] %in% PGS$FISnumber) {
    PGS$Fath_EUR_1KG_Outlier[i] <- PGS$EUR_1KG_Outlier[which(PGS$FISnumber==PGS$FATH[i])]
  }
}

#save a copy of the new file as intermediate file
save(PGS, file="data/intermediate/OffSpringAndParental_PGSs_HeightBMISmokingSchizophreniaEASocialdeprivation.RData")

################################################################################
#
# Merge DNA methylation covariate, PGS & basic covariate files per data set 
# (ACTION pilot, NTR ACTION, Avera Group2)
#
################################################################################

#to simplify merging, ensure the name of the person identifier is the same in all files
names(pilot)[names(pilot) == "fisnumber"] <- "FISnumber"
names(ACTION)[names(ACTION) == "SubjectFISNumber"] <- "FISnumber"
names(Avera)[names(Avera) == "fisnumber"] <- "FISnumber"
names(ped)[names(ped) == "fisnumber"] <- "FISnumber"
names(Pheno)[names(Pheno) == "FISNumber"] <- "FISnumber"

#to simplify merging, ensure age at DNA sampling is the same in all files
names(Avera)[names(Avera) == "DNA_Sampling_Age"] <- "age_at_collection"

#create merged ACTION pilot file
dat.pilot <- Reduce(function(x,y) merge(x = x, y = y, by = "FISnumber",all.x=T), 
                    list(pilot,ped[,c("FISnumber","ACTION_EPIC_buccal_Pilot")],MRG16,PGS[,-c(8:10)],Pheno[,c(1,15:ncol(Pheno))]))

#create merged NTR ACTION file
dat.ACTION <- Reduce(function(x,y) merge(x = x, y = y, by = "FISnumber",all.x=T), 
                     list(ACTION,ped[,c("FISnumber","ACTION_EPIC_buccal")],MRG16,PGS[,-c(8:10)],Pheno[,c(1,15:ncol(Pheno))]))

#create merged Avera Group 2 file
dat.Avera <- Reduce(function(x,y) merge(x = x, y = y, by = "FISnumber",all.x=T), 
                    list(Avera,ped[,c("FISnumber","AVERA_EPIC_2022_Group2_buccal_afterQC")],MRG16,PGS[,-c(8:10)],Pheno[,c(1,15:ncol(Pheno))]))


################################################################################
#
# Exclude participants per dataset (ACTION pilot, NTR ACTION, Avera Group2)
#   0. NTR ACTION Biomarker Study: participant is not part of NTR
#   1. Consent for participation has been withdrawn after generation of omics data
#   2. Removal of duplicate samples
#   3. Child has no PGS data
#   4. Child or parent(s) of child are identified as ancestry outlier
#   5. Avera group 2: remove participants >13 years of age
#
################################################################################

## PILOT ##

# 1. no consent
#are there any participants without consent? 
sum(dat.pilot$ACTION_EPIC_buccal_Pilot=="NO",na.rm = T) #no

# 2. remove duplicate samples
#are there any duplicate samples?
sum(duplicated(dat.pilot$FISnumber),na.rm = T) #2 duplicate measures
#remove duplicate samples
dat.pilot <- dat.pilot[!duplicated(dat.pilot$FISnumber), ] #remove
dim(dat.pilot) #sanity check - indeed 2 rows removed

# 3. child has no PGS
#are there children without PGSs?
sum(is.na(dat.pilot$P_inf_SCORE_Height_MRG16_LDp1),na.rm = T) #2 children without PGSs
#remove children without PGSs
dat.pilot <- dat.pilot[-which(is.na(dat.pilot$P_inf_SCORE_Height_MRG16_LDp1)), ] #remove
dim(dat.pilot) #sanity check - indeed 2 rows removed
sum(as.data.frame(table(dat.pilot$FamilyNumber))$Freq==2) #N complete twin pairs with PGSs = 48
sum(!is.na(unique(as.numeric(dat.pilot$MOTH[which(!is.na(dat.pilot$Moth_P_inf_SCORE_Height_MRG16_LDp1))])))) #43 mothers have PGS 
sum(!is.na(unique(as.numeric(dat.pilot$FATH[which(!is.na(dat.pilot$Fath_P_inf_SCORE_Height_MRG16_LDp1))])))) #36 fathers have PGS

# 4. remove ancestry outliers (parents + offspring)
#how many children are identified as ancestry outliers?
sum(dat.pilot$EUR_1KG_Outlier==1,na.rm = T) #2 children identified as ancestry outliers
#remove children identified as ancestry outliers
dat.pilot <- dat.pilot[-which(dat.pilot$EUR_1KG_Outlier==1), ] #remove
dim(dat.pilot) #sanity check - indeed 2 rows removed
sum(as.data.frame(table(dat.pilot$FamilyNumber))$Freq==2) #N complete non-ancestry outlier twin pairs with PGSs = 47

#after removal of ancestry outliers in the offpsring generation, do we retain any outlier in the parental generation?
sum(dat.pilot$Moth_EUR_1KG_Outlier==1,na.rm = T) #0 mothers identified as ancestry outliers
sum(dat.pilot$Fath_EUR_1KG_Outlier==1,na.rm = T) #0 fathers identified as ancestry outliers

#save a copy of the new file as intermediate file
save(dat.pilot, file="data/intermediate/ACTIONpilot_PGSandCovs.RData")

## ACTION ##

# 0. remove Curium participants
#how many participants in the ACTION Biomarker Study are from LUMC-Curium and not NTR?
sum(dat.ACTION$cohort=="curium") #187 participants
#remove Curium participants
dat.ACTION <- dat.ACTION[-which(dat.ACTION$cohort=="curium"), ] #remove
dim(dat.ACTION) #sanity check - indeed 187 rows removed

# 1. no consent
#are there any participants without consent? 
sum(is.na(dat.ACTION$ACTION_EPIC_buccal)) #2 individuals
#remove children without consent
dat.ACTION <- dat.ACTION[-which(is.na(dat.ACTION$ACTION_EPIC_buccal)), ] #remove
dim(dat.ACTION) #sanity check - indeed 2 rows removed

# 2. remove duplicate samples
#are there any duplicate samples?
sum(duplicated(dat.ACTION$FISnumber),na.rm = T) #2 duplicate measures
#remove duplicate samples
dat.ACTION <- dat.ACTION[!duplicated(dat.ACTION$FISnumber), ] #remove
dim(dat.ACTION) #sanity check - indeed 2 rows removed

# 3. child has no PGS
#are there children without PGSs?
sum(is.na(dat.ACTION$P_inf_SCORE_Height_MRG16_LDp1),na.rm = T) #12 children without PGSs
#remove children without PGSs
dat.ACTION <- dat.ACTION[-which(is.na(dat.ACTION$P_inf_SCORE_Height_MRG16_LDp1)), ] #remove
dim(dat.ACTION) #sanity check - indeed 12 rows removed
sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4) #N complete twin pairs with PGSs = 558
sum(!is.na(unique(as.numeric(dat.ACTION$MOTH[which(!is.na(dat.ACTION$Moth_P_inf_SCORE_Height_MRG16_LDp1))])))) #629 mothers have PGS 
sum(!is.na(unique(as.numeric(dat.ACTION$FATH[which(!is.na(dat.ACTION$Fath_P_inf_SCORE_Height_MRG16_LDp1))])))) #550 mothers have PGS

# 4. remove ancestry outliers (parents + offspring)
#how many children are identified as ancestry outliers?
sum(dat.ACTION$EUR_1KG_Outlier==1,na.rm = T) #67 children identified as ancestry outliers
#remove children identified as ancestry outliers
dat.ACTION <- dat.ACTION[-which(dat.ACTION$EUR_1KG_Outlier==1), ] #remove
dim(dat.ACTION) #sanity check - indeed 67 rows removed
sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4) #N complete non-ancestry outlier twin pairs with PGSs = 528

#after removal of ancestry outliers in the offpsring generation, do we retain any outlier in the parental generation?
sum(dat.ACTION$Moth_EUR_1KG_Outlier==1,na.rm = T) #10 twins with mothers identified as ancestry outliers
sum(dat.ACTION$Fath_EUR_1KG_Outlier==1,na.rm = T) #3 twins with fathers identified as ancestry outliers
#remove mothers identified as ancestry outliers
dat.ACTION <- dat.ACTION[-which(dat.ACTION$Moth_EUR_1KG_Outlier==1), ] #remove
dim(dat.ACTION) #sanity check - indeed 10 rows removed
sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4) #N complete non-ancestry outlier twin pairs with PGSs = 524
#after removal of ancestry outliers in the offpsring generation & twins whose mothers are ancestry outliers, do we retain any outlier in the parental generation?
sum(dat.ACTION$Fath_EUR_1KG_Outlier==1,na.rm = T) #3 twins with fathers identified as ancestry outliers
#remove mothers identified as ancestry outliers
dat.ACTION <- dat.ACTION[-which(dat.ACTION$Fath_EUR_1KG_Outlier==1), ] #remove
dim(dat.ACTION) #sanity check - indeed 3 rows removed
sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4) #N complete non-ancestry outlier twin pairs with PGSs = 523

#save a copy of the new file as intermediate file
save(dat.ACTION, file="data/intermediate/NTRACTION_PGSandCovs.RData")

## AVERA ##

# N type of participants prior to exclussion
table(dat.Avera$buccalsubgroup)

# 1. no consent
#are there any participants without consent? 
sum(is.na(dat.Avera$AVERA_EPIC_2022_Group2_buccal_afterQC)) #0 individuals

# 2. remove duplicate samples
#are there any duplicate samples?
sum(duplicated(dat.Avera$FISnumber),na.rm = T) #3 duplicate measures
#remove duplicate samples
dat.Avera <- dat.Avera[!duplicated(dat.Avera$FISnumber), ] #remove
dim(dat.Avera) #sanity check - indeed 3 rows removed

# 3. child has no PGS
#are there children without PGSs?
sum(is.na(dat.Avera$P_inf_SCORE_Height_MRG16_LDp1),na.rm = T) #7 children without PGSs
#remove children without PGSs
dat.Avera <- dat.Avera[-which(is.na(dat.Avera$P_inf_SCORE_Height_MRG16_LDp1)), ] #remove
dim(dat.Avera) #sanity check - indeed 7 rows removed
table(dat.Avera$buccalsubgroup) # overview of children with PGSs per subgroup
sum(!is.na(unique(as.numeric(dat.Avera$MOTH[which(!is.na(dat.Avera$Moth_P_inf_SCORE_Height_MRG16_LDp1))])))) #186 mothers have PGS 
sum(!is.na(unique(as.numeric(dat.Avera$FATH[which(!is.na(dat.Avera$Fath_P_inf_SCORE_Height_MRG16_LDp1))])))) #166 fathers have PGS

# 4. remove ancestry outliers (parents + offspring)
#how many children are identified as ancestry outliers?
sum(dat.Avera$EUR_1KG_Outlier==1,na.rm = T) #12 children identified as ancestry outliers
#remove children identified as ancestry outliers
dat.Avera <- dat.Avera[-which(dat.Avera$EUR_1KG_Outlier==1), ] #remove
dim(dat.Avera) #sanity check - indeed 12 rows removed

#after removal of ancestry outliers in the offpsring generation, do we retain any outlier in the parental generation?
sum(dat.Avera$Moth_EUR_1KG_Outlier==1,na.rm = T) #1 mothers identified as ancestry outliers
sum(dat.Avera$Fath_EUR_1KG_Outlier==1,na.rm = T) #0 fathers identified as ancestry outliers
#remove mothers identified as ancestry outliers
dat.Avera <- dat.Avera[-which(dat.Avera$Moth_EUR_1KG_Outlier==1), ] #remove
dim(dat.Avera) #sanity check - indeed 1 rows removed

# 5. remove participants >18 years of age
# how many children are >18 years of age?
sum(dat.Avera$age_at_collection>18) #65 children
#remove participants >18 years of age
dat.Avera <- dat.Avera[-which(dat.Avera$age_at_collection>18), ] #remove
dim(dat.Avera) #sanity check - indeed 1 rows removed

#save a copy of the new file as intermediate file
save(dat.Avera, file="data/intermediate/AveraGroup2_PGSandCovs.RData")

################################################################################
#
# Obtain descriptive statistics for each dataset + totals
#
################################################################################

#create empty dataframe
overviewmaxN <- as.data.frame(matrix(NA,nrow = 4,ncol = 0))
#add dataset indicitor
overviewmaxN$dataset <- c("dataset 1", "dataset 2", "dataset 3", "total")
# add N 
overviewmaxN$N_total <- c(nrow(dat.pilot),
                          nrow(dat.ACTION),
                          nrow(dat.Avera),
                          nrow(dat.pilot)+nrow(dat.ACTION)+nrow(dat.Avera))

# Add age information per dataset
overviewmaxN$M_SD_range_Age <- c(paste0(round(mean(dat.pilot$age_at_collection,na.rm = T),1)," (", 
                                        round(sd(dat.pilot$age_at_collection,na.rm = T),1),")", " [",
                                        round(min(dat.pilot$age_at_collection,na.rm = T),1), " - ", 
                                        round(max(dat.pilot$age_at_collection, na.rm = T),1), "]"),
                                 paste0(round(mean(dat.ACTION$age_at_collection,na.rm = T),1)," (", 
                                        round(sd(dat.ACTION$age_at_collection,na.rm = T),1),")", " [",
                                        round(min(dat.ACTION$age_at_collection,na.rm = T),1), " - ", 
                                        round(max(dat.ACTION$age_at_collection, na.rm = T),1), "]"),
                                 paste0(round(mean(dat.Avera$age_at_collection,na.rm = T),1)," (", 
                                        round(sd(dat.Avera$age_at_collection,na.rm = T),1),")", " [",
                                        round(min(dat.Avera$age_at_collection,na.rm = T),1), " - ", 
                                        round(max(dat.Avera$age_at_collection, na.rm = T),1), "]"),
                                 paste0(round(mean(c(dat.pilot$age_at_collection,dat.ACTION$age_at_collection,dat.Avera$age_at_collection),na.rm = T),1)," (", 
                                        round(sd(c(dat.pilot$age_at_collection,dat.ACTION$age_at_collection,dat.Avera$age_at_collection),na.rm = T),1),")", " [",
                                        round(min(c(dat.pilot$age_at_collection,dat.ACTION$age_at_collection,dat.Avera$age_at_collection),na.rm = T),1), " - ", 
                                        round(max(c(dat.pilot$age_at_collection,dat.ACTION$age_at_collection,dat.Avera$age_at_collection), na.rm = T),1), "]"))

# add gender information per dataset
overviewmaxN$N_males <- c(paste0(sum(dat.pilot$SEX==1), " (", 
                                 round(prop.table(table(dat.pilot$SEX))[1]*100,1), "%)"),
                          paste0(sum(dat.ACTION$SEX==1), " (", 
                                 round(prop.table(table(dat.ACTION$SEX))[1]*100,1), "%)"),
                          paste0(sum(dat.Avera$SEX==1), " (", 
                                 round(prop.table(table(dat.Avera$SEX))[1]*100,1), "%)"),
                          paste0(sum(c(dat.pilot$SEX,dat.ACTION$SEX,dat.Avera$SEX)==1), " (", 
                                 round(prop.table(table(c(dat.pilot$SEX,dat.ACTION$SEX,dat.Avera$SEX)))[1]*100,1), "%)"))
# add N twins per dataset
overviewmaxN$N_twin <- c(sum(dat.pilot$multiple_type==2,na.rm = T),
                         sum(dat.ACTION$multiple_type=="twin",na.rm = T),
                         sum(dat.Avera$multiple_type=="twin",na.rm = T),
                         sum(dat.pilot$multiple_type==2,na.rm = T)+sum(dat.ACTION$multiple_type=="twin",na.rm = T)+
                           sum(dat.Avera$multiple_type=="twin",na.rm = T))
# add N complete twin pairs
overviewmaxN$N_twin_cpair <- c(sum(as.data.frame(table(dat.pilot$FamilyNumber))$Freq==2),
                               sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+
                                 sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4),
                               sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="twin")]))$Freq==2),
                               sum(as.data.frame(table(dat.pilot$FamilyNumber))$Freq==2)+
                                 sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+
                                 sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4)+
                                 sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==2)+sum(as.data.frame(table(dat.ACTION$FamilyNumber))$Freq==4))

# Add age information per dataset
overviewmaxN$twin_M_SD_range_Age <- c(paste0(round(mean(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],na.rm = T),1)," (", 
                                        round(sd(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],na.rm = T),1),")", " [",
                                        round(min(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],na.rm = T),1), " - ", 
                                        round(max(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)], na.rm = T),1), "]"),
                                 paste0(round(mean(dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],na.rm = T),1)," (", 
                                        round(sd(dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],na.rm = T),1),")", " [",
                                        round(min(dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],na.rm = T),1), " - ", 
                                        round(max(dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")], na.rm = T),1), "]"),
                                 paste0(round(mean(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")],na.rm = T),1)," (", 
                                        round(sd(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")],na.rm = T),1),")", " [",
                                        round(min(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")],na.rm = T),1), " - ", 
                                        round(max(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")], na.rm = T),1), "]"),
                                 paste0(round(mean(c(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],
                                                     dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],
                                                     dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")]),na.rm = T),1)," (", 
                                        round(sd(c(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],
                                                   dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],
                                                   dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")]),na.rm = T),1),")", " [",
                                        round(min(c(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],
                                                    dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],
                                                    dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")]),na.rm = T),1), " - ", 
                                        round(max(c(dat.pilot$age_at_collection[which(dat.pilot$multiple_type==2)],
                                                    dat.ACTION$age_at_collection[which(dat.ACTION$multiple_type=="twin")],
                                                    dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin")]), na.rm = T),1), "]"))

# add sex information per data set 
overviewmaxN$twin_N_males <- c(paste0(sum(dat.pilot$SEX[which(dat.pilot$multiple_type==2)]==1), " (", 
                                 round(prop.table(table(dat.pilot$SEX[which(dat.pilot$multiple_type==2)]))[1]*100,1), "%)"),
                          paste0(sum(dat.ACTION$SEX[which(dat.ACTION$multiple_type=="twin")]==1), " (", 
                                 round(prop.table(table(dat.ACTION$SEX[which(dat.ACTION$multiple_type=="twin")]))[1]*100,1), "%)"),
                          paste0(sum(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin")]==1), " (", 
                                 round(prop.table(table(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin")]))[1]*100,1), "%)"),
                          paste0(sum(c(dat.pilot$SEX[which(dat.pilot$multiple_type==2)],dat.ACTION$SEX[which(dat.ACTION$multiple_type=="twin")],
                                       dat.Avera$SEX[which(dat.Avera$multiple_type=="twin")])==1), " (", 
                                 round(prop.table(table(c(dat.pilot$SEX[which(dat.pilot$multiple_type==2)],
                                                          dat.ACTION$SEX[which(dat.ACTION$multiple_type=="twin")],
                                                          dat.Avera$SEX[which(dat.Avera$multiple_type=="twin")])))[1]*100,1), "%)"))

# add N MZ per dataset
overviewmaxN$N_MZ <- c(sum(dat.pilot$twzyg==1|dat.pilot$twzyg==3),
                        sum(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF"),
                        sum(dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="MZM"|
                              dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="MZF"),
                        sum(dat.pilot$twzyg==1|dat.pilot$twzyg==3)+
                        sum(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")+
                        sum(dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="MZM"|
                              dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="MZF"))
# add N MZ complete pairs
overviewmaxN$N_MZ_cpair <- c(sum(as.data.frame(table(dat.pilot$FamilyNumber[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)]))$Freq==2),
                       sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")]))$Freq==2)+
                                           sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="MZM"|
                                                                                                   dat.ACTION$twzyg=="MZF")]))$Freq==4),
                       sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="twin" & 
                                                                              (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]))$Freq==2),
                       sum(as.data.frame(table(dat.pilot$FamilyNumber[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)]))$Freq==2)+
                       sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")]))$Freq==2)+
                         sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="MZM"|
                                                                                 dat.ACTION$twzyg=="MZF")]))$Freq==4)+
                       sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="twin" & 
                                                                              (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]))$Freq==2))
# add age information
overviewmaxN$MZ_M_SD_range_Age <- c(paste0(round(mean(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],na.rm = T),1)," (", 
                                           round(sd(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],na.rm = T),1),")", " [",
                                           round(min(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],na.rm = T),1), " - ", 
                                           round(max(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)], na.rm = T),1), "]"),
                                    paste0(round(mean(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],na.rm = T),1)," (", 
                                           round(sd(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],na.rm = T),1),")", " [",
                                           round(min(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],na.rm = T),1), " - ", 
                                           round(max(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")], na.rm = T),1), "]"),
                                    paste0(round(mean(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                          (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))],na.rm = T),1)," (", 
                                           round(sd(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                        (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))],na.rm = T),1),")", " [",
                                           round(min(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                         (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))],na.rm = T),1), " - ", 
                                           round(max(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                         (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))], na.rm = T),1), "]"),
                                    paste0(round(mean(c(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                                        dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                                        dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                            (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]),na.rm = T),1)," (", 
                                           round(sd(c(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                                      dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                                      dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                          (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]),na.rm = T),1),")", " [",
                                           round(min(c(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                                       dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                                       dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                           (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]),na.rm = T),1), " - ", 
                                           round(max(c(dat.pilot$age_at_collection[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                                       dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                                       dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                           (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]), na.rm = T),1), "]"))

# add sex information per data set 
overviewmaxN$MZ_N_males <- c(paste0(sum(dat.pilot$SEX[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)]==1), " (", 
                                    round(prop.table(table(dat.pilot$SEX[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)]))[1]*100,1), "%)"),
                             paste0(sum(dat.ACTION$SEX[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")]==1), " (", 
                                    round(prop.table(table(dat.ACTION$SEX[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")]))[1]*100,1), "%)"),
                             paste0(sum(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                              (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]==1), " (", 
                                    round(prop.table(table(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                                 (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))]))[1]*100,1), "%)"),
                             paste0(sum(c(dat.pilot$SEX[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                          dat.ACTION$SEX[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                          dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))])==1), " (", 
                                    round(prop.table(table(c(dat.pilot$SEX[which(dat.pilot$twzyg==1|dat.pilot$twzyg==3)],
                                                             dat.ACTION$SEX[which(dat.ACTION$twzyg=="MZM"|dat.ACTION$twzyg=="MZF")],
                                                             dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                                   (dat.Avera$twzyg=="MZM"|dat.Avera$twzyg=="MZF"))])))[1]*100,1), "%)"))


# add N DZ per dataset
overviewmaxN$N_DZ <- c(0,
                       sum(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm"),
                       sum(dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DZM"|
                             dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DZF"|
                             dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DOSfm"|
                             dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DOSmf"),
                       0+sum(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")+
                         sum(dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DZM"|
                               dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DZF"|
                               dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DOSfm"|
                               dat.Avera$twzyg[which(dat.Avera$multiple_type=="twin")]=="DOSmf"))
# add N DZ complete pairs
overviewmaxN$N_DZ_cpair <- c(0,
                             sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                     dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")]))$Freq==2),
                             sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="twin" & 
                                                                                    (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                       dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]))$Freq==2),
                             0+sum(as.data.frame(table(dat.ACTION$FamilyNumber[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                       dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")]))$Freq==2)+
                             sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="twin" & 
                                                                                    (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                       dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]))$Freq==2))

# add age information
overviewmaxN$DZ_M_SD_range_Age <- c(paste0("NA"),
                                    paste0(round(mean(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                           dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],na.rm = T),1)," (", 
                                           round(sd(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                         dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],na.rm = T),1),")", " [",
                                           round(min(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                          dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],na.rm = T),1), " - ", 
                                           round(max(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                          dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")], na.rm = T),1), "]"),
                                    paste0(round(mean(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                          (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                             dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))],na.rm = T),1)," (", 
                                           round(sd(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                        (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                           dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))],na.rm = T),1),")", " [",
                                           round(min(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                         (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                            dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))],na.rm = T),1), " - ", 
                                           round(max(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                         (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                            dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))], na.rm = T),1), "]"),
                                    paste0(round(mean(c(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                             dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                                        dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                            (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                               dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]),na.rm = T),1)," (", 
                                           round(sd(c(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                           dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                                      dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                          (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                             dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]),na.rm = T),1),")", " [",
                                           round(min(c(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                            dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                                       dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                           (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                              dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]),na.rm = T),1), " - ", 
                                           round(max(c(dat.ACTION$age_at_collection[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                            dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                                       dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="twin" & 
                                                                                           (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                              dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]), na.rm = T),1), "]"))

# add sex information per data set 
overviewmaxN$DZ_N_males <- c(paste0("NA"),
                             paste0(sum(dat.ACTION$SEX[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                               dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")]==1), " (", 
                                    round(prop.table(table(dat.ACTION$SEX[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                  dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")]))[1]*100,1), "%)"),
                             paste0(sum(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                              (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                 dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]==1), " (", 
                                    round(prop.table(table(dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                                 (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                    dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))]))[1]*100,1), "%)"),
                             paste0(sum(c(dat.ACTION$SEX[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                 dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                          dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                   dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))])==1), " (", 
                                    round(prop.table(table(c(dat.ACTION$SEX[which(dat.ACTION$twzyg=="DZM"|dat.ACTION$twzyg=="DZF"|
                                                                                    dat.ACTION$twzyg=="DOSmf"|dat.ACTION$twzyg=="DOSfm")],
                                                             dat.Avera$SEX[which(dat.Avera$multiple_type=="twin" & 
                                                                                   (dat.Avera$twzyg=="DZM"|dat.Avera$twzyg=="DZF"|
                                                                                      dat.Avera$twzyg=="DOSmf"|dat.Avera$twzyg=="DOSfm"))])))[1]*100,1), "%)"))

# N triplets
overviewmaxN$N_triplets <- c(0,0,sum(dat.Avera$multiple_type=="triplet",na.rm = T),
                             0+0+sum(dat.Avera$multiple_type=="triplet",na.rm = T))
# add N complete triplet sets
overviewmaxN$N_triplets_cset <- c(0,0,sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="triplet")]))$Freq==3),
                                  0+0+sum(as.data.frame(table(dat.Avera$FamilyNumber[which(dat.Avera$multiple_type=="triplet")]))$Freq==3))
# Add age information per dataset
overviewmaxN$trip_M_SD_range_Age <- c(paste0("NA"),
                                      paste0("NA"),
                                      paste0(round(mean(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")],na.rm = T),1)," (", 
                                             round(sd(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")],na.rm = T),1),")", " [",
                                             round(min(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")],na.rm = T),1), " - ", 
                                             round(max(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")], na.rm = T),1), "]"),
                                      paste0(round(mean(c(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")]),na.rm = T),1)," (", 
                                             round(sd(c(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")]),na.rm = T),1),")", " [",
                                             round(min(c(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")]),na.rm = T),1), " - ", 
                                             round(max(c(dat.Avera$age_at_collection[which(dat.Avera$multiple_type=="triplet")]), na.rm = T),1), "]"))

# add sex information per data set 
overviewmaxN$trip_N_males <- c(paste0("NA"),
                               paste0("NA"),
                               paste0(sum(dat.Avera$SEX[which(dat.Avera$multiple_type=="triplet")]==1), " (", 
                                      round(prop.table(table(dat.Avera$SEX[which(dat.Avera$multiple_type=="triplet")]))[1]*100,1), "%)"),
                               paste0(sum(c(dat.Avera$SEX[which(dat.Avera$multiple_type=="triplet")])==1), " (", 
                                      round(prop.table(table(c(dat.Avera$SEX[which(dat.Avera$multiple_type=="triplet")])))[1]*100,1), "%)"))

# N siblings
overviewmaxN$N_siblings <- c(0,0,sum(dat.Avera$role=="sibling"),
                             0+0+sum(dat.Avera$role=="sibling"))
# Add age information per dataset
overviewmaxN$sib_M_SD_range_Age <- c(paste0("NA"),
                                      paste0("NA"),
                                      paste0(round(mean(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")],na.rm = T),1)," (", 
                                             round(sd(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")],na.rm = T),1),")", " [",
                                             round(min(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")],na.rm = T),1), " - ", 
                                             round(max(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")], na.rm = T),1), "]"),
                                      paste0(round(mean(c(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")]),na.rm = T),1)," (", 
                                             round(sd(c(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")]),na.rm = T),1),")", " [",
                                             round(min(c(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")]),na.rm = T),1), " - ", 
                                             round(max(c(dat.Avera$age_at_collection[which(dat.Avera$role=="sibling")]), na.rm = T),1), "]"))

# add sex information per data set 
overviewmaxN$sib_N_males <- c(paste0("NA"),
                               paste0("NA"),
                               paste0(sum(dat.Avera$SEX[which(dat.Avera$role=="sibling")]==1), " (", 
                                      round(prop.table(table(dat.Avera$SEX[which(dat.Avera$role=="sibling")]))[1]*100,1), "%)"),
                               paste0(sum(c(dat.Avera$SEX[which(dat.Avera$role=="sibling")])==1), " (", 
                                      round(prop.table(table(c(dat.Avera$SEX[which(dat.Avera$role=="sibling")])))[1]*100,1), "%)"))


#write to file
write.table(overviewmaxN, file=paste0("data/intermediate/", as.character(Sys.Date()), "_EPIC_buccal_Descriptives.txt"), row.names=F, col.names = T, sep="\t",quote = F)
