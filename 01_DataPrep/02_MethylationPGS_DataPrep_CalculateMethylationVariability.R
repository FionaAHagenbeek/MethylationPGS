################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: Calculate variability (SD) of the residualized DNA methylation probes
# after combining all datasets and removing all outliers (>3 IQR)
#
# data: 1. imputed NTR ACTION Methylation data all batches  + covariates
# (.RData) 2. imputed ACTION pilot data  + covariates (.Rdata) 3. imputed Avera
# data + covariates (Rdata) 4. EPIC annotation files
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
#install.packages("remotes")
# library(remotes)
# install_github("perishky/ewaff")
library(ewaff)
#install.packages("doMC")
library(doMC)
registerDoMC(15)  #number of cpu's to use


################################################################################
#
# Load data
#
################################################################################

# DNA methylation covariate files NTR ACTION (.Rdata)
covariates.ACTION <- get(load("/home/fhagenbeek/data/MethylationPGS/data/Covariates_NTR_Curium_V3.RData"))
covariates.pilot <- get(load("/home/fhagenbeek/data/MethylationPGS/data/Covariates_ACTION_EPICarray_pilot_DSR_3454_21032023.RData"))
covariates.Avera <- get(load("/home/fhagenbeek/data/MethylationPGS/data/Covariates_EPIC_Avera_Group2_Buccal_NTR_DSR_3454_09032023.RData")) 
rm(covariates) # remove duplicate covariates 
#rename Avera covariates
names(covariates.Avera) <- c("fisnumber","IdatName","Sample_Well","Sample_Plate","Sentrix_ID","Sentrix_Position","duplicates",
                             "PedigreeNumber","FamilyNumber","Extension","role","multiple_type","birthorder","Sex","trans",
                             "confirmed","twzyg","famzyg","mult_ext","ANTR","YNTR","zyg12","zyg13","zyg23","zygdna12","zygdna13",
                             "zygdna23","buccalsubgroup","Epi","Fib","B","NK","CD4T","CD8T","Mono","Neutro","Eosino","Array_rownum",
                             "age_at_collection")

# Buccal cell counts - NTR ACTION
load("/home/fhagenbeek/data/MethylationPGS/data/buccalestimates_EPIDISH_FN5pcs_BIOS.Robj")
EPIdish_WBC <- as.data.frame(EPIdish_WBC) #as data frame

#read in imputed methylation betas 
betas.pilot <- get(load("/data/PUBLIC/Methylation/ACTION_EPICarray_pilot/ACTION.EPIC.betas_pilot_nomissings.RData")) #load ACTION pilot
betas.ACTION <- get(load("/data/PUBLIC/Methylation/ACTION_EPICarray/ACTION.EPIC.betas_NTR_nomissings.RData")) #load betas NTR ACTION
betas.Avera <- get(load("/data/PUBLIC/Methylation/Avera_2022_EPICarray/Group2_Buccal/Betas_Avera2022_Group2_Buccal_012023_nomissings.RData")) #load betas Avera Group 2
rm(Beta.FunNorm) # remove duplicate betas 

#load annotation files 
crossreactive <- get(load("/home/fhagenbeek/data/MethylationPGS/data/crossreactive47.RData"))
DNApolymorphismtarget <- get(load("/home/fhagenbeek/data/MethylationPGS/data/DNApolymorphismtargetCpG_EUR01.RData"))
DNApolymorphismSBE <- get(load("/home/fhagenbeek/data/MethylationPGS/data/DNApolymorphismSBE_EUR01.RData"))
epicanno <- get(load("/home/fhagenbeek/data/MethylationPGS/data/anno_epic_072017.RData"))

################################################################################
#
# probe filtering - remove cross-reactive probes, probes in a SNP and XY-probes
# and remove outliers (>3 IQR)
#
################################################################################

## keep overlapping CpG sites ## 

# pilot
OL <- intersect(rownames(betas.pilot),epicanno$name) #obtain overlap
rownames(epicanno) <- epicanno$name #give the annotation row names the probe names
epicanno <- epicanno[OL,] #only retain overlapping CpGs in annotation file
betas.pilot <- betas.pilot[OL,] #only retain overlapping CpGs in the file with DNA methylation B-values

# NTR ACTION
OL <- intersect(rownames(betas.ACTION),epicanno$name) #obtain overlap
rownames(epicanno) <- epicanno$name #give the annotation row names the probe names
epicanno <- epicanno[OL,] #only retain overlapping CpGs in annotation file
betas.ACTION <- betas.ACTION[OL,] #only retain overlapping CpGs in the file with DNA methylation B-values

# Avera
OL <- intersect(rownames(betas.Avera),epicanno$name) #obtain overlap
rownames(epicanno) <- epicanno$name #give the annotation row names the probe names
epicanno <- epicanno[OL,] #only retain overlapping CpGs in annotation file
betas.Avera <- betas.Avera[OL,] #only retain overlapping CpGs in the file with DNA methylation B-values

## remove cross-reactive probes and probes overlapping SNPs ##

# pilot
rm <- c(crossreactive, DNApolymorphismtarget,DNApolymorphismSBE) #create a vector with all CpGs to be removed
betas.pilot <- betas.pilot[which(!rownames(betas.pilot) %in% rm),] #only retain CpGs not in the to-be-removed list

# NTR ACTION
rm <- c(crossreactive, DNApolymorphismtarget,DNApolymorphismSBE) #create a vector with all CpGs to be removed
betas.ACTION <- betas.ACTION[which(!rownames(betas.ACTION) %in% rm),] #only retain CpGs not in the to-be-removed list

# Avera
rm <- c(crossreactive, DNApolymorphismtarget,DNApolymorphismSBE) #create a vector with all CpGs to be removed
betas.Avera <- betas.Avera[which(!rownames(betas.Avera) %in% rm),] #only retain CpGs not in the to-be-removed list

## remove XY chromosomes ##

# pilot
epicanno <- epicanno[which(epicanno$chromosome=="chrX" | epicanno$chromosome=="chrY"),] #only retain X and Y CHRs in the annotation file
if(length(rownames(epicanno))>0) {
  betas.pilot <- betas.pilot[which(!rownames(betas.pilot) %in% rownames(epicanno)),] #only retain DNA methylation B-values not in X or Y CHR
}

# NTR ACTION
epicanno <- epicanno[which(epicanno$chromosome=="chrX" | epicanno$chromosome=="chrY"),] #only retain X and Y CHRs in the annotation file
if(length(rownames(epicanno))>0) {
  betas.ACTION <- betas.ACTION[which(!rownames(betas.ACTION) %in% rownames(epicanno)),] #only retain DNA methylation B-values not in X or Y CHR
}

# Avera
epicanno <- epicanno[which(epicanno$chromosome=="chrX" | epicanno$chromosome=="chrY"),] #only retain X and Y CHRs in the annotation file
if(length(rownames(epicanno))>0) {
  betas.Avera <- betas.Avera[which(!rownames(betas.Avera) %in% rownames(epicanno)),] #only retain DNA methylation B-values not in X or Y CHR
}

## remove outliers ##
methylation.pilot <- ewaff.handle.outliers(betas.pilot, method="iqr")[[1]] # pilot
methylation.ACTION <- ewaff.handle.outliers(betas.ACTION, method="iqr")[[1]] # NTR ACTION
methylation.Avera <- ewaff.handle.outliers(betas.Avera, method="iqr")[[1]] # Avera


################################################################################
#
# combine datasets 
#
################################################################################

# NTR ACTION merge buccal cell counts to covariates; only retain individuals
# with DNA methylation & buccal cell counts
EPIdish_WBC <- EPIdish_WBC[intersect(rownames(EPIdish_WBC),covariates.ACTION$IdatName),]
# add column with ID
EPIdish_WBC$IdatName <- rownames(EPIdish_WBC)
# merge to covariates
covariates.ACTION <- merge(covariates.ACTION, EPIdish_WBC, by = "IdatName")

# combine covariates
covariates <- rbind(covariates.pilot[,c("IdatName","age_at_collection","Sex","Sample_Plate","Array_rownum","Epi","NK")],
                          covariates.ACTION[,c("IdatName","age_at_collection","Sex","Sample_Plate","Array_rownum","Epi","NK")],
                          covariates.Avera[,c("IdatName","age_at_collection","Sex","Sample_Plate","Array_rownum","Epi","NK")])

# remove betas not present in all datasets
OL <- intersect(rownames(methylation.pilot),rownames(methylation.ACTION)) #obtain overlap
OL2 <- intersect(OL,rownames(methylation.Avera)) #obtain overlap
methylation.pilot <- methylation.pilot[OL2,]
methylation.ACTION <- methylation.ACTION[OL2,]
methylation.Avera <- methylation.Avera[OL2,]

# combine betas
methylation <- cbind(methylation.pilot,methylation.ACTION,methylation.Avera)

# how many CpGs across all three datasets?
nrow(methylation)


################################################################################
#
# Calculate DNA methylation probe variability
#
################################################################################

# only retain individuals with DNA methylation & covariate data
covariates <- covariates[which(covariates$IdatName %in% colnames(methylation)),]

# select variables
age <- covariates$age_at_collection
sex <- as.numeric(factor(covariates$Sex, levels = c("M","F")))-1
Sample_Plate <- as.factor(covariates$Sample_Plate)
Array_rownum <- covariates$Array_rownum
Epi <- covariates$Epi
NK <- covariates$NK

#create 1 dataframe with betas and covariates
batch <- data.frame(as.data.frame(t(methylation)),age,sex,Sample_Plate,Array_rownum,Epi,NK)

# create empty  matrix to store CpG name + SD raw value 
CpG <- as.data.frame(matrix(NA,ncol=0,nrow=nrow(methylation))) # create empty dataframe
CpG$CpG <- rownames(methylation) # add CpG name column 

# function to calculate the variability (SDs) of the residualized DNA methylation beta values
var.resid = function(dat,probe) { 
  out=sd(resid(lm(as.formula(paste0(probe," ~ age + sex + Sample_Plate + Array_rownum + Epi + NK")), data = dat, na.action=na.exclude)),na.rm = T)
  return(out)
}

#run function to calculate the variability (SD) of each CpG
var.resid.CpG  <- foreach(i=1:nrow(methylation)) %dopar% {
  var.resid(dat = batch, probe = CpG$CpG[i])
}

#add residual sd column
CpG$ressd <- c(unlist(var.resid.CpG))


#write output to file
save(CpG, file=paste("/home/fhagenbeek/data/MethylationPGS/data/intermediate/", as.character(Sys.Date()), "_EPIC_Buccal_SDprobes_alldats.RData",sep=""))
