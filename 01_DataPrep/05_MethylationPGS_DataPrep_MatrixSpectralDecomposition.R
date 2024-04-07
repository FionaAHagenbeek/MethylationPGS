################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: matrix spectral decomposition of top 10% most variable DNA probes
#
# data: 1. NTR ACTION DNA methylation betas (imputed) (.RData) 2. top 10% most
# variable DNA Methylation probes (.RData)
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
#install.packages("dplyr")
library(dplyr)


################################################################################
#
# Load & prepare data
#
################################################################################

#read in imputed methylation betas 
betas.ACTION <- get(load("/data/PUBLIC/Methylation/ACTION_EPICarray/ACTION.EPIC.betas_NTR_nomissings.RData")) #load betas NTR ACTION

# top 10% most variable DNA methylation betas  (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/2024-01-03_EPIC_Buccal_Top10VariableCpGs_alldats.RData")

# retain only overlapping betas
OL.ACTION <- intersect(rownames(betas.ACTION),CpG.top10.anno$CpG)
betas.ACTION <- betas.ACTION[OL.ACTION,] #betas NTR ACTION

#create dataframes with only complete cases
completedat <- na.omit(t(betas.ACTION)) 

#create vectors with variable names 
CpGs <- colnames(completedat)


################################################################################
#
# Matrix Spectral decomposition (MEffLi Method)
#
# code from BBMRI Metabolomics workshop 
#
# "Here's a function for determining MEffLi. Note that for this method, you have
# to apply some operations on the data, ensuring that there are no missing
# values. Argument Columns are, in this case the list of metabolites present in
# the filtered data."
#
################################################################################

# Matrix Spectral Decomposition function
DetermineMEffLi <- function(Data,Columns){
  MeanCenteredAutoScaledData <- Data
  for(i in 1:length(Columns)){
    Entry <- Columns[i]
    Mean <- mean(Data[,Entry])
    Std <- sd(Data[,Entry])
    MeanCenteredAutoScaledData[,Entry] <- (MeanCenteredAutoScaledData[,Entry] - Mean)/Std
  }
  
  
  SVD <- svd(x=as.matrix(MeanCenteredAutoScaledData[,Columns]))
  EigenValues <- SVD$d^2 / dim(SVD$u)[1]
  M <- length(EigenValues)
  L <- M - 1
  Var <- var(EigenValues)
  MEff <- M*(1.0-(L*Var/M^2))
  IntEVals <- as.numeric(EigenValues>=1.0)
  NonIntEVals <- EigenValues - floor(EigenValues)
  MEffLi <- ceiling(sum(IntEVals) + sum(NonIntEVals))
  MEffLi
}

#running MSD
MEffLiC <- DetermineMEffLi(Data=completedat, Columns=CpGs) 
print (MEffLiC) #1850

#bonferonni correction 
print(0.05/MEffLiC) #2.702703e-05
