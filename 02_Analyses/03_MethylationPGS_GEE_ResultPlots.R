################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: Create annotated table with significant CpGs, Manhattan plots and Main
# Figure 2 (overview of significantly associated CpGs)
#
# data: 1. results mega-analysis all traits (.txt) 2. EPIC annotation file
# (.RData)
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
#install.packages("data.table")
library(data.table)
#install.packages("qqman")
library(qqman)
#install.packages("qvalue")
library(qvalue)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("tidyverse")
library(tidyverse)


################################################################################
#
# Load & prepare data
#
################################################################################

#load annotation files 
epicanno <- get(load("/home/fhagenbeek/data/MethylationPGS/data/anno_epic_072017.RData"))
# remove X and Y 
epicanno <- epicanno[-which(epicanno$chromosome=="chrX" | epicanno$chromosome=="chrY"),]
# remove chr from CHR and make numeric
epicanno$chromosome <- gsub("chr","",epicanno$chromosome)
epicanno$chromosome <- as.numeric(epicanno$chromosome)

## Read in results overlapping CpGs all three datasets and merge annotation ##
# read in mega-analyses results
mega <- fread("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/2024-01-05_EPIC_Buccal_megaanalysis_GEE_results_alltraits.txt", data.table = F)
mega <- merge(mega,epicanno[,c("name","chromosome","position")], by.x = "trait", by.y = "name", all.x = T)

# remove CpGs that didn't converge 
mega <- mega[-which(mega$error==104),]

write.table(mega, file=paste0("/home/fhagenbeek/data/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                    "_EPIC_Buccal_megaanalysis_GEE_results_alltraits_noerror.txt"), row.names=F, col.names = T, sep="\t",quote = T) 


################################################################################
#
# Create table with number of significant hits
# set genome-significant threshold to bonferonni corrected for number of
# independent CpGs (0.05/1,850) = 2.702703e-05 and calculate FDR-significance
#
################################################################################

# Create table with number of significant hits
signtable <- as.data.frame(matrix(NA,ncol=0,nrow=6*3))
signtable$trait <- c("BMI - offspring", "BMI - maternal", "BMI - paternal",
                     "Height - offspring", "Height - maternal", "Height - paternal",
                     "Schizophrenia - offspring", "Schizophrenia - maternal", "Schizophrenia - paternal",
                     "Smoking Initiation - offspring", "Smoking Initiation - maternal", "Smoking Initiation - paternal",
                     "Educational Attainment - offspring", "Educational Attainment - maternal", "Educational Attainment - paternal",
                     "Social Deprivation - offspring", "Social Deprivation - maternal", "Social Deprivation - paternal")
signtable$mega_N_bonfsign <- c(sum(mega$BMI_p<=2.702703e-05),sum(mega$Moth_BMI_p<=2.702703e-05),sum(mega$Fath_BMI_p<=2.702703e-05),
                               sum(mega$Height_p<=2.702703e-05),sum(mega$Moth_Height_p<=2.702703e-05),sum(mega$Fath_Height_p<=2.702703e-05),
                               sum(mega$Schizophrenia_p<=2.702703e-05),sum(mega$Moth_Schizophrenia_p<=2.702703e-05),sum(mega$Fath_Schizophrenia_p<=2.702703e-05),
                               sum(mega$Smoking_SmokingInitiation_p<=2.702703e-05),sum(mega$Moth_Smoking_SmokingInitiation_p<=2.702703e-05),sum(mega$Fath_Smoking_SmokingInitiation_p<=2.702703e-05),
                               sum(mega$EducationalAttainment_p<=2.702703e-05),sum(mega$Moth_EducationalAttainment_p<=2.702703e-05),sum(mega$Fath_EducationalAttainment_p<=2.702703e-05),
                               sum(mega$SocialDeprivation_p<=2.702703e-05),sum(mega$Moth_SocialDeprivation_p<=2.702703e-05),sum(mega$Fath_SocialDeprivation_p<=2.702703e-05))
signtable$mega_N_FDRsign <- c(sum(qvalue(mega$BMI_p)$qvalue < 0.05),sum(qvalue(mega$Moth_BMI_p)$qvalue < 0.05),sum(qvalue(mega$Fath_BMI_p)$qvalue < 0.05),
                              sum(qvalue(mega$Height_p)$qvalue < 0.05),sum(qvalue(mega$Moth_Height_p)$qvalue < 0.05),sum(qvalue(mega$Fath_Height_p)$qvalue < 0.05),
                              sum(qvalue(mega$Schizophrenia_p)$qvalue < 0.05),sum(qvalue(mega$Moth_Schizophrenia_p)$qvalue < 0.05),sum(qvalue(mega$Fath_Schizophrenia_p)$qvalue < 0.05),
                              sum(qvalue(mega$Smoking_SmokingInitiation_p)$qvalue < 0.05),sum(qvalue(mega$Moth_Smoking_SmokingInitiation_p)$qvalue < 0.05),sum(qvalue(mega$Fath_Smoking_SmokingInitiation_p)$qvalue < 0.05),
                              sum(qvalue(mega$EducationalAttainment_p)$qvalue < 0.05),sum(qvalue(mega$Moth_EducationalAttainment_p)$qvalue < 0.05),sum(qvalue(mega$Fath_EducationalAttainment_p)$qvalue < 0.05),
                              sum(qvalue(mega$SocialDeprivation_p)$qvalue < 0.05),sum(qvalue(mega$Moth_SocialDeprivation_p)$qvalue < 0.05),sum(qvalue(mega$Fath_SocialDeprivation_p)$qvalue < 0.05))

# write each of the results for each of the phenotypes to file
write.table(signtable, file=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                   "_EPIC_Buccal_mega_alltraits_results_significancetable.txt"), row.names=F, col.names = T, sep="\t",quote = T)

# Define significance threshold
significance_threshold <- 2.702703e-05

# Initialize an empty dataframe to store significant outcomes and predictors
significant_results <- data.frame(CpG = character(),
                                  PGS = character(),
                                  stringsAsFactors = FALSE)

# Loop through each row of the dataframe
for (i in 1:nrow(mega)) {
  # Extract p-values for the predictors
  p_values <- mega[i, c("Height_p","BMI_p","Smoking_SmokingInitiation_p","Schizophrenia_p","EducationalAttainment_p","SocialDeprivation_p",
                        "Moth_Height_p","Moth_BMI_p","Moth_Smoking_SmokingInitiation_p","Moth_Schizophrenia_p","Moth_EducationalAttainment_p",
                        "Moth_SocialDeprivation_p","Fath_Height_p","Fath_BMI_p","Fath_Smoking_SmokingInitiation_p","Fath_Schizophrenia_p",
                        "Fath_EducationalAttainment_p","Fath_SocialDeprivation_p")]
  
  # Check if any p-value is below the significance threshold
  if (any(p_values < significance_threshold)) {
    # Extract significant predictors
    significant_predictors <- names(p_values)[which(p_values < significance_threshold)]
    
    # Rename significant predictors
    significant_predictors <- gsub("^Moth_", "Maternal_", significant_predictors)
    significant_predictors <- gsub("^Fath_", "Paternal_", significant_predictors)
    significant_predictors <- gsub("_p$", "", significant_predictors)
    significant_predictors <- ifelse(!grepl("^Maternal_|^Paternal_", significant_predictors), paste0("Offspring_", significant_predictors), significant_predictors)
    
    # Store significant outcomes and predictors in the dataframe
    significant_results <- rbind(significant_results, data.frame(CpG = mega$trait[i],
                                                                 PGS = significant_predictors,
                                                                 stringsAsFactors = FALSE))
  }
}

# add annotation
significant_results <- merge(significant_results,epicanno[,c("name","chromosome","position","gene.symbol")], by.x = "CpG", by.y = "name", all.x = T)

# write each of the results for each of the phenotypes to file
write.table(significant_results, file=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                   "_EPIC_Buccal_mega_alltraits_results_significancetable_annotated.txt"), row.names=F, col.names = T, sep="\t",quote = T)



################################################################################
#
# Create MH plots
# set genome-significant threshold to bonferonni corrected for number of
# independent CpGs (0.05/1,850) = 2.702703e-05
#
################################################################################

## for each trait plot the MH per trait for the offspring and paternal PGSs ##

# BMI
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                       "_GEEalltraits_results_MH_BMI_megaanalysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "BMI_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_BMI_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_BMI_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Paternal")
dev.off()

# EA
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                       "_GEEalltraits_results_MH_EA_megaanalysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "EducationalAttainment_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_EducationalAttainment_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_EducationalAttainment_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Paternal")

dev.off()

# Height
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                       "_GEEalltraits_results_MH_Height_megaanlysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "Height_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_Height_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_Height_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 6),
          main = "Paternal")
dev.off()

# SCZ
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                       "_GEEalltraits_results_MH_SCZ_megaanalysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "Schizophrenia_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 8),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_Schizophrenia_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 8),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_Schizophrenia_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 8),
          main = "Paternal")
dev.off()

# SMOKE
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                       "_GEEalltraits_results_MH_SMOKE_megaanalysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "Smoking_SmokingInitiation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 4),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_Smoking_SmokingInitiation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 4),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_Smoking_SmokingInitiation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 4),
          main = "Paternal")
dev.off()

# SOCDEP
tiff(filename = paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                       "_GEEalltraits_results_MH_SOCDEP_megaanalysis.tiff"),
     width=2600,height=3774,units='px',res=600,pointsize=6)
par(mfrow=c(3,1))
manhattan(mega,chr = "chromosome", bp = "position", p = "SocialDeprivation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Offspring")
manhattan(mega,chr = "chromosome", bp = "position", p = "Moth_SocialDeprivation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Maternal")
manhattan(mega,chr = "chromosome", bp = "position", p = "Fath_SocialDeprivation_p", snp = "trait",
          suggestiveline=FALSE, genomewideline = -log10(2.702703e-05),ylim = c(0, 5),
          main = "Paternal")
dev.off()


################################################################################
#
# Create forest plot style figure of significant CpGs
#
################################################################################

# create dataframe with any significant CpGs - maternal or paternal or offspring
mega_both <- mega[which(mega$Height_p<=2.702703e-05 | mega$BMI_p<=2.702703e-05 | mega$Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Schizophrenia_p<=2.702703e-05 | mega$EducationalAttainment_p<=2.702703e-05 | mega$SocialDeprivation_p<=2.702703e-05 | 
                          mega$Moth_Height_p<=2.702703e-05 | mega$Moth_BMI_p<=2.702703e-05 | mega$Moth_Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Moth_Schizophrenia_p<=2.702703e-05 | mega$Moth_EducationalAttainment_p<=2.702703e-05 | mega$Moth_SocialDeprivation_p<=2.702703e-05 |
                          mega$Fath_Height_p<=2.702703e-05 | mega$Fath_BMI_p<=2.702703e-05 | mega$Fath_Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Fath_Schizophrenia_p<=2.702703e-05 | mega$Fath_EducationalAttainment_p<=2.702703e-05 | mega$Fath_SocialDeprivation_p<=2.702703e-05),]

# reformat output of the models for plotting 
dfc <- data.frame(beta = unlist(mega_both[,c("Height_beta","BMI_beta","Smoking_SmokingInitiation_beta",
                                             "Schizophrenia_beta","EducationalAttainment_beta","SocialDeprivation_beta",
                                             "Moth_Height_beta","Moth_BMI_beta","Moth_Smoking_SmokingInitiation_beta",
                                                 "Moth_Schizophrenia_beta","Moth_EducationalAttainment_beta","Moth_SocialDeprivation_beta","Fath_Height_beta",
                                                 "Fath_BMI_beta","Fath_Smoking_SmokingInitiation_beta","Fath_Schizophrenia_beta","Fath_EducationalAttainment_beta",
                                                 "Fath_SocialDeprivation_beta")]),
                  se = unlist(mega_both[,c("Height_se","BMI_se","Smoking_SmokingInitiation_se",
                                           "Schizophrenia_se","EducationalAttainment_se","SocialDeprivation_se",
                                           "Moth_Height_se","Moth_BMI_se","Moth_Smoking_SmokingInitiation_se",
                                               "Moth_Schizophrenia_se","Moth_EducationalAttainment_se","Moth_SocialDeprivation_se","Fath_Height_se",
                                               "Fath_BMI_se","Fath_Smoking_SmokingInitiation_se","Fath_Schizophrenia_se","Fath_EducationalAttainment_se",
                                               "Fath_SocialDeprivation_se")]),
                  p = unlist(mega_both[,c("Height_p","BMI_p","Smoking_SmokingInitiation_p",
                                          "Schizophrenia_p","EducationalAttainment_p","SocialDeprivation_p",
                                          "Moth_Height_p","Moth_BMI_p","Moth_Smoking_SmokingInitiation_p",
                                              "Moth_Schizophrenia_p","Moth_EducationalAttainment_p","Moth_SocialDeprivation_p","Fath_Height_p",
                                              "Fath_BMI_p","Fath_Smoking_SmokingInitiation_p","Fath_Schizophrenia_p","Fath_EducationalAttainment_p",
                                              "Fath_SocialDeprivation_p")]),
                  trait = c(rep(mega_both$trait, 18)),
                  PGS = c(rep(c("Height","BMI","Smoking Initiation","Schizophrenia","Educational Attainment","Social Deprivation"),
                              each=nrow(mega_both))),
                  Parent = c(rep(c("Offspring","Maternal","Paternal"),
                                 each=6*nrow(mega_both))))
# Trait and parent levels as factor to plot EA in order of magnitude
dfc$trait <- factor(dfc$trait, levels = paste(unique(dfc$trait)))
dfc$Parent <- factor(dfc$Parent, levels = c("Maternal","Offspring","Paternal"))

# calculate lower and upper bound of 95% CI
dfc$lb <- dfc$beta-1.96*dfc$se
dfc$ub <- dfc$beta+1.96*dfc$se

# set all non-signicant beta's and se's to NA
dfc$beta[which(dfc$p>2.702703e-05)] <- NA
dfc$se[which(dfc$p>2.702703e-05)] <- NA
dfc$lb[which(dfc$p>2.702703e-05)] <- NA
dfc$ub[which(dfc$p>2.702703e-05)] <- NA

# color blind friendly pallette
cbPalette <- c("#009E73","#56B4E9","#F0E442","#CC79A7","#000000","#E69F00","#0072B2","#D55E00")

# plot 
betaplot3 <- ggplot(dfc, aes(x = beta, y = trait, 
                             xmin=lb, xmax=ub,
                             group = PGS, color = PGS, shape = Parent)) + 
  geom_point(position = position_dodge(1),size=3.5) + 
  geom_linerange(position = position_dodge(1),size=1) +
  theme_classic() + scale_color_manual(values = cbPalette) + 
  geom_vline(xintercept = 0) +
  theme(legend.position = "bottom",legend.title = element_blank(),text=element_text(size=14)) + 
  labs(x="Beta 95% CI", y=NULL)


# save figure as pdfa.
ggsave(filename=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/",as.character(Sys.Date()),
                       "_GEEalltraits_results_FIG_BetaCoefficients_both.pdf"),
       plot = betaplot3, device = pdf, dpi = 300, width = 350, height = 250, units = "mm")
