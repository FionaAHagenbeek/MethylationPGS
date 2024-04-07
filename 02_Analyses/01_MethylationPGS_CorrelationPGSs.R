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
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("corrplot")
library(corrplot)


################################################################################
#
# Load data
#
################################################################################

# residualized and imputed DNA methylation betas and PGSs (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/final/2024-01-04_ImputedResidualDNAmPGSandCovs_alldats.RData")


################################################################################
#
# Correlation polygenic scores
#
################################################################################

# calculate correlations
PGScordat <- cor(combi2[,c("P_inf_SCORE_Height_MRG16_LDp1","P_inf_SCORE_BMI_MRG16_LDp1",                           
                        "P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","P_inf_SCORE_Schizophrenia_MRG16_LDp1",                
                        "P_inf_SCORE_EducationalAttainment_MRG16_LDp1","P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                        "Moth_P_inf_SCORE_Height_MRG16_LDp1","Moth_P_inf_SCORE_BMI_MRG16_LDp1",
                        "Moth_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Moth_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                        "Moth_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Moth_P_inf_SCORE_SocialDeprivation_MRG16_LDp1",
                        "Fath_P_inf_SCORE_Height_MRG16_LDp1","Fath_P_inf_SCORE_BMI_MRG16_LDp1",
                        "Fath_P_inf_SCORE_Smoking_SmokingInitiation_MRG16_LDp1","Fath_P_inf_SCORE_Schizophrenia_MRG16_LDp1",
                        "Fath_P_inf_SCORE_EducationalAttainment_MRG16_LDp1","Fath_P_inf_SCORE_SocialDeprivation_MRG16_LDp1")])
# rename columns and rows
colnames(PGScordat) <- c("Offspring Height","Offspring BMI","Offspring Smoking Initiation","Offspring Schizophrenia",
                      "Offspring Educational Attainment","Offspring Social Deprivation","Maternal Height","Maternal BMI",
                      "Maternal Smoking Initiation","Maternal Schizophrenia","Maternal Educational Attainment",
                      "Maternal Social Deprivation","Paternal Height","Paternal BMI","Paternal Smoking Initiation",
                      "Paternal Schizophrenia","Paternal Educational Attainment","Paternal Social Deprivation")
rownames(PGScordat) <- c("Offspring Height","Offspring BMI","Offspring Smoking Initiation","Offspring Schizophrenia",
                      "Offspring Educational Attainment","Offspring Social Deprivation","Maternal Height","Maternal BMI",
                      "Maternal Smoking Initiation","Maternal Schizophrenia","Maternal Educational Attainment",
                      "Maternal Social Deprivation","Paternal Height","Paternal BMI","Paternal Smoking Initiation",
                      "Paternal Schizophrenia","Paternal Educational Attainment","Paternal Social Deprivation")

# write correlations to file 
write.table(PGScordat, file=paste0("/home/fhagenbeek/data/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()), 
                                "_GEEalltraits_results_Correlations_PGSs.txt"), 
            row.names=F, col.names = T, sep="\t",quote = T)


################################################################################
#
# Figure correlation polygenic scores
#
################################################################################

# Order vector
order_vector <- c("Offspring Schizophrenia", "Maternal Schizophrenia", "Paternal Schizophrenia",
                  "Offspring Smoking Initiation", "Maternal Smoking Initiation", "Paternal Smoking Initiation",
                  "Offspring Educational Attainment", "Maternal Educational Attainment", "Paternal Educational Attainment",
                  "Offspring Social Deprivation", "Maternal Social Deprivation", "Paternal Social Deprivation",
                  "Offspring BMI", "Maternal BMI", "Paternal BMI",
                  "Offspring Height", "Maternal Height", "Paternal Height")

# function to reorder correlation matrix
reorder <- function(CorMat){
  CorMar <- CorMat[rev(order_vector), rev(order_vector)]
}

# functions to get upper and lower triangles
get_upper_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}
get_lower_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}

# reorder correlation matrix and get upper/lower triagle
CorMat <- reorder(PGScordat)
upper_tri <- get_upper_tri(CorMat)
lower_tri <- get_lower_tri(CorMat)
meltNum <- melt(lower_tri, na.rm = T)
meltColor <- melt(upper_tri, na.rm = T)

# plot
PGScor <- ggplot() +
  labs(x = NULL, y = NULL) +
  geom_tile(data = meltColor, 
            mapping = aes(Var2, Var1, 
                          fill = value)) +
  geom_text(data = meltNum,
            mapping = aes(Var2, Var1,
                          label = round(value, digit = 2))) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1,1), name = "Pearson\nCorrelation") +
  theme(axis.text.x = element_text(angle=90,hjust = 0),
        text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  coord_fixed()

# save plot
ggsave(filename=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/",as.character(Sys.Date()),
                       "_GEEalltraits_results_Correlations_PGS.pdf"),
       plot = PGScor, device = pdf, dpi = 600, width = 350, height = 350, units = "mm")


################################################################################
#
# Summarize correlation polygenic scores
#
################################################################################

# melt complete correlation set
meltcor <- melt(PGScordat)

## range, mean, and median correlations - same trait ##
# subset data
same <- meltcor[which((meltcor$Var1=="Offspring Schizophrenia" & (meltcor$Var2=="Maternal Schizophrenia" | meltcor$Var2=="Paternal Schizophrenia")) | 
                        (meltcor$Var1=="Offspring Smoking Initiation" & (meltcor$Var2=="Maternal Smoking Initiation" | meltcor$Var2=="Paternal Smoking Initiation")) | 
                        (meltcor$Var1=="Offspring Educational Attainment" & (meltcor$Var2=="Maternal Educational Attainment" | meltcor$Var2=="Paternal Educational Attainment")) | 
                        (meltcor$Var1=="Offspring Social Deprivation" & (meltcor$Var2=="Maternal Social Deprivation" | meltcor$Var2=="Paternal Social Deprivation")) | 
                        (meltcor$Var1=="Offspring BMI" & (meltcor$Var2=="Maternal BMI" | meltcor$Var2=="Paternal BMI")) | 
                        (meltcor$Var1=="Offspring Height" & (meltcor$Var2=="Maternal Height" | meltcor$Var2=="Paternal Height"))),]
# range, mean, and median correlations
round(range(same$value),2)
round(mean(same$value),2)
round(median(same$value),2)

## range, mean, and median correlations - same trait parents ##
# subset data
parents <- meltcor[which((meltcor$Var1=="Maternal Schizophrenia" & meltcor$Var2=="Paternal Schizophrenia") | 
                        (meltcor$Var1=="Maternal Smoking Initiation" & meltcor$Var2=="Paternal Smoking Initiation") | 
                        (meltcor$Var1=="Maternal Educational Attainment" & meltcor$Var2=="Paternal Educational Attainment") | 
                        (meltcor$Var1=="Maternal Social Deprivation" & meltcor$Var2=="Paternal Social Deprivation") | 
                        (meltcor$Var1=="Maternal BMI" & meltcor$Var2=="Paternal BMI") | 
                        (meltcor$Var1=="Maternal Height" & meltcor$Var2=="Paternal Height")),]
# range, mean, and median correlations
round(range(parents$value),2)
round(mean(parents$value),2)
round(median(parents$value),2)

# range, mean, and median correlations - different traits
# subset data
diff <- meltcor[-which((meltcor$Var1=="Maternal Schizophrenia" & meltcor$Var2=="Paternal Schizophrenia") | 
                           (meltcor$Var1=="Maternal Smoking Initiation" & meltcor$Var2=="Paternal Smoking Initiation") | 
                           (meltcor$Var1=="Maternal Educational Attainment" & meltcor$Var2=="Paternal Educational Attainment") | 
                           (meltcor$Var1=="Maternal Social Deprivation" & meltcor$Var2=="Paternal Social Deprivation") | 
                           (meltcor$Var1=="Maternal BMI" & meltcor$Var2=="Paternal BMI") | 
                           (meltcor$Var1=="Maternal Height" & meltcor$Var2=="Paternal Height") |
                         (meltcor$Var2=="Maternal Schizophrenia" & meltcor$Var1=="Paternal Schizophrenia") | 
                         (meltcor$Var2=="Maternal Smoking Initiation" & meltcor$Var1=="Paternal Smoking Initiation") | 
                         (meltcor$Var2=="Maternal Educational Attainment" & meltcor$Var1=="Paternal Educational Attainment") | 
                         (meltcor$Var2=="Maternal Social Deprivation" & meltcor$Var1=="Paternal Social Deprivation") | 
                         (meltcor$Var2=="Maternal BMI" & meltcor$Var1=="Paternal BMI") | 
                         (meltcor$Var2=="Maternal Height" & meltcor$Var1=="Paternal Height") |
                        (meltcor$Var1=="Offspring Schizophrenia" & (meltcor$Var2=="Maternal Schizophrenia" | meltcor$Var2=="Paternal Schizophrenia")) | 
                        (meltcor$Var1=="Offspring Smoking Initiation" & (meltcor$Var2=="Maternal Smoking Initiation" | meltcor$Var2=="Paternal Smoking Initiation")) | 
                        (meltcor$Var1=="Offspring Educational Attainment" & (meltcor$Var2=="Maternal Educational Attainment" | meltcor$Var2=="Paternal Educational Attainment")) | 
                        (meltcor$Var1=="Offspring Social Deprivation" & (meltcor$Var2=="Maternal Social Deprivation" | meltcor$Var2=="Paternal Social Deprivation")) | 
                        (meltcor$Var1=="Offspring BMI" & (meltcor$Var2=="Maternal BMI" | meltcor$Var2=="Paternal BMI")) | 
                        (meltcor$Var1=="Offspring Height" & (meltcor$Var2=="Maternal Height" | meltcor$Var2=="Paternal Height")) |
                         (meltcor$Var2=="Offspring Schizophrenia" & (meltcor$Var1=="Maternal Schizophrenia" | meltcor$Var1=="Paternal Schizophrenia")) | 
                         (meltcor$Var2=="Offspring Smoking Initiation" & (meltcor$Var1=="Maternal Smoking Initiation" | meltcor$Var1=="Paternal Smoking Initiation")) | 
                         (meltcor$Var2=="Offspring Educational Attainment" & (meltcor$Var1=="Maternal Educational Attainment" | meltcor$Var1=="Paternal Educational Attainment")) | 
                         (meltcor$Var2=="Offspring Social Deprivation" & (meltcor$Var1=="Maternal Social Deprivation" | meltcor$Var1=="Paternal Social Deprivation")) | 
                         (meltcor$Var2=="Offspring BMI" & (meltcor$Var1=="Maternal BMI" | meltcor$Var1=="Paternal BMI")) | 
                         (meltcor$Var2=="Offspring Height" & (meltcor$Var1=="Maternal Height" | meltcor$Var1=="Paternal Height")) |
                        meltcor$Var1==meltcor$Var2),]
# range, mean, and median correlations
round(range(diff$value),2)
round(mean(diff$value),3)
round(median(diff$value),2)
# strongest correlation
diff[which(abs(diff$value) %in% max(abs(diff$value))),]

