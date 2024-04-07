################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: compare significant CpGs to hits previous EWAS studies; compare mQTLs
# for significant CpGs with hits previous GWAS studies; and look-up mQTLs for
# significant CpGs in the buccal-brain overlap data
#
# data: 1. results mega-analysis all traits (.txt) 2. EPIC annotation file
# (.RData) 3. list of mQTLs associated with the significant CpGs (.csv) 4. GWAS
# summary statistics for schizophrenia, EA, and height (.tsv/.txt) 5. file with
# correlations between buccal and brain data
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
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("dplyr")
# install ggplot extension for creating VennDiagrams and require cowplot for combining plots
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(cowplot)

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

# create dataframe with any significant CpGs - maternal or paternal or offspring
mega_both <- mega[which(mega$Height_p<=2.702703e-05 | mega$BMI_p<=2.702703e-05 | mega$Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Schizophrenia_p<=2.702703e-05 | mega$EducationalAttainment_p<=2.702703e-05 | mega$SocialDeprivation_p<=2.702703e-05 | 
                          mega$Moth_Height_p<=2.702703e-05 | mega$Moth_BMI_p<=2.702703e-05 | mega$Moth_Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Moth_Schizophrenia_p<=2.702703e-05 | mega$Moth_EducationalAttainment_p<=2.702703e-05 | mega$Moth_SocialDeprivation_p<=2.702703e-05 |
                          mega$Fath_Height_p<=2.702703e-05 | mega$Fath_BMI_p<=2.702703e-05 | mega$Fath_Smoking_SmokingInitiation_p<=2.702703e-05 |
                          mega$Fath_Schizophrenia_p<=2.702703e-05 | mega$Fath_EducationalAttainment_p<=2.702703e-05 | mega$Fath_SocialDeprivation_p<=2.702703e-05),]

# read in mQTLs
mQTLs <- fread("/data/fhagenbeek/MethylationPGS/data/mQTLs.csv",data.table = FALSE)
# add CHR:BP as markername
mQTLs$MarkerName <- gsub("chr","",gsub(":INDEL","", gsub(":SNP","",mQTLs$name)))

# read in GWAS summary statistics 
SCZ <- fread("/data/fhagenbeek/MethylationPGS/data/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv",data.table = FALSE) 
Height <- fread("/data/fhagenbeek/MethylationPGS/data/Meta-analysis_Wood_et_al+UKBiobank_2018.txt",data.table = FALSE)
EA <- fread("/data/fhagenbeek/MethylationPGS/data/GWAS_EA_excl23andMe.txt",data.table = FALSE) 

# only retain genome-wide significant SNPs
SCZ <- SCZ[which(SCZ$PVAL<=5E-08),]
Height <- Height[which(Height$P<=5E-08),]
EA <- EA[which(EA$Pval<=5E-08),]

# add CHR:BP as MarkerName
SCZ$MarkerName <- paste(SCZ$CHROM,SCZ$POS,sep=":")
Height$MarkerName <- paste(Height$CHR,Height$POS,sep=":")
colnames(EA)[1] <- "rsid"
EA$MarkerName <- paste(EA$CHR,EA$POS,sep = ":")

# read in file with correlations between buccal and brain data
BBdat <- fread("/data/fhagenbeek/MethylationPGS/data/sommerer_et_al.results_summary_complete.tsv",data.table = FALSE)


################################################################################
#
# Are any of the significant CpGs present in previous EWAS studies?
#
################################################################################

# list of SCZ PGS CpGs from Hannon et al. 2016 (DOI 10.1186/s13059-016-1041-x)
Hannon <- c("cg05110828","cg03879918","cg14452952","cg14583365","cg01373089","cg26935265",
            "cg01999476","cg07172007","cg11594131","cg21852757","cg01651886","cg12414301",
            "cg25876057","cg04738464","cg17288324","cg12598198","cg26719024","cg01908272",
            "cg14927519","cg04846720","cg12555369","cg22715764","cg19529472","cg00329052",
            "cg04708391","cg24715815","cg17278466","cg14410503","cg13763120","cg03906681",
            "cg26254193","cg12064008","cg15070798","cg04195807","cg15928016","cg14632485",
            "cg07024458","cg24580784","cg05495029","cg18712348","cg14904733","cg19998012",
            "cg03994274","cg25920214","cg14599440","cg21821108","cg14266410","cg02679336",
            "cg02155405","cg06874767","cg08643855","cg13101865","cg13217184","cg12425283",
            "cg21525572","cg03074175","cg01235387","cg01450228","cg01989731","cg14060402",
            "cg17402245","cg06715755","cg11119301","cg07094993","cg13053563","cg19777067",
            "cg04384987","cg21484863","cg03108296","cg16171687","cg27093949","cg02912112",
            "cg18394881","cg22589169","cg01418910","cg21163714","cg11614536","cg26153715",
            "cg25700513","cg02493904","cg09846625","cg09405083","cg02631092","cg14571714",
            "cg20782653","cg24943198","cg19430423","cg03228312","cg07488338","cg09803164",
            "cg13782380","cg16739796","cg15094668","cg14021555","cg13177860","cg20725261",
            "cg11417675","cg26974111","cg22251684","cg04130549","cg14879574","cg09779405",
            "cg15941057","cg03451959","cg22721796","cg23743778","cg27174978","cg05335107",
            "cg12766287","cg11245719","cg17827670","cg27115821","cg02652569","cg21314839",
            "cg16692534","cg05580512","cg21859992","cg00504075","cg21334191","cg14494854",
            "cg10578782","cg19847436","cg07402939","cg22750364","cg07790870","cg20741169",
            "cg02739853","cg23966361","cg11961845","cg14523725","cg00499139","cg06521280",
            "cg04283218","cg12378950","cg09656541","cg21901307","cg18061532","cg14736087",
            "cg09130043","cg08940827","cg27193366","cg18958579","cg21953428","cg00827239",
            "cg26327548","cg24531714","cg27652795","cg03328526","cg06180389","cg21950166",
            "cg10500909","cg24114458","cg21697944","cg13390332","cg03483900","cg25250431")
# are any of the CpGs from Hannon found in our study?
sum(Hannon %in% mega_both$trait) # = 0, no overlap

# list of Smoking PGS CpGs from Yu et al. 2022 (https://doi.org/10.1080/15592294.2022.2088038) = nominally significant (p<0.05)
Yu_smoke <- c("cg10520740","cg24379915","cg24719910","cg18369034","cg01097768","cg13782301",
              "cg10255761","cg03227037","cg00177243","cg00071265","cg09084200","cg08331398",
              "cg00541303","cg19593285","cg04232972","cg19070894","cg23230929","cg01940273",
              "cg01995548","cg26132737","cg24158160","cg09570614","cg14030904","cg05575921",
              "cg08257009","cg04881228","cg21566642","cg12086028","cg06126421","cg26669717",
              "cg05951221","cg01765406","cg22132788","cg26788216","cg13127741","cg02787737",
              "cg00639656","cg00540464","cg00207731","cg10399789","cg03109660","cg23222488",
              "cg06644428","cg09935388","cg14796406","cg05875421","cg13001142","cg19918734",
              "cg12803068","cg24688690","cg23771366")
# are any of the smoking CpGs from Yu found in our study?
sum(Yu_smoke %in% mega_both$trait) # = 0, no overlap

# list of BMI PGS CpGs from Yu et al. 2022 (https://doi.org/10.1080/15592294.2022.2088038) = nominally significant (p<0.05)
Yu_BMI <- c("cg09613192","cg11202345","cg08972190","cg26257082","cg08305942",
            "cg26253134","cg01798813")
# are any of the BMI CpGs from Yu found in our study?
sum(Yu_BMI %in% mega_both$trait) # = 0, no overlap

# list of SCZ PGS CpGs from Viana et al. 2017 (10.1093/hmg/ddw373) 
Viana <- c("cg18847009","cg26893445","cg12595281","cg20640266","cg27150552","cg05209768",
           "cg07793808","cg10218777","cg01682070","cg11786558","cg26053083","cg01022840",
           "cg08478539","cg23788334","cg16904520")
# are any of the SCZ CpGs from Viana found in our study?
sum(Viana %in% mega_both$trait) # = 0, no overlap

# list of SCZ PGS CpGs from Kiltschewskij et al. 2023 (https://doi.org/10.1016/j.biopsych.2023.07.010) - suggestively significant threshold of 6.72x10-5
Kilt <- c("cg23296010","cg02786025","cg08977179","cg27153316","cg01236038","cg24253714",
          "cg15209053","cg26775289","cg05534464","cg25978835","cg16579770","cg24904935",
          "cg03040554","cg09349940","cg01978221","cg13435832","cg26369912","cg03341317",
          "cg18131836","cg19882427","cg14126918","cg17341955","cg27382744")
# are any of the SCZ CpGs from Viana found in our study?
sum(Kilt %in% mega_both$trait) # = 0, no overlap


################################################################################
#
# Are any of the identified mQTLs present in previous GWAS studies?
#
################################################################################

# overlap mQTLs and significant GWAS SNPs
mSCZ <- dplyr::inner_join(mQTLs[which(mQTLs$PGS=="Maternal Schizophrenia"),],SCZ,by="MarkerName") #0
oSCZ <- dplyr::inner_join(mQTLs[which(mQTLs$PGS=="Offspring Schizophrenia"),],SCZ,by="MarkerName") #0
oEA <- dplyr::inner_join(mQTLs[which(mQTLs$PGS=="Offspring Educational Attainment"),],EA,by="MarkerName") #0
pHeight <- dplyr::inner_join(mQTLs[which(mQTLs$PGS=="Paternal Height"),],Height,by="MarkerName") #56
oHeight <- dplyr::inner_join(mQTLs[which(mQTLs$PGS=="Offspring Height"),],Height,by="MarkerName") #17

# save overlap to txt file
write.table(pHeight, file=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                    "_Overlap_mQTLs_SNPs_PaternalHeight.txt"), row.names=F, col.names = T, sep="\t",quote = T) 
write.table(oHeight, file=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                    "_Overlap_mQTLs_SNPs_OffspringHeight.txt"), row.names=F, col.names = T, sep="\t",quote = T) 


# adjust data for plotting - both paternal and offpsring only include 1 CpG for which mQTLs were found
pHeight <- dplyr::full_join(mQTLs[which(mQTLs$PGS=="Paternal Height"),],Height,by="MarkerName") 
ph <- list(mQTL = unique(pHeight$MarkerName[which(!is.na(pHeight$cpg))]),
           SNP = unique(pHeight$MarkerName[which(!is.na(pHeight$SNP))]))
#
oHeight <- dplyr::full_join(mQTLs[which(mQTLs$PGS=="Offspring Height"),],Height,by="MarkerName") 
oh <- list(mQTL = unique(oHeight$MarkerName[which(!is.na(oHeight$cpg))]),
           SNP = unique(oHeight$MarkerName[which(!is.na(oHeight$SNP))]))

# create plots
paternalplot <- ggvenn(ph, fill_color = cbPalette[3:4],
                       stroke_size = 0.5, set_name_size = 4) +  ggtitle("(a) Paternal")
offspringplot <- ggvenn(oh, fill_color = cbPalette[3:4],
                        stroke_size = 0.5, set_name_size = 4)  +  ggtitle("(b) Offspring")
# combine plots
venns <- plot_grid(paternalplot,offspringplot,nrow = 2)
# save plots
ggsave(filename=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/",as.character(Sys.Date()),
                       "_GEEalltraits_results_FIG_mQTLSNP_Venn.pdf"),
       plot = venns, device = pdf, dpi = 300, width = 350, height = 250, units = "mm")


################################################################################
#
# Are any of the identified mQTLs present in the buccal-brain overlap data?
#
################################################################################

# dataset with significant associations
sign_cpgs <- dfc[which(!is.na(dfc$beta)),]
sign_cpgs$PGS <- paste(sign_cpgs$Parent,sign_cpgs$PGS,sep = " ")

# retain overlap between significant CpGs and buccal-brain correlation data (49/51 CpGs in buccal-brain data)
sign_BBdat <- dplyr::inner_join(sign_cpgs,BBdat, by = c("trait" = "CpG"))

# save overlap to txt file
write.table(sign_BBdat, file=paste0("/data/fhagenbeek/MethylationPGS/output/GEE_alltraits/", as.character(Sys.Date()),
                                 "_Overlap_signCpGs_buccalbraincorrelations.txt"), row.names=F, col.names = T, sep="\t",quote = T) 

## offspring SCZ ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
oSCZ.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Offspring Schizophrenia")])))
oSCZ.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Schizophrenia")]
oSCZ.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Schizophrenia")]
oSCZ.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Offspring Schizophrenia")]
# range, mean, and median correlations
round(range(oSCZ.Cor$value),2)
round(mean(oSCZ.Cor$value),2)
round(median(oSCZ.Cor$value),2)

## Maternal SCZ ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
mSCZ.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Maternal Schizophrenia")])))
mSCZ.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Maternal Schizophrenia")]
mSCZ.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Maternal Schizophrenia")]
mSCZ.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Maternal Schizophrenia")]
# range, mean, and median correlations
round(range(mSCZ.Cor$value),2)
round(mean(mSCZ.Cor$value),2)
round(median(mSCZ.Cor$value),2)

## offspring EA ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
oEA.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Offspring Educational Attainment")])))
oEA.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Educational Attainment")]
oEA.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Educational Attainment")]
oEA.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Offspring Educational Attainment")]
# range, mean, and median correlations
round(range(oEA.Cor$value),2)
round(mean(oEA.Cor$value),2)
round(median(oEA.Cor$value),2)

## Paternal EA ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
pEA.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Paternal Educational Attainment")])))
pEA.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Paternal Educational Attainment")]
pEA.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Paternal Educational Attainment")]
pEA.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Paternal Educational Attainment")]
# range, mean, and median correlations
round(range(pEA.Cor$value),2)
round(mean(pEA.Cor$value),2)
round(median(pEA.Cor$value),2)

## offspring Height ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
oHeight.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Offspring Height")])))
oHeight.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Height")]
oHeight.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Offspring Height")]
oHeight.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Offspring Height")]
# range, mean, and median correlations
round(range(oHeight.Cor$value),2)
round(mean(oHeight.Cor$value),2)
round(median(oHeight.Cor$value),2)

## Maternal Height ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
mHeight.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Maternal Height")])))
mHeight.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Maternal Height")]
mHeight.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Maternal Height")]
mHeight.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Maternal Height")]
# range, mean, and median correlations
round(range(mHeight.Cor$value),2)
round(mean(mHeight.Cor$value),2)
round(median(mHeight.Cor$value),2)

## Paternal Height ## 
# reorder significant CpGs with blood-brain correlation into "melted correlation matrix"
pHeight.Cor <- as.data.frame(matrix(NA,ncol=0,nrow = length(sign_BBdat$PGS[which(sign_BBdat$PGS=="Paternal Height")])))
pHeight.Cor$V1 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Paternal Height")]
pHeight.Cor$V2 <- sign_BBdat$trait[which(sign_BBdat$PGS=="Paternal Height")]
pHeight.Cor$value <- sign_BBdat$R[which(sign_BBdat$PGS=="Paternal Height")]
# range, mean, and median correlations
round(range(pHeight.Cor$value),2)
round(mean(pHeight.Cor$value),2)
round(median(pHeight.Cor$value),2)
