################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: Select top 10% most variable (>SD) of the residualized imputed DNA
# methylation probes after removing outliers
#
# data: 1. SD of DNA Methylation probes (.RData) 2. annotation file EPIC data
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

################################################################################
#
# Load data
#
################################################################################

# SD DNA methylation betas ACTION by batch (.Rdata)
load("/home/fhagenbeek/data/MethylationPGS/data/intermediate/2023-12-08_EPIC_Buccal_SDprobes_alldats.RData")

#load annotation files 
epicanno <- get(load("/home/fhagenbeek/data/MethylationPGS/data/anno_epic_072017.RData"))

################################################################################
#
# Order SD of CpGs by largest to smallest SD
# write top 10% most variable probes to file
#
################################################################################

#order by SD
CpG <- CpG[with(CpG, order(ressd,decreasing = TRUE)),]

#how many CpGs are 10%?
ten <- quantile(1:nrow(CpG),probs = 0.1)

#select top 10% not imputed raw probes
CpG.top10 <- CpG[1:ten,]

#annotate top10% most variable probes
CpG.top10.anno <- merge(CpG.top10,epicanno[,c("name","chromosome","position","gene.symbol","gene.accession","gene.region","cpg.island.name",
                                              "relation.to.island")], by.x = "CpG", by.y = "name", all.x = T)

#write output to file
save(CpG.top10.anno, file=paste("/home/fhagenbeek/data/MethylationPGS/data/intermediate/", as.character(Sys.Date()), "_EPIC_Buccal_Top10VariableCpGs_alldats.RData",sep=""))
