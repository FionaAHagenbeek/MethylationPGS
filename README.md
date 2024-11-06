# Intergenerational transmission of complex traits on the offspring methylome
Analysis scripts for the project "Intergenerational transmission of complex traits on the offspring methylome".

### Citation
Hagenbeek FA, et al. Intergenerational transmission of complex traits on the offspring methylome. [medRxiv](https://doi.org/10.1101/2024.04.15.24305824)

### Dependencies  
These scripts assume you have R v4.2.2 or higher installed. Required R packages: foreign, tidyverse/dplyr, foreach, doParallel, remotes, ewaff, RColorBrewer, plyr, ggplot2, reshape2, corrplot, gee, data.table, qqman, ggvenn, MASS, psych, lavaan, and cowplot.  

### Contact
Please contact Fiona Hagenbeek (fiona.hagenbeek@helsinki.fi) if you have any questions.

# 1. Data Preparation
The pipeline for DNA methylation–array analysis developed by the Biobank-based Integrative Omics Study (BIOS) consortium is available at [https://molepi.github.io/DNAmArray_workflow/](https://molepi.github.io/DNAmArray_workflow/).  
[01_MethylationPGS_DataPrep_ParticipantSelection_CreateDescriptives.R](01_DataPrep/01_MethylationPGS_DataPrep_ParticipantSelection_CreateDescriptives.R): combines the covariate and polygenic score files in each of the three datasets, performs participant exclusion, combines the three datasets, and creates the descriptive table(s).  
[02_MethylationPGS_DataPrep_CalculateMethylationVariability.R](01_DataPrep/02_MethylationPGS_DataPrep_CalculateMethylationVariability.R): combines the DNA methylation files of the three datasets and calculates the variability (SD) of the overlapping DNA methylation sites.  
[03_MethylationPGS_DataPrep_SelectTop10MostVariableMethylation.R](01_DataPrep/03_MethylationPGS_DataPrep_SelectTop10MostVariableMethylation.R): extracts the top 10% most variable DNA methylation sites and creates a figure that plots the variability by decile.  
[04_MethylationPGS_DataPrep_CreateFinalAnalysisDataset.R](01_DataPrep/04_MethylationPGS_DataPrep_CreateFinalAnalysisDataset.R): combines the covariates, polygenic scores, and DNA methylation sites for each of the three datasets into a single file for analysis.  
[05_MethylationPGS_DataPrep_MatrixSpectralDecomposition](01_DataPrep/05_MethylationPGS_DataPrep_MatrixSpectralDecomposition.R): script applies Matrix Spectral Decomposition (MSD) to the top 10% most variable DNA methylation sites in the NTR-ACTION Biomarker Study to determine the number of independent linear combinations of the DNA methylation sites to determine the Bonferonni significance threshold to apply in the analyses. 

# 2. Analyses
[01_MethylationPGS_CorrelationPGSs.R](02_Analyses/01_MethylationPGS_CorrelationPGSs.R): calculates the correlations between the polygenic scores and creates a plot where the upper triangle is a correlation heatmap and the lower triangle the corresponding correlations.  
[02_MethylationPGS_GEE_megaanalysis.R](02_Analyses/02_MethylationPGS_GEE_megaanalysis.R): runs the GEE models to assess the association of the parent and offspring polygenic scores for all complex traits simultaneously with each of the DNA methylation sites, while correcting for familial clustering.  
[03_MethylationPGS_GEE_ResultPlots.R](02_Analyses/03_MethylationPGS_GEE_ResultPlots.R): creates the annotated overview of all DNA methylation sites significantly associated with direct or indirect genetic effects and summarizes this in a figure and the Manhattan plots.  
[04_MethylationPGS_LookUps.R](02_Analyses/04_MethylationPGS_LookUps.R): compares significant CpGs associated with direct or indirect genetic effects to hits observed in previous EWAS studies, compares mQTLs associated with the significant CpGs associated with direct or indirect genetic effects to hits observed in previous GWAS studies, and looks the mQTLs up in a dataset describing the correlation between buccal and brain DNA methylation sites.

# 3. Simulation
[01_MethylationPGS_Simulation.R](03_Simulation/01_MethylationPGS_Simulation.R): simulates parental and offspring polygenic scores to compare the estimate for indirect genetic effects using the transmitted/non-transmitted model and the genetic parent-offspring model (see **Supplementary Note 3** in the preprint).
[02_MethylationPGS_MeanImputation_Simulations.R](03_Simulation/02_MethylationPGS_MeanImputation_Simulations.R): simulates parental and offspring polygenic scores where the percentage of missing parental polygenic scores can be set to compare the power given mean imputation of missing parental polygenic scores versus full information maximum likelihood (FIML).
