################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: simulate phenotypes in the offspring and offspring and parental
# polygenic scores for 2,000 parent-offspring trios (5,000 replications).
# Compare the effect with and without mean imputation for different percentages
# of parental polygenic scores.
#
# Author: C.V. Dolan (c.v.dolan@vu.nl)
#
################################################################################
#
# load required R-packages
# install.packages("MASS","psych","lavaan")
#
library(MASS)   # for the function mvrnorm - generate mv normal data
library(psych)  # function describe for summary stats
#
# if analyses in the simulation are to be performed with full information
# maximum likelihood (FIML; advised when not imputing missing polygenic scores)
# also require lavaan to do structural equation modelling (SEM)
#
dosem=FALSE # specify whether or not analyses require FIML
#           # the use of FIML provides full information estimates of the parameters
#           # and full information likelihood ratio tests of null hypotheses concerning the parameters
#           # these results may be of interest to gauge the effect of mean-imputation on the power 
#           # the comparison of interest being: power given mean-imputation vs power given FIML.
#
if (dosem) {
  library(lavaan)
}
#
# set number of replications and number of families
Nrep=5000  # number of replications 5000 is arbitrary but enough to get stable results of the simulation
Nfam=2000 # number of families 2000 is arbitrary but large enough to get stable results
#
# Create relatedness matrix with coefficients of relatedness and random mating
# (don't change)
#
RA=matrix(c(
  1,0,.5,
  0,1,.5,
  .5,.5,1),3,3)
#
# Set latent variance for latent genetic (A) and unique environmental (E)
# effects (don't change)
sa2 = 1 #  latent variance A 
se2 = 1 #  latent variance E
#
# proportion of explained variance polygenic score (change as you like)
R2PRS = .1 # i.e., R2PRS * sa2 is the PRS variance # the setting .1 implies that 10% of the A variance is accounted for by the PRS.
#
#
# phenotype variance decomposition in the offspring, given bParent = 0 (i.e., no
# Cultural transmission/genetic nurture)
# This is a simple AE twin model if bParent is zero
h2 = .5 # heritability, assuming no Cultural transmission
a_=sqrt(h2) # genetic effects
e_=sqrt(1-h2) # unique evironmental effects
#
# True effect of indirect genetic effects/cultural transmission/genetic nurture
# (change as you like)
#
bParent=c(.2) # c(0,.1,.2) #  arbritrary . effect size is given in the output called results
#
# Percentage of missingness in mother and father (change as you like)
# The percentages of interest were .08 (mother) and .19 (father). The values 0 are present for checking 
# the type I error rate (probability of rejecting bParent = 0, given in truth bParent = 0)
#
#                                                                                            setting of interest would be...
Perc_misM = c(0,.08,.19) # c(0,.10,.20,.30) # percentage missing data in parents arbitrary # Perc_misM = c(0,.08) 
Perc_misF = c(0,.08,.19) # c(0,.10,.20,.30) # percentage missing data in parents arbitrary # Perc_misF = c(0,.19)
#
#
# the design involves 2 factors bParent and Prec_mis: i.e., combinations of the
# true effect cultural transmission/genetic nurture/indirect genetic effects and
# percentages of missing maternal and paternal polygenic scores
#
ncell = length(bParent)*length(Perc_misM)*length(Perc_misF) # number of cells in factors design
keep=matrix(0, Nrep*ncell, 15) # to keep results, parameter estimates
keepSEM=matrix(0,Nrep*ncell, 5) # to keep results, parameter estimates (in case dosem==T)
#
# Set latent variance polygenic scores and residuals
#
SA_PRS=RA*R2PRS*sa2 # latent variance polygenic score = relatedness matrix by explained variance PRS by latent variance genetic effects
SA_res = RA*(1-R2PRS)*sa2 # residual variance = relatedness matrix by one minus the explained variance PRS by latent variance genetic effects
#
# counter
ii=0 
#
# Simulate the data and run the analyses
for (bP in bParent) { # loop across the different true effect sizes for the cultural transmission/genetic nurture/indirect genetic effects
  for (pmisM in Perc_misM) { # loop across the different percentages of missing maternal polygenic scores
    for (pmisF in Perc_misF) { # loop across the different percentages of missing paternal polygenic scores
      # 
      # replications
      for (irep in 1:Nrep) {   # loop across the number of repications
        #
        ii=ii+1 # set counter
        print(c(ii,ncell)) # print which counter is being run for how many combinations of the true effect cultural transmission/genetic nurture/indirect genetic effects and percentages of missing maternal and paternal polygenic scores
        # simulate data
        PRSdat = mvrnorm(Nfam, rep(0,3), Sigma= SA_PRS)   # simulate PRS scores
        Aresdat = mvrnorm(Nfam, rep(0,3), Sigma= SA_res)  # simulate residual A scores 
        A = PRSdat + Aresdat   # Simulate additive genetic scores
        E=rnorm(Nfam,0,1) # simulate unique environment scores
        bm=bf=bP # simulate effect cultural transmission/genetic nurture/indirect genetic effects to be equal in both parents
# build the regression model: 
       Ph = a_*A[,3] + e_*E +   bm*A[,1] + bf*A[,2] # simulate phenotype in offspring
#       
        # Create percentages of missing polygenic scores (pmis) in fathers and
        # mothers
#
        MismatM = matrix(runif(Nfam*1),Nfam,1) # matrix to contain the missing polygenic scores for mothers
        MismatF = matrix(runif(Nfam*1),Nfam,1) # matrix to contain the missing polygenic scores for fathers
        PRSdatmis = PRSdat # set up the polygenic score data (now equals the complete polygenic scores, will later be partially substituted with specified percentage missing/imputed polygenic scores)
        PRSdatmisMx = PRSdat # set up the polygenic score data (now equals the complete polygenic scores, will later be partially substituted with specified percentage missing/imputed polygenic scores)
        
        # mean impute percentage of missing parental polygenic scores
        PRSdatmis[MismatM[,1]<pmisM,1]=0 # impute mean for missing maternal polygenic scores (mean is 0 when polygenic scores are standardized to have a mean of 0 and sd of 1)
        PRSdatmis[MismatF[,1]<pmisF,2]=0 # impute mean for missing paternal polygenic scores (mean is 0 when polygenic scores are standardized to have a mean of 0 and sd of 1)
 #       
        # do not impute percentage of missing parental polygenic scores but keep
        # these as missing values (in case dosem==T)
        PRSdatmisMx[MismatM[,1]<pmisM,1]=NA # set prespecified percentage of maternal polygenic scores to missing 
        PRSdatmisMx[MismatF[,1]<pmisF,2]=NA # set prespecified percentage of paternal polygenic scores to missing
  #      
        # Create dataframe to retain the offspring phenotype (Ph) and polygenic
        # scores without mean imputation (in case dosem==T)
        PRSdatmisMx = as.data.frame( cbind(Ph,PRSdatmisMx)) # dataframe combining phenotype and polygenic scores
        colnames(PRSdatmisMx) = varnames =c('Ph','PRSm','PRSf','PRSc') # rename columns of new dataframe with phenotype and polygenic scores
  #      
        # specify the regressions for analyzing the association of the offspring
        # and parental polygenic scores with offspring phenotype when missing
        # parental polygenic scores have been mean imputed
        r1 = summary(lm(Ph ~ PRSdat[,3])) # r1: get the effect size of the association of the offspring polygenic score and offspring phenotype
        r2 =  summary(lm(Ph ~ PRSdat[,1]   +PRSdat[,2]   +PRSdat[,3])) # r2: get the effect size and parameters of the association of the offspring and parental polygenic scores (without missingness) and offspring phenotype
        r2m = summary(lm(Ph ~ PRSdatmis[,1]+PRSdatmis[,2]+PRSdatmis[,3])) #r3: get the effect size and parameters of the association of the offspring and parental polygenic scores (mean imputation missing parental polygenic scores) and offspring phenotype
  #      
        # Extract estimates and parameters
        keep[ii,1]=pmisM # percentage of missing maternal polygenic scores
        keep[ii,2]=pmisF # percentage of missing paternal polygenic scores
        keep[ii,3]=bP # true effect cultural transmission/genetic nurture/indirect genetic effect
        #
        keep[ii, 1+3] = r1$r.squared # proportion of explained variance model r1 (offspring phenotype ~ offspring polygenic score)
        keep[ii, 2+3] = r2$r.squared # proportion of explained variance model r2 (offspring phenotype ~ offspring + parental polygenic scores without missing data)
        keep[ii, 3+3] = r2$coefficients[2,1] # estimate maternal polygenic score from model r2
        keep[ii, 4+3] = r2$coefficients[3,1] # estimate paternal polygenic score from model r2
        keep[ii, 5+3] = r2$coefficients[4,1] # estimate offspring polygenic score from model r2
        keep[ii, 6+3] = r2m$coefficients[2,1] # estimate maternal polygenic score from model r3 (offspring phenotype ~ offspring + parental polygenic scores with mean imputation of missing parental polygenic scores)
        keep[ii, 7+3] = r2m$coefficients[3,1] # estimate paternal polygenic score from model r3
        keep[ii, 8+3] = r2m$coefficients[4,1] # estimate offspring polygenic score from model r3
        keep[ii, 9+3] =   as.numeric(r2$coefficients[2,4]<.05) # set p-value for the estimate of the maternal polygenic score from model r2 to zero if <0.05
        keep[ii, 10+3] =  as.numeric(r2$coefficients[3,4]<.05) # set p-value for the estimate of the paternal polygenic score from model r2 to zero if <0.05
        keep[ii, 11+3] = as.numeric(r2m$coefficients[2,4]<.05) # set p-value for the estimate of the maternal polygenic score from model r3 to zero if <0.05
        keep[ii, 12+3] = as.numeric(r2m$coefficients[3,4]<.05) # set p-value for the estimate of the paternal polygenic score from model r3 to zero if <0.05
   #     
        # combine offspring phenotype and offspring and parental polygenic
        # scores in a single dataframe
        PRSdat = as.data.frame( cbind(Ph,PRSdat)) # dataframe combining phenotype and polygenic scores
        colnames(PRSdat) = varnames =c('Ph','PRSm','PRSf','PRSc') # rename columns of new dataframe with phenotype and polygenic scores
 #       
        # use SEM (in lavaan) to run the regression analyses with FIML when
        # missing parental polygenic scores were not mean imputed
        if (dosem) {
          # model 1: offspring phenotype is predicted by offspring and parental polygenic scores 
          model1 <- '
  # specify regressions for model 1
    Ph_ =~ 1*Ph # latent offspring phenotype defined by the observed offspring phenotype with a fixed factor loading of 1
    PRSm_ =~ 1*PRSm # latent maternal polygenic score defined by the observed maternal polygenic score with a fixed factor loading of 1
    PRSf_ =~ 1*PRSf # latent paternal polygenic score defined by the observed paternal polygenic score with a fixed factor loading of 1
    PRSc_ =~ 1*PRSc # latent offspring polygenic score defined by the observed offspring polygenic score with a fixed facotr loading of 1
    Ph_ ~ PRSm_ + PRSf_ + PRSc_ # regression equation where the latent offspring phenotype is regressed on the latent offspring and parental polygenic scores
'
          fit1 <- sem(model1, missing="ML", data =  PRSdatmisMx) # run SEM model 1 with FIML
          #summary(fit1, standardized = FALSE)
          
          # model 2: offspring phenotype is predicted by offspring and parental polygenic scores where the maternal polygenic scores are constrained to zero
          model2 <- '
  # specify regressions for model 2
    Ph_ =~ 1*Ph # latent offspring phenotype defined by the observed offspring phenotype with a fixed factor loading of 1
    PRSm_ =~ 1*PRSm # latent maternal polygenic score defined by the observed maternal polygenic score with a fixed factor loading of 1
    PRSf_ =~ 1*PRSf # latent paternal polygenic score defined by the observed paternal polygenic score with a fixed factor loading of 1
    PRSc_ =~ 1*PRSc # latent offspring polygenic score defined by the observed offspring polygenic score with a fixed facotr loading of 1
    Ph_ ~ 0*PRSm_ + PRSf_ + PRSc_  # regression equation where the latent offspring phenotype is regressed on the latent offspring and parental polygenic scores and the maternal polygenic score is constrained to zero
'
          fit2 <- sem(model2, missing="ML", data =  PRSdatmisMx) # run SEM model 2 with FIML
          #summary(fit2, standardized = FALSE)
          
          # model 3: offspring phenotype is predicted by offspring and parental polygenic scores where the paternal polygenic scores are constrained to zero
          model3 <- '
  # specify regressions model 3
    Ph_ =~ 1*Ph # latent offspring phenotype defined by the observed offspring phenotype with a fixed factor loading of 1
    PRSm_ =~ 1*PRSm # latent maternal polygenic score defined by the observed maternal polygenic score with a fixed factor loading of 1
    PRSf_ =~ 1*PRSf # latent paternal polygenic score defined by the observed paternal polygenic score with a fixed factor loading of 1
    PRSc_ =~ 1*PRSc # latent offspring polygenic score defined by the observed offspring polygenic score with a fixed facotr loading of 1
    Ph_ ~ PRSm_ + 0*PRSf_ + PRSc_ # regression equation where the latent offspring phenotype is regressed on the latent offspring and parental polygenic scores and the paternal polygenic score is constrained to zero
'
          fit3 <- sem(model3, missing="ML", data =  PRSdatmisMx) # run SEM model 3 with FIML
          #summary(fit3, standardized = FALSE)
 #         
          # extract estimates and parameters from the SEM model 1
          keepSEM[ii,1:3]=summary(fit1)$pe[5:7,5]
 #         
          # likelihood ratio tests to compare the model fits
          ptestm=anova(fit1, fit2)[2,8] # testing whether the maternal polygenic score contributes significantly to the offspring phenotype
          ptestf=anova(fit1, fit3)[2,8] # testing whether the paternal polygenic score contributes significantly to the offspring phenotype
 #         
          # add the p-values from the likelihood ratio tests to the parameters
          # of SEM model 1
          keepSEM[ii,4:5]=c(ptestm, ptestf)
        } # dosem
        #
      } # replication
    } } } # end design: factors bP pmism pmisf

# analyse the results across the different combinations of true cultural
# transmission/genetic nurture/indirect genetic effects and percentages of
# missing parental polygenic scores

# initialize result matrix 
result=matrix(0,ncell, 15) # create matrix to contain results, populate with 0
colnames(result) = c("pmisM","pmisF","bP","r2ch","r2all","gm","gf","gc","gm_m","gf_m","gc_m",'mhit','fhit','mhit_m','fhit_m') # name result columns
ii=0 # initialize counter

# obtain results in loop
for (bP in bParent) { # loop across the different true cultural transmission/genetic nurture/indirect genetic effects
  for (pmisM in Perc_misM) { # loop across the different percentages of missing maternal polygenic scores
    for (pmisF in Perc_misF) { # loop across the different percentages of missing paternal polygenic scores
      
      ii=ii+1 # increment index/counter
      sel=(keep[,1]==pmisM & keep[,2]==pmisF & keep[,3]==bP) # creates a logical vector to select rows from the 'keep' matrix that match the cells in the current design
      tmp=keep[sel,] # filter 'keep' matrix using the logical vector 'sel' and stores this in a temporary matrix 'tmp'
      for (i in 1:15) { # loop across the result columns
        result[ii,i] = mean(tmp[,i]) # calculate the mean of the results (i.e., across the replicates)
      }
    }}}	
#
if (dosem) {
  resultSEM=matrix(0,ncell, 5) # initialize result matrix and populate with 0
  for (i in 1:3) { # loop across the columns
    resultSEM[1,i] = mean(keepSEM[,i]) # calculate the mean of the results (i.e., across the replicates)
  }
  resultSEM[1,4]=sum(keepSEM[,4]<.05)/Nrep # proportion of p-values <0.05 across the replicates where maternal polygenic scores significantly contribute to the offspring phenotype
  resultSEM[1,5]=sum(keepSEM[,5]<.05)/Nrep	# proportion of p-values <0.05 across the replicates where paternal polygenic scores significantly contribute to the offspring phenotype
  resultSEM=round(resultSEM,3) # round the simulation results without mean imputation to 3 decimals
  print(resultSEM) # print the simulation results without mean imputation
}

# round the simulation results with mean imputation to 3 decimals
result=round(result,3)

# print the simulation results with mean imputation
result