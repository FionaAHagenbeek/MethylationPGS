################################################################################
#
# project: Direct & Indirect Genetic Effects On Buccal DNA Methylation levels
#
# script: simulate phenotypes in the offspring and transmitted and
# non-transmitted polygenic scores for 100,000 parent-offspring trios based on
# the relation in the path diagram in Supplementary Figure 9 and compare the
# results of the transmitted/non-trasnmitted design to the genetic
# parent-offspring model for estimating genotype-environment covariance.
#
# Author: C.V. Dolan (c.v.dolan@vu.nl)
#
################################################################################

## simulate genetic and phenotypic data ##

# set allele frequecies (p and q) to .5 
p=.5; q=1-p

# number of simulated individuals 
N=100000

# simulate maternal and paternal transmitted and non-transmitted alleles 
ATM=sample(c(0,1),N,replace=T,prob=c(p,q)) # maternal transmitted
ANTM=sample(c(0,1),N,replace=T,prob=c(p,q)) # maternal non-transmitted
ATF=sample(c(0,1),N,replace=T,prob=c(p,q)) # paternal transmitted
ANTF=sample(c(0,1),N,replace=T,prob=c(p,q)) # paternal non-transmitted

# create simulated parent and offspring genotypes (sum of transmitted and non-transmitted)
GM=ATM+ANTM # maternal genotype
GF=ATF+ANTF # paternal genotype
GTkid=ATM+ATF # offspring genotype (transmitted from parents to offspring)
GNTkid=ANTM+ANTF # genotype based on alleles not transmitted from parents to offspring

# set direct genetic effect on the phenotype (parameter a) to .5 and the path
# from the parental phenotype to the offspring phenotype (path coefficient g,
# source of the genotype-shared environment covariance) to .2
a=.5
g=.2

# simulate the offspring phenotype
Ph=a*GTkid + rnorm(N,0,1) + g*a*GM + g*a*GF

## compare estimates of genotype-environment covariance as estimated with the transmitted/non-transmitted and genetic parent-offspring models ##

# run regression model for the transmitted/non-transmitted design
summary(lm(Ph~GTkid+GNTkid))

# run regression model for the genetic parent-offspring model
summary(lm(Ph~GTkid+I(GM+GF)))

# covariance (s1) and correlation (R1) between the simulated phenotype and transmitted/non-transmitted genotypes
S1=cov(cbind(Ph, GTkid, GNTkid))
R1=cov2cor(S1) 

# covariance (s2) and correlation (R2) between the simulated phenotype and parent-offspring genotypes
S2=cov(cbind(Ph, GTkid, GM+GF))
R2=cov2cor(S2)

# check the correlation between the genotypes
check=cor(cbind(GTkid, GM, GF))
