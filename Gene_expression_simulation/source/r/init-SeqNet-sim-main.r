#############################################################################################################
# 
# This is a main init script for an RNA-seq simulation project;
#
# Program:  init-SeqNet-sim-main.r
# Version:  
# Author:   zshahjahan
# Purpose:  Initializes important global variables. This includes knitR setup.
# Input:    
# Output:  	N/A
#############################################################################################################


##
##The raw data was downloaded from this link: https://www.kaggle.com/datasets/lachmann12/human-liver-rnaseq-gene-expression-903-samples/data
##Manually decompressed and stored in the dtaDir location


## reset variables
rm(list=ls(all=TRUE))

#Directories section
projDir = file.path("C:/Users/zayed/OneDrive/Desktop/Research_v2.0/SeqNet_playthru")
srcDir  = file.path(projDir,"source")
workDir = file.path(srcDir,"r")
funDir  = file.path(workDir,"functions")
dtaDir  = file.path(projDir,"data")
## NOTE: subdirectories under the analysis directory should correspond to the
## modules from which the objects, tables or models are being created from
## however, the exception to this should be any reference genome, transcriptome
## or database that has been downloaded and processed and may be used by multiple
## programs
anaDir  = file.path(projDir,"analysis")
resDir  = file.path(projDir,"report")

tblDir  = file.path(resDir,"tab")
pltDir  = file.path(resDir,"fig")

#Working directory should be the r folder under source dir
setwd(workDir)
#Source functions from script containing functions in the function directory
source(file.path(funDir,"init-functions.r"))

#My system using MRAN instead of CRAN, this is causing dependency issues
#trying this code will apparently help. UPDATE: and it does :)
# local({r = getOption("repos")
# r["CRAN"] = "http://cran.r-project.org"
# options(repos=r)})
options(repos = c(CRAN = "https://cloud.r-project.org"))

#Vector of software packages required for this project; modify this whenever
#a new package is to be used.

#NOTE: tidyverse is a real pain to work with, just load dplyr instead and work
#with that
libs = c("dplyr",
         "SeqNet",
         "data.table",
         "edgeR",
         "GenomicFeatures",
         "rtracklayer",
         "stringr",
         "gridExtra",
         "cowplot",
         "lme4",
         "MASS",
         "rmutil",
         "glmnet",
         "ggplot2",
         "Cairo")

## load libraries using the shared R library
# Use sapply to check and install packages
#sapply(libs, check_and_install)

## load libraries using the shared R library
lapply(libs, require, character.only=T)


#Set a random seed
set.seed(80085)

#subset about 1000 random rows (genes) and 500 columns (samples)
#Earlier I arbitrarily chose the numbers below, but I am going pick the sample 
#size and gene size numbers to reflect the sample and gene numbers from the 
#paper by Goll et al. 2018

## Initialize analysis objects
#new idea, use common processes of the liver and group these together and simulated
#co-expression
#TG always 
liverPtwVect=c("bile",
                 "bilirubin",
                 "Cytochrome",
                 "conjugation")
#Specify study population covariate groups using the following analysis variables
covariateObjList = list(treatSite = c("Site A","Site B","Site C","Site D"),
                        treatType = c("TRT","PLAC"),
                        sexVar    = c("M","F"),
                        tranplType = c("HEST","SOTL"),
                        ethnVar    = c("ethn_1","ethn_2","ethn_3","ethn_4"))

#Construct analysis-by objects. These are the conditions used
#Flag to perform TMM-normalization
tmmFlag = c("no-TMM"=0,"TMM"=1)
#Flag to analyze reference gene expression or simulated gene expression
#Names of vectors correspond to the naming convention from SeqNet::gen_rnaseq()
refSim = c("reference"=0,"x"=1)

#Usually the metadata should be loaded when the main init script is called
#but this script is an exception

#Set up hyperparameters related to tuning the linear model to be used
alphaVals = seq(0,1,0.01)
nBoot=500
nFold=15
percBoot=0.25 #percentage of 

#Obtain standardization flags for the gene expression data
stdFlags = c('org','std');
stdLabls = c('Original Variables','Standardized Variables');



